# frozen_string_literal: true

require_relative "genetic_code"

module Hgvs2vcf
  # Parse the reference part and the variant description of an HGVS string.
  # Returns a Hash describing the reference and the edit; coordinate resolution
  # and sequence lookups happen later in the pipeline.
  module Parser
    ParseError = Class.new(StandardError)

    # Reference identifier classification.
    def self.classify_reference(ref)
      case ref
      when /\A(NM_|NR_|XM_|XR_)\d+(\.\d+)?\z/ then :refseq_transcript
      when /\A(NP_|XP_)\d+(\.\d+)?\z/         then :refseq_protein
      when /\AENST\d+(\.\d+)?\z/              then :ensembl_transcript
      when /\AENSP\d+(\.\d+)?\z/              then :ensembl_protein
      when /\A(NC_|NG_)\d+(\.\d+)?\z/         then :genomic
      else :gene_symbol
      end
    end

    def self.split_version(accession)
      base, ver = accession.split(".", 2)
      [base, ver&.to_i]
    end

    def self.parse(hgvs)
      raise ParseError, "empty HGVS" if hgvs.nil? || hgvs.strip.empty?

      ref, desc = hgvs.strip.split(":", 2)
      raise ParseError, "missing ':' in #{hgvs.inspect}" if desc.nil?

      # "NM_004006.2(DMD):c.1234A>G" — the parenthesised gene symbol is a comment
      # on the accession, not part of it.
      gene_hint = ref[/\(([^()]+)\)\z/, 1]
      ref = ref.sub(/\([^()]+\)\z/, "")

      kind = desc[0] # p / c / n / g
      base, version = split_version(ref)

      result = {
        input: hgvs,
        reference: ref,
        reference_base: base,
        reference_version: version,
        reference_type: classify_reference(ref),
        gene_hint: gene_hint,
        coordinate_type: kind
      }

      case kind
      when "p" then result.merge!(parse_protein(desc))
      when "c", "n" then result.merge!(parse_cdna(desc))
      when "g" then result.merge!(parse_genomic(desc))
      else raise ParseError, "unsupported coordinate type #{kind.inspect}"
      end
      result
    end

    # --- amino-acid tokens ----------------------------------------------------
    AA3 = "Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|" \
          "Pro|Ser|Thr|Trp|Tyr|Val|Ter|Sec|Xaa"
    AA1 = "[ARNDCQEGHILKMFPSTWYVUX*]"
    AA = "(?:#{AA3}|#{AA1})"
    AA_RE = /#{AA}/

    # --- protein --------------------------------------------------------------
    # Substitution / extension are resolvable to nucleotides; deletion, duplication,
    # insertion and deletion-insertion are resolvable up to synonymous-codon choice;
    # frameshift, silent and "unknown" forms are parsed but carry no unique variant.
    P_SUB    = /\A(?<ref>#{AA})(?<pos>\d+)(?<alt>#{AA})\z/
    P_SILENT = /\A(?<ref>#{AA})(?<pos>\d+)=\z/
    P_DELINS = /\A(?<ref1>#{AA})(?<pos1>\d+)(?:_(?<ref2>#{AA})(?<pos2>\d+))?delins(?<alt>(?:#{AA})+)\z/
    P_DEL    = /\A(?<ref1>#{AA})(?<pos1>\d+)(?:_(?<ref2>#{AA})(?<pos2>\d+))?del\z/
    P_DUP    = /\A(?<ref1>#{AA})(?<pos1>\d+)(?:_(?<ref2>#{AA})(?<pos2>\d+))?dup\z/
    P_INS    = /\A(?<ref1>#{AA})(?<pos1>\d+)_(?<ref2>#{AA})(?<pos2>\d+)ins(?<alt>(?:#{AA})+)\z/
    P_FS     = /\A(?<ref>#{AA})(?<pos>\d+)(?<alt>#{AA})?fs(?:(?:Ter|\*)(?<len>\d+|\?))?\z/
    P_EXT    = /\A(?<ref>#{AA})(?<pos>\d+)(?<alt>#{AA})?ext(?:Ter|\*)?(?<len>-?\d+|\?)?\z/
    P_MET1   = /\A(?:Met|M)1\?\z/

    def self.parse_protein(desc)
      body = desc.sub(/\Ap\./, "")
      # p.(Glu504Lys) — parentheses mark a *predicted* consequence.
      predicted = !(body =~ /\A\(.*\)\z/).nil?
      body = body[1..-2] if predicted
      base = { predicted: predicted }

      return base.merge(edit: :protein_unknown, reason: "p.? — consequence unknown") if body == "?"
      return base.merge(edit: :protein_no_protein, reason: "p.0 — no protein produced") if body == "0" || body == "0?"
      if P_MET1.match(body)
        return base.merge(edit: :protein_start_lost, reason: "p.Met1? — translation initiation lost")
      end

      if (m = P_SILENT.match(body))
        # p.Glu504= keeps the residue but not necessarily the codon: it stands for
        # any synonymous change, which back-translation enumerates like a substitution.
        aa = aa1(m[:ref])
        return base.merge(edit: :protein_silent, aa_ref: aa, aa_pos: m[:pos].to_i, aa_alt: aa)
      end
      if (m = P_DELINS.match(body))
        return base.merge(protein_range(m).merge(edit: :protein_delins, alt_aas: split_aas(m[:alt])))
      end
      if (m = P_DEL.match(body))
        return base.merge(protein_range(m).merge(edit: :protein_del))
      end
      if (m = P_DUP.match(body))
        return base.merge(protein_range(m).merge(edit: :protein_dup))
      end
      if (m = P_INS.match(body))
        return base.merge(protein_range(m).merge(edit: :protein_ins, alt_aas: split_aas(m[:alt])))
      end
      if (m = P_FS.match(body))
        return base.merge(edit: :protein_frameshift, aa_ref: aa1(m[:ref]), aa_pos: m[:pos].to_i,
                          aa_alt: m[:alt] && aa1(m[:alt]))
      end
      if (m = P_EXT.match(body))
        return base.merge(edit: :protein_extension, aa_ref: aa1(m[:ref]), aa_pos: m[:pos].to_i,
                          aa_alt: m[:alt] && aa1(m[:alt]))
      end
      if (m = P_SUB.match(body))
        return base.merge(edit: :protein_substitution, aa_ref: aa1(m[:ref]),
                          aa_pos: m[:pos].to_i, aa_alt: aa1(m[:alt]))
      end

      raise ParseError, "unsupported protein change #{desc.inspect}"
    end

    # Shared "AA pos [_ AA pos]" span for protein del/dup/ins/delins.
    def self.protein_range(m)
      { aa_ref_start: aa1(m[:ref1]), aa_start: m[:pos1].to_i,
        aa_ref_end: m[:ref2] ? aa1(m[:ref2]) : aa1(m[:ref1]),
        aa_end: m[:pos2] ? m[:pos2].to_i : m[:pos1].to_i }
    end

    def self.split_aas(str)
      str.scan(AA_RE).map { |a| aa1(a) }
    end

    def self.aa1(aa)
      GeneticCode.aa1(aa)
    rescue ArgumentError => e
      raise ParseError, e.message
    end

    # --- cdna: c.1510G>A, c.-30G>A, c.*30del, c.1510+5G>A, c.76_78del, ... ------
    # A position is [-|*]?<datum>[+|-<intron offset>]:
    #   c.100     CDS                c.-30    5' UTR        c.*30    3' UTR
    #   c.100+5   intron, 5' side    c.100-5  intron, 3' side
    POS = "[-*]?\\d+(?:[+-]\\d+)?"
    POS_RE = /\A(?<sign>[-*])?(?<datum>\d+)(?<offset>[+-]\d+)?\z/
    SNV    = /\A(?<p>#{POS})(?<ref>[ACGT])>(?<alt>[ACGT])\z/
    DELINS = /\A(?<p1>#{POS})(?:_(?<p2>#{POS}))?delins(?<alt>[ACGT]+)\z/
    INS    = /\A(?<p1>#{POS})_(?<p2>#{POS})ins(?<alt>[ACGT]+)\z/
    INV    = /\A(?<p1>#{POS})_(?<p2>#{POS})inv(?<seq>[ACGT]*)\z/
    DUP    = /\A(?<p1>#{POS})(?:_(?<p2>#{POS}))?dup(?<seq>[ACGT]*)\z/
    DEL    = /\A(?<p1>#{POS})(?:_(?<p2>#{POS}))?del(?<seq>[ACGT]*)\z/

    def self.parse_cdna(desc)
      body = desc.sub(/\A[cn]\./, "")
      if (m = SNV.match(body))
        return { edit: :snv, pos: parse_pos(m[:p]), ref: m[:ref], alt: m[:alt] }
      end
      if (m = DELINS.match(body))
        return { edit: :delins, pos_start: parse_pos(m[:p1]), pos_end: parse_pos(m[:p2] || m[:p1]), alt: m[:alt] }
      end
      if (m = INS.match(body))
        return { edit: :ins, pos_start: parse_pos(m[:p1]), pos_end: parse_pos(m[:p2]), alt: m[:alt] }
      end
      if (m = INV.match(body))
        return { edit: :inv, pos_start: parse_pos(m[:p1]), pos_end: parse_pos(m[:p2]), seq: m[:seq] }
      end
      if (m = DUP.match(body))
        return { edit: :dup, pos_start: parse_pos(m[:p1]), pos_end: parse_pos(m[:p2] || m[:p1]), seq: m[:seq] }
      end
      if (m = DEL.match(body))
        return { edit: :del, pos_start: parse_pos(m[:p1]), pos_end: parse_pos(m[:p2] || m[:p1]), seq: m[:seq] }
      end

      raise ParseError, "unsupported cDNA change #{desc.inspect}"
    end

    def self.parse_pos(str)
      m = POS_RE.match(str) or raise ParseError, "bad position #{str.inspect}"

      { datum: m[:datum].to_i,
        offset: (m[:offset] || "0").to_i,
        cds_type: case m[:sign]
                  when "-" then :utr5
                  when "*" then :utr3
                  else :cds
                  end }
    end

    def self.parse_genomic(desc)
      body = desc.sub(/\Ag\./, "")
      if (m = /\A(?<pos>\d+)(?<ref>[ACGT])>(?<alt>[ACGT])\z/.match(body))
        return { edit: :snv, pos: { datum: m[:pos].to_i, offset: 0, cds_type: :cds }, ref: m[:ref], alt: m[:alt] }
      end

      raise ParseError, "unsupported genomic change #{desc.inspect}"
    end
  end
end
