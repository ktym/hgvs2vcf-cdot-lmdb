# frozen_string_literal: true

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

      kind = desc[0] # p / c / n / g
      base, version = split_version(ref)

      result = {
        input: hgvs,
        reference: ref,
        reference_base: base,
        reference_version: version,
        reference_type: classify_reference(ref),
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

    # --- protein: p.Glu504Lys, p.E504K, p.Arg97* / p.Arg97Ter ------------------
    PROTEIN_SUB = /\Ap\.\(?
      (?<ref>[A-Z][a-z]{2}|[A-Z])
      (?<pos>\d+)
      (?<alt>[A-Z][a-z]{2}|[A-Z]|\*)
    \)?\z/x

    def self.parse_protein(desc)
      m = PROTEIN_SUB.match(desc)
      raise ParseError, "unsupported protein change #{desc.inspect} (only substitutions here)" unless m

      {
        edit: :protein_substitution,
        aa_ref: GeneticCode.aa1(m[:ref]),
        aa_pos: m[:pos].to_i,
        aa_alt: m[:alt] == "*" ? "*" : GeneticCode.aa1(m[:alt])
      }
    end

    # --- cdna: c.1510G>A, c.1510+5G>A, c.76_78del, c.76dup, c.76_77insACT ------
    POS = /(?<datum>\d+)(?<offset>[+-]\d+)?/ # CDS datum position + optional intronic offset
    SNV = /\Ac?\.?(?<p>#{POS})(?<ref>[ACGT])>(?<alt>[ACGT])\z/
    DEL = /\Ac?\.?(?<p1>#{POS})(?:_(?<p2>#{POS}))?del(?<seq>[ACGT]*)\z/
    DUP = /\Ac?\.?(?<p1>#{POS})(?:_(?<p2>#{POS}))?dup(?<seq>[ACGT]*)\z/
    INS = /\Ac?\.?(?<p1>#{POS})_(?<p2>#{POS})ins(?<seq>[ACGT]+)\z/
    DELINS = /\Ac?\.?(?<p1>#{POS})(?:_(?<p2>#{POS}))?delins(?<seq>[ACGT]+)\z/

    def self.parse_cdna(desc)
      body = desc.sub(/\A[cn]\./, "")
      if (m = SNV.match("c." + body))
        return { edit: :snv, pos: parse_pos(m[:p]), ref: m[:ref], alt: m[:alt] }
      end
      if (m = DELINS.match("c." + body))
        return { edit: :delins, pos_start: parse_pos(m[:p1]), pos_end: parse_pos(m[:p2] || m[:p1]), alt: m[:seq] }
      end
      if (m = INS.match("c." + body))
        return { edit: :ins, pos_start: parse_pos(m[:p1]), pos_end: parse_pos(m[:p2]), alt: m[:seq] }
      end
      if (m = DUP.match("c." + body))
        return { edit: :dup, pos_start: parse_pos(m[:p1]), pos_end: parse_pos(m[:p2] || m[:p1]), seq: m[:seq] }
      end
      if (m = DEL.match("c." + body))
        return { edit: :del, pos_start: parse_pos(m[:p1]), pos_end: parse_pos(m[:p2] || m[:p1]), seq: m[:seq] }
      end

      raise ParseError, "unsupported cDNA change #{desc.inspect}"
    end

    def self.parse_pos(str)
      m = /\A(?<datum>\d+)(?<offset>[+-]\d+)?\z/.match(str)
      { datum: m[:datum].to_i, offset: (m[:offset] || "0").to_i }
    end

    def self.parse_genomic(desc)
      body = desc.sub(/\Ag\./, "")
      if (m = /\A(?<pos>\d+)(?<ref>[ACGT])>(?<alt>[ACGT])\z/.match(body))
        return { edit: :snv, pos: { datum: m[:pos].to_i, offset: 0 }, ref: m[:ref], alt: m[:alt] }
      end

      raise ParseError, "unsupported genomic change #{desc.inspect}"
    end
  end
end
