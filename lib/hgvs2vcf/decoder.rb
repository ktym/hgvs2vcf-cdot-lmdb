# frozen_string_literal: true

require_relative "genetic_code"
require_relative "parser"
require_relative "protein"
require_relative "mapper"
require_relative "transcripts"
require_relative "sequence"

module Hgvs2vcf
  # Ties parsing, transcript resolution, coordinate projection and sequence
  # lookups into a single HGVS -> VCF conversion.
  class Decoder
    COMPLEMENT = { "A" => "T", "C" => "G", "G" => "C", "T" => "A", "N" => "N" }.freeze

    # HGVSp forms that are syntactically fine but carry no unique nucleotide change.
    UNRESOLVABLE_PROTEIN = {
      protein_frameshift: "frameshift: any of many indels can produce it; " \
                          "HGVSp does not identify the nucleotide change",
      protein_unknown: "p.? — consequence unknown",
      protein_no_protein: "p.0 — no protein produced; no nucleotide change is implied",
      protein_start_lost: "p.Met1? — initiation lost; the nucleotide change is not determined"
    }.freeze

    def initialize(transcripts:, sequence:)
      @tx = transcripts
      @seq = sequence
    end

    # Returns a hash with :input, :vcf (array of {chrom,pos,ref,alt}), :warnings,
    # :transcript, :ambiguous.
    def decode(hgvs)
      parsed = Parser.parse(hgvs)
      warnings = []

      if %i[refseq_protein ensembl_protein].include?(parsed[:reference_type])
        raise Parser::ParseError, "protein-accession HGVSp (NP_/ENSP) needs a protein->transcript map; " \
                                  "supply gene symbol or transcript for now"
      end
      if (reason = UNRESOLVABLE_PROTEIN[parsed[:edit]])
        raise Parser::ParseError, "#{parsed[:input]}: #{reason}"
      end

      tx_id = @tx.resolve_transcript_id(
        reference_type: parsed[:reference_type], reference: parsed[:reference],
        reference_base: parsed[:reference_base], reference_version: parsed[:reference_version]
      )
      tx = @tx.normalize(tx_id)
      mapper = Mapper.new(tx)

      vcf, ambiguous =
        case parsed[:edit]
        when :snv, :del, :dup, :ins, :delins, :inv
          [[cdna_to_vcf(parsed, mapper, warnings)], false]
        when :protein_substitution, :protein_extension, :protein_silent
          protein_substitution_to_vcf(parsed, tx, mapper, warnings)
        when :protein_del, :protein_dup, :protein_ins, :protein_delins
          protein_indel_to_vcf(parsed, tx, mapper, warnings)
        else
          raise Parser::ParseError, "edit type #{parsed[:edit]} is not supported"
        end

      { input: hgvs, transcript: tx_id, gene: tx[:gene], vcf: vcf,
        ambiguous: ambiguous, warnings: warnings }
    rescue Protein::AmbiguityError => e
      raise Parser::ParseError, e.message
    end

    private

    # --- coordinate helpers ---------------------------------------------------

    # `n.` numbers raw transcript positions; `c.` numbers the CDS with -N/*N for UTRs.
    def cds_type_for(parsed, pos)
      parsed[:coordinate_type] == "n" ? :cdna : pos[:cds_type]
    end

    def g_pos(parsed, mapper, pos)
      mapper.c_to_genomic(datum: pos[:datum], offset: pos[:offset], cds_type: cds_type_for(parsed, pos))
    end

    # A transcript range projected onto the genome: [contig, start0, end0] with
    # start0 <= end0 regardless of strand. On a spliced range the genomic interval
    # also covers the introns in between, which is what a genomic deletion of that
    # range actually removes — flagged in the warnings.
    def g_range(parsed, mapper, pos_start, pos_end, warnings)
      a = g_pos(parsed, mapper, pos_start)
      b = g_pos(parsed, mapper, pos_end)
      s0, e0 = [a[:pos0], b[:pos0]].minmax

      unless a[:intronic] || b[:intronic]
        tx_len = (mapper.cdna_coord(datum: pos_end[:datum], cds_type: cds_type_for(parsed, pos_end)) -
                  mapper.cdna_coord(datum: pos_start[:datum], cds_type: cds_type_for(parsed, pos_start))).abs + 1
        if e0 - s0 + 1 != tx_len
          warnings << "range spans a splice junction: the genomic interval also covers " \
                      "#{(e0 - s0 + 1) - tx_len} intronic bases"
        end
      end
      [a[:contig], s0, e0]
    end

    def revcomp(seq)
      seq.reverse.chars.map { |c| COMPLEMENT.fetch(c) }.join
    end

    # transcript-oriented -> genomic-oriented allele
    def to_genomic_strand(mapper, seq)
      mapper.plus? ? seq : revcomp(seq)
    end

    def vcf_record(contig, pos0, ref, alt)
      pos0, ref, alt = Normalize.left_align(@seq, contig, pos0, ref, alt)
      { chrom: contig, pos: pos0 + 1, ref: ref, alt: alt }
    end

    # transcript-oriented base at a cDNA-space position (reads genome, flips for minus strand)
    def transcript_base(mapper, c_datum, cds_type: :cds)
      g = mapper.c_to_genomic(datum: c_datum, cds_type: cds_type)
      b = @seq.base(g[:contig], g[:pos0])
      mapper.plus? ? b : COMPLEMENT.fetch(b)
    end

    # --- c. / n. edits --------------------------------------------------------

    def cdna_to_vcf(parsed, mapper, warnings)
      case parsed[:edit]
      when :snv    then snv_to_vcf(parsed, mapper, warnings)
      when :del    then del_to_vcf(parsed, mapper, warnings)
      when :dup    then dup_to_vcf(parsed, mapper, warnings)
      when :ins    then ins_to_vcf(parsed, mapper, warnings)
      when :delins then delins_to_vcf(parsed, mapper, warnings)
      when :inv    then inv_to_vcf(parsed, mapper, warnings)
      end
    end

    def snv_to_vcf(parsed, mapper, warnings)
      g = g_pos(parsed, mapper, parsed[:pos])
      g_ref = to_genomic_strand(mapper, parsed[:ref])
      g_alt = to_genomic_strand(mapper, parsed[:alt])

      observed = @seq.base(g[:contig], g[:pos0])
      if observed != g_ref
        warnings << "reference mismatch: HGVS implies #{g_ref} at #{g[:contig]}:#{g[:pos0] + 1} " \
                    "but genome has #{observed}"
      end

      vcf_record(g[:contig], g[:pos0], g_ref, g_alt)
    end

    def del_to_vcf(parsed, mapper, warnings)
      contig, s0, e0 = g_range(parsed, mapper, parsed[:pos_start], parsed[:pos_end], warnings)
      deleted = @seq.fetch(contig, s0, e0 + 1)
      check_stated_sequence(mapper, deleted, parsed[:seq], contig, s0, warnings)
      vcf_record(contig, s0, deleted, "")
    end

    def dup_to_vcf(parsed, mapper, warnings)
      contig, s0, e0 = g_range(parsed, mapper, parsed[:pos_start], parsed[:pos_end], warnings)
      duplicated = @seq.fetch(contig, s0, e0 + 1)
      check_stated_sequence(mapper, duplicated, parsed[:seq], contig, s0, warnings)
      # A duplicated genomic segment is inserted immediately 3' of itself on the
      # genome; left-alignment then moves it to the canonical VCF position.
      vcf_record(contig, e0 + 1, "", duplicated)
    end

    def ins_to_vcf(parsed, mapper, _warnings)
      a = g_pos(parsed, mapper, parsed[:pos_start])
      b = g_pos(parsed, mapper, parsed[:pos_end])
      unless (a[:pos0] - b[:pos0]).abs == 1
        raise Mapper::MapError,
              "insertion flanks are not adjacent on the genome " \
              "(#{a[:pos0] + 1} and #{b[:pos0] + 1}); the insertion point is ambiguous"
      end

      vcf_record(a[:contig], [a[:pos0], b[:pos0]].min + 1, "", to_genomic_strand(mapper, parsed[:alt]))
    end

    # An inversion is the reverse complement of the genomic segment itself, so it
    # needs no strand handling — only the segment's own orientation flips.
    def inv_to_vcf(parsed, mapper, warnings)
      contig, s0, e0 = g_range(parsed, mapper, parsed[:pos_start], parsed[:pos_end], warnings)
      inverted = @seq.fetch(contig, s0, e0 + 1)
      check_stated_sequence(mapper, inverted, parsed[:seq], contig, s0, warnings)
      vcf_record(contig, s0, inverted, revcomp(inverted))
    end

    def delins_to_vcf(parsed, mapper, warnings)
      contig, s0, e0 = g_range(parsed, mapper, parsed[:pos_start], parsed[:pos_end], warnings)
      vcf_record(contig, s0, @seq.fetch(contig, s0, e0 + 1), to_genomic_strand(mapper, parsed[:alt]))
    end

    # `c.76_78delACT` states the deleted bases in transcript orientation; verify them.
    def check_stated_sequence(mapper, genomic_seq, stated, contig, s0, warnings)
      return if stated.nil? || stated.empty?

      expected = mapper.plus? ? genomic_seq : revcomp(genomic_seq)
      return if stated == expected

      warnings << "reference mismatch: HGVS states #{stated} but the transcript has " \
                  "#{expected} at #{contig}:#{s0 + 1}"
    end

    # --- p. edits -------------------------------------------------------------

    def codon_at(mapper, aa_pos)
      first = 3 * (aa_pos - 1) + 1
      (0..2).map { |i| transcript_base(mapper, first + i) }.join
    end

    # Only the two residues HGVS names explicitly can be checked against the CDS.
    def check_reference_aas(mapper, parsed, warnings)
      [[parsed[:aa_start], parsed[:aa_ref_start]], [parsed[:aa_end], parsed[:aa_ref_end]]].uniq.each do |pos, aa|
        got = GeneticCode.translate_codon(codon_at(mapper, pos))
        next if got == aa

        warnings << "reference mismatch: codon #{pos} encodes #{got || '?'} but HGVS states #{aa}"
      end
    end

    def protein_substitution_to_vcf(parsed, tx, mapper, warnings)
      raise Parser::ParseError, "HGVSp requires a coding transcript" unless tx[:coding]

      if parsed[:edit] == :protein_extension
        unless parsed[:aa_alt]
          raise Parser::ParseError, "extension without a replacement residue does not identify a nucleotide change"
        end

        warnings << "extension: only the stop-codon change is emitted; the new reading frame is not verified"
      end
      if parsed[:edit] == :protein_silent && GeneticCode.codons_for(parsed[:aa_ref]).size == 1
        raise Parser::ParseError,
              "#{parsed[:aa_ref]} has a single codon; a silent change at residue #{parsed[:aa_pos]} is impossible"
      end

      aa_pos = parsed[:aa_pos]
      ref_codon = codon_at(mapper, aa_pos)
      bt = Protein.back_translate(ref_codon, parsed[:aa_ref], aa_pos, parsed[:aa_alt])
      warnings << "reference codon #{ref_codon} does not encode #{parsed[:aa_ref]}" unless bt[:reference_matches]
      if bt[:ambiguous]
        warnings << "amino-acid change is ambiguous at the nucleotide level (#{bt[:candidates].size} candidates)"
      end
      warnings << "amino-acid change requires a multi-nucleotide variant" if bt[:requires_mnv]

      vcf = bt[:candidates].map do |cand|
        edit = Protein.candidate_to_cdna_edit(cand)
        if edit[:edit] == :snv
          snv_to_vcf({ pos: { datum: edit[:c_pos], offset: 0, cds_type: :cds }, ref: edit[:ref], alt: edit[:alt] },
                     mapper, warnings)
        else
          # codon-level delins -> map both ends, emit MNV (left-align handles trimming)
          span = { pos_start: { datum: edit[:c_start], offset: 0, cds_type: :cds },
                   pos_end: { datum: edit[:c_end], offset: 0, cds_type: :cds },
                   alt: edit[:alt] }
          delins_to_vcf(span, mapper, warnings)
        end
      end
      [vcf, bt[:ambiguous]]
    end

    # p.Glu504del / dup are exact at the nucleotide level (whole codons). p.…ins /
    # delins are not: every inserted residue can come from several codons, so we
    # enumerate the alternatives and report the result as ambiguous.
    def protein_indel_to_vcf(parsed, tx, mapper, warnings)
      raise Parser::ParseError, "HGVSp requires a coding transcript" unless tx[:coding]

      aa_start = parsed[:aa_start]
      aa_end = parsed[:aa_end]
      raise Parser::ParseError, "residue range #{aa_start}_#{aa_end} is reversed" if aa_end < aa_start

      case parsed[:edit]
      when :protein_del, :protein_dup
        check_reference_aas(mapper, parsed, warnings)
        c_start, c_end = Protein.codon_span(aa_start, aa_end)
        span = { pos_start: { datum: c_start, offset: 0, cds_type: :cds },
                 pos_end: { datum: c_end, offset: 0, cds_type: :cds }, seq: nil }
        rec = parsed[:edit] == :protein_del ? del_to_vcf(span, mapper, warnings) : dup_to_vcf(span, mapper, warnings)
        [[rec], false]
      when :protein_ins
        unless aa_end == aa_start + 1
          raise Parser::ParseError, "insertion flanks must be consecutive residues (got #{aa_start}_#{aa_end})"
        end

        check_reference_aas(mapper, parsed, warnings)
        candidates = Protein.codon_candidates(parsed[:alt_aas])
        warnings << ambiguity_note(candidates, parsed[:alt_aas]) if candidates.size > 1
        c_end_of_first = 3 * aa_start
        vcf = candidates.map do |nt|
          ins_to_vcf({ pos_start: { datum: c_end_of_first, offset: 0, cds_type: :cds },
                       pos_end: { datum: c_end_of_first + 1, offset: 0, cds_type: :cds },
                       alt: nt }, mapper, warnings)
        end
        [vcf, candidates.size > 1]
      when :protein_delins
        check_reference_aas(mapper, parsed, warnings)
        candidates = Protein.codon_candidates(parsed[:alt_aas])
        warnings << ambiguity_note(candidates, parsed[:alt_aas]) if candidates.size > 1
        c_start, c_end = Protein.codon_span(aa_start, aa_end)
        vcf = candidates.map do |nt|
          delins_to_vcf({ pos_start: { datum: c_start, offset: 0, cds_type: :cds },
                          pos_end: { datum: c_end, offset: 0, cds_type: :cds },
                          alt: nt }, mapper, warnings)
        end
        [vcf, candidates.size > 1]
      end
    end

    def ambiguity_note(candidates, aas)
      "#{aas.join} can be encoded #{candidates.size} ways; all synonymous alternatives are returned"
    end
  end
end
