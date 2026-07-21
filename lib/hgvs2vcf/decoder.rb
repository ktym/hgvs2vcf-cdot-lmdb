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

    def initialize(transcripts:, sequence:)
      @tx = transcripts
      @seq = sequence
    end

    # Returns a hash with :input, :vcf (array of {chrom,pos,ref,alt}), :warnings,
    # :transcript, :ambiguous.
    def decode(hgvs)
      parsed = Parser.parse(hgvs)
      warnings = []

      if parsed[:reference_type] == :refseq_protein || parsed[:reference_type] == :ensembl_protein
        raise Parser::ParseError, "protein-accession HGVSp (NP_/ENSP) needs a protein->transcript map; " \
                                  "supply gene symbol or transcript for now"
      end

      tx_id = @tx.resolve_transcript_id(
        reference_type: parsed[:reference_type], reference: parsed[:reference],
        reference_base: parsed[:reference_base], reference_version: parsed[:reference_version]
      )
      tx = @tx.normalize(tx_id)
      mapper = Mapper.new(tx)

      vcf, ambiguous =
        case parsed[:edit]
        when :snv          then [[snv_to_vcf(parsed, tx, mapper, warnings)], false]
        when :protein_substitution then protein_to_vcf(parsed, tx, mapper, warnings)
        else
          raise Parser::ParseError, "edit type #{parsed[:edit]} not yet implemented in this prototype"
        end

      { input: hgvs, transcript: tx_id, gene: tx[:gene], vcf: vcf,
        ambiguous: ambiguous, warnings: warnings }
    end

    private

    # transcript-oriented base at a c. datum (reads genome, flips for minus strand)
    def transcript_base(tx, mapper, c_datum)
      g = mapper.c_to_genomic(datum: c_datum)
      b = @seq.base(g[:contig], g[:pos0])
      mapper.plus? ? b : COMPLEMENT.fetch(b)
    end

    def snv_to_vcf(parsed, tx, mapper, warnings)
      pos = parsed[:pos]
      g = mapper.c_to_genomic(datum: pos[:datum], offset: pos[:offset])
      ref_t = parsed[:ref]
      alt_t = parsed[:alt]
      # transcript orientation -> genomic orientation
      g_ref = mapper.plus? ? ref_t : COMPLEMENT.fetch(ref_t)
      g_alt = mapper.plus? ? alt_t : COMPLEMENT.fetch(alt_t)

      observed = @seq.base(g[:contig], g[:pos0])
      warnings << "reference mismatch: HGVS implies #{g_ref} at #{g[:contig]}:#{g[:pos0] + 1} but genome has #{observed}" if observed != g_ref

      pos0, ref, alt = Normalize.left_align(@seq, g[:contig], g[:pos0], g_ref, g_alt)
      { chrom: g[:contig], pos: pos0 + 1, ref: ref, alt: alt }
    end

    def protein_to_vcf(parsed, tx, mapper, warnings)
      raise Parser::ParseError, "HGVSp requires a coding transcript" unless tx[:coding]

      aa_pos = parsed[:aa_pos]
      cds_first = 3 * (aa_pos - 1) + 1
      ref_codon = (0..2).map { |i| transcript_base(tx, mapper, cds_first + i) }.join
      bt = Protein.back_translate(ref_codon, parsed[:aa_ref], aa_pos, parsed[:aa_alt])
      warnings << "reference codon #{ref_codon} does not encode #{parsed[:aa_ref]}" unless bt[:reference_matches]
      warnings << "amino-acid change is ambiguous at the nucleotide level (#{bt[:candidates].size} candidates)" if bt[:ambiguous]
      warnings << "amino-acid change requires a multi-nucleotide variant" if bt[:requires_mnv]

      vcf = bt[:candidates].map do |cand|
        edit = Protein.candidate_to_cdna_edit(cand)
        if edit[:edit] == :snv
          snv_to_vcf({ pos: { datum: edit[:c_pos], offset: 0 }, ref: edit[:ref], alt: edit[:alt] },
                     tx, mapper, warnings)
        else
          # codon-level delins -> map both ends, emit MNV (left-align handles trimming)
          g = mapper.c_to_genomic(datum: edit[:c_start])
          g_ref = mapper.plus? ? edit[:ref] : edit[:ref].reverse.chars.map { |c| COMPLEMENT.fetch(c) }.join
          g_alt = mapper.plus? ? edit[:alt] : edit[:alt].reverse.chars.map { |c| COMPLEMENT.fetch(c) }.join
          start0 = mapper.plus? ? g[:pos0] : g[:pos0] - 2
          pos0, ref, alt = Normalize.left_align(@seq, g[:contig], start0, g_ref, g_alt)
          { chrom: g[:contig], pos: pos0 + 1, ref: ref, alt: alt }
        end
      end
      [vcf, bt[:ambiguous]]
    end
  end
end
