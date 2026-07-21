# frozen_string_literal: true

module Hgvs2vcf
  # HGVSp -> candidate cDNA edits.
  #
  # An amino-acid substitution is inherently ambiguous at the nucleotide level:
  # the target amino acid may have several synonymous codons, each reachable from
  # the reference codon by a different set of nucleotide changes. We resolve as
  # far as the data allows (the reference codon is known from the transcript CDS)
  # and return every *minimal* candidate, flagging ambiguity for the caller.
  module Protein
    Candidate = Struct.new(:codon_pos, :ref_codon, :alt_codon, :changes, :n_changes, keyword_init: true) do
      # changes: array of [codon_index(0..2), ref_base, alt_base]
      def snv?
        n_changes == 1
      end
    end

    module_function

    # ref_codon: the reference codon (from transcript CDS), aa_pos: 1-based residue index.
    # Returns { candidates: [...], ambiguous: bool, reference_matches: bool }
    def back_translate(ref_codon, aa_ref, aa_pos, aa_alt)
      ref_codon = ref_codon.upcase
      reference_matches = GeneticCode.translate_codon(ref_codon) == aa_ref
      target_codons = aa_alt == "*" ? GeneticCode.codons_for("*") : GeneticCode.codons_for(aa_alt)

      scored = target_codons.map do |alt_codon|
        changes = (0..2).filter_map do |i|
          [i, ref_codon[i], alt_codon[i]] if ref_codon[i] != alt_codon[i]
        end
        Candidate.new(
          codon_pos: aa_pos, ref_codon: ref_codon, alt_codon: alt_codon,
          changes: changes, n_changes: changes.size
        )
      end.reject { |c| c.n_changes.zero? } # skip synonymous no-ops

      min = scored.map(&:n_changes).min
      minimal = scored.select { |c| c.n_changes == min }

      {
        candidates: minimal,
        ambiguous: minimal.size > 1,
        requires_mnv: min && min > 1,
        reference_matches: reference_matches
      }
    end

    # Convert a candidate to a cDNA-space edit relative to the codon's first CDS base.
    # cds_first = c-position of the codon's first base = 3*(aa_pos-1)+1
    def candidate_to_cdna_edit(candidate)
      cds_first = 3 * (candidate.codon_pos - 1) + 1
      if candidate.snv?
        i, ref, alt = candidate.changes.first
        { edit: :snv, c_pos: cds_first + i, ref: ref, alt: alt }
      else
        # Contiguous or spread changes -> represent as delins over the codon span.
        { edit: :delins, c_start: cds_first, c_end: cds_first + 2,
          ref: candidate.ref_codon, alt: candidate.alt_codon }
      end
    end
  end
end
