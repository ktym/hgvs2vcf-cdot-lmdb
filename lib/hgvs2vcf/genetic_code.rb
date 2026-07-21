# frozen_string_literal: true

module Hgvs2vcf
  # Standard genetic code (transl_table=1) plus helpers for the HGVSp path.
  module GeneticCode
    CODON_TO_AA = {
      "TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
      "CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
      "ATT" => "I", "ATC" => "I", "ATA" => "I", "ATG" => "M",
      "GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
      "TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
      "CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
      "ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
      "GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
      "TAT" => "Y", "TAC" => "Y", "TAA" => "*", "TAG" => "*",
      "CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
      "AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
      "GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
      "TGT" => "C", "TGC" => "C", "TGA" => "*", "TGG" => "W",
      "CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
      "AGT" => "S", "AGC" => "S", "AGA" => "R", "AGG" => "R",
      "GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G"
    }.freeze

    # amino acid (1-letter) -> [codons]
    AA_TO_CODONS = CODON_TO_AA.each_with_object(Hash.new { |h, k| h[k] = [] }) do |(codon, aa), acc|
      acc[aa] << codon
    end.freeze

    # three-letter <-> one-letter (Ter/* included)
    THREE_TO_ONE = {
      "Ala" => "A", "Arg" => "R", "Asn" => "N", "Asp" => "D", "Cys" => "C",
      "Gln" => "Q", "Glu" => "E", "Gly" => "G", "His" => "H", "Ile" => "I",
      "Leu" => "L", "Lys" => "K", "Met" => "M", "Phe" => "F", "Pro" => "P",
      "Ser" => "S", "Thr" => "T", "Trp" => "W", "Tyr" => "Y", "Val" => "V",
      "Ter" => "*", "Sec" => "U"
    }.freeze

    module_function

    def translate_codon(codon)
      CODON_TO_AA[codon.upcase]
    end

    def codons_for(aa_one_letter)
      AA_TO_CODONS[aa_one_letter] || []
    end

    def aa1(aa)
      return aa if aa.length == 1

      THREE_TO_ONE.fetch(aa) { raise ArgumentError, "unknown amino acid: #{aa}" }
    end
  end
end
