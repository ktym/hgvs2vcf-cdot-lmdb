# frozen_string_literal: true

module Hgvs2vcf
  # Random access to a reference FASTA via its .fai index (samtools faidx).
  # For production use a bgzipped .fa.gz + .fai + .gzi and mmap; this plain-FASTA
  # reader is enough to build/verify VCF alleles and to left-align indels.
  class Sequence
    def initialize(fasta_path, fai_path = nil)
      @path = fasta_path
      @fai = load_fai(fai_path || "#{fasta_path}.fai")
      @io = File.open(fasta_path, "rb")
    end

    # 0-based, half-open [start0, end0) -> upstream bases (uppercase).
    def fetch(contig, start0, end0)
      idx = @fai.fetch(contig) { raise "contig #{contig.inspect} not in .fai" }
      len = end0 - start0
      return "" if len <= 0

      line_bases = idx[:line_bases]
      line_width = idx[:line_width]
      # byte offset of the first requested base
      first = idx[:offset] + (start0 / line_bases) * line_width + (start0 % line_bases)
      # read enough bytes including newlines
      n_lines = (start0 % line_bases + len) / line_bases
      to_read = len + n_lines * (line_width - line_bases) + line_width
      @io.seek(first)
      @io.read(to_read).delete("\n\r")[0, len].upcase
    end

    def base(contig, pos0)
      fetch(contig, pos0, pos0 + 1)
    end

    private

    def load_fai(path)
      raise "missing FASTA index #{path}; run: samtools faidx <fasta>" unless File.exist?(path)

      File.foreach(path).each_with_object({}) do |line, h|
        name, length, offset, line_bases, line_width = line.chomp.split("\t")
        h[name] = {
          length: length.to_i, offset: offset.to_i,
          line_bases: line_bases.to_i, line_width: line_width.to_i
        }
      end
    end
  end

  # VCF normalization: parsimony + left-alignment against the reference
  # (equivalent to `bcftools norm`). ref/alt are the raw alleles at pos0 (0-based).
  module Normalize
    module_function

    # Returns [pos0, ref, alt] left-aligned and trimmed.
    def left_align(seq, contig, pos0, ref, alt)
      ref = ref.dup.upcase
      alt = alt.dup.upcase

      # 1. trim common suffix (keep >=1 base each)
      while ref.length > 1 && alt.length > 1 && ref[-1] == alt[-1]
        ref.chop!
        alt.chop!
      end
      # 2. shift left while alleles share no anchor and can roll left
      while ref.empty? || alt.empty? || ref[0] == alt[0]
        break if pos0.zero?

        if ref.empty? || alt.empty? || (ref[-1] == alt[-1] && ref.length.positive? && alt.length.positive?)
          left = seq.base(contig, pos0 - 1)
          ref = left + ref
          alt = left + alt
          pos0 -= 1
          # re-trim shared suffix created by the roll
          while ref.length > 1 && alt.length > 1 && ref[-1] == alt[-1]
            ref.chop!
            alt.chop!
          end
        else
          break
        end
      end
      # 3. trim common leading base for pure SNV/MNV where possible
      while ref.length > 1 && alt.length > 1 && ref[0] == alt[0]
        ref = ref[1..]
        alt = alt[1..]
        pos0 += 1
      end
      [pos0, ref, alt]
    end
  end
end
