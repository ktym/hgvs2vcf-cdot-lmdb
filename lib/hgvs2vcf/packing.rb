# frozen_string_literal: true

require "stringio"

module Hgvs2vcf
  # Byte layout for a per-transcript coordinate record. This is the contract the
  # Ruby builder/reader and the future Rust reader both implement. All integers
  # little-endian; strings are length-prefixed UTF-8/ASCII, no null terminator.
  #
  #   magic        u8    format version (=1)
  #   strand       u8    '+' (0x2B) or '-' (0x2D)
  #   coding       u8    1 = coding (CDS present), 0 = non-coding
  #   contig_id    u16   index into the interned contig table (kv: meta:contigs)
  #   cds_start_i  u32   transcript coord (1-based) of first CDS base (c.1); 0 if non-coding
  #   cds_end_i    u32   transcript coord of last CDS base;                  0 if non-coding
  #   gene_len     u16
  #   gene         [gene_len] bytes
  #   n_exons      u16
  #   n_exons ×:
  #     g_start    u32   0-based genomic start (half-open)
  #     g_end      u32   0-based genomic end (exclusive)
  #     cdna_start u32   1-based transcript coord of exon start
  #     cdna_end   u32   1-based transcript coord of exon end
  #     gap_len    u16
  #     gap        [gap_len] bytes  (cDNA_match GAP CIGAR, usually empty)
  module Packing
    MAGIC = 1
    HEAD = "CCCS<L<L<S<"   # magic,strand,coding,contig_id,cds_start_i,cds_end_i,gene_len
    EXON = "L<L<L<L<S<"    # g_start,g_end,cdna_start,cdna_end,gap_len

    module_function

    def pack_transcript(strand:, contig_id:, coding:, cds_start_i:, cds_end_i:, gene:, exons:)
      gene_b = gene.to_s.b
      out = [MAGIC, strand.ord, coding ? 1 : 0, contig_id,
             coding ? cds_start_i : 0, coding ? cds_end_i : 0, gene_b.bytesize].pack(HEAD)
      out << gene_b
      out << [exons.size].pack("S<")
      exons.each do |g_start, g_end, cdna_start, cdna_end, gap|
        gap_b = (gap || "").b
        out << [g_start, g_end, cdna_start, cdna_end, gap_b.bytesize].pack(EXON)
        out << gap_b
      end
      out.b
    end

    def unpack_transcript(bytes)
      io = StringIO.new(bytes.b)
      # HEAD = CCC S< L< L< S< = 1+1+1+2+4+4+2 = 15 bytes
      magic, strand, coding, contig_id, cds_start_i, cds_end_i, gene_len =
        io.read(15).unpack(HEAD)
      raise "bad record magic #{magic}" unless magic == MAGIC

      gene = io.read(gene_len)
      (n,) = io.read(2).unpack("S<")
      exons = Array.new(n) do
        g_start, g_end, cdna_start, cdna_end, gap_len = io.read(18).unpack(EXON)
        gap = gap_len.zero? ? nil : io.read(gap_len)
        # reconstruct in the 6-tuple shape Mapper expects (exon_no placeholder at idx 2)
        [g_start, g_end, 0, cdna_start, cdna_end, gap]
      end
      {
        strand: strand.chr, coding: coding == 1, contig_id: contig_id,
        cds_start_i: cds_start_i, cds_end_i: cds_end_i, gene: gene, exons: exons
      }
    end
  end
end
