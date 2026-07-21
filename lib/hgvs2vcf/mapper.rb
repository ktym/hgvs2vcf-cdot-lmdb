# frozen_string_literal: true

module Hgvs2vcf
  # Projects transcript (cDNA / c.) coordinates onto the genome using a cdot
  # transcript record. Honors strand and the per-exon transcript<->genome
  # alignment "gap" (CIGAR), which is where RefSeq transcripts that don't align
  # perfectly to the reference genome are handled correctly.
  #
  # cdot exon tuple: [g_start(0-based), g_end(excl), exon_number, cdna_start(1-based),
  #                   cdna_end(1-based), gap|null]
  class Mapper
    MapError = Class.new(StandardError)

    # tx: a normalized transcript record (see Transcripts#normalize) with keys:
    #   :contig, :strand ("+"/"-"), :exons (as above), :cds_start_i (transcript
    #   coord of first CDS base, 1-based), :cds_end_i (transcript coord of last CDS base)
    def initialize(tx)
      @tx = tx
      @exons = tx[:exons].sort_by { |e| e[3] } # ascending cDNA order
    end

    # c. datum (+ intron offset) -> { contig:, pos0: (0-based), intronic: bool }
    # ref/alt base orientation is handled by the caller via #strand.
    def c_to_genomic(datum:, offset: 0, cds_type: :cds)
      cdna = cds_to_cdna(datum, cds_type)
      if offset.zero?
        { contig: @tx[:contig], pos0: cdna_to_genomic(cdna), intronic: false }
      else
        # Intronic: land on the exon boundary base, then step `offset` into the
        # intron in genomic space (sign flips on the minus strand).
        boundary0 = cdna_to_genomic(cdna)
        step = plus? ? offset : -offset
        { contig: @tx[:contig], pos0: boundary0 + step, intronic: true }
      end
    end

    def strand
      @tx[:strand]
    end

    def plus?
      @tx[:strand] == "+"
    end

    private

    # Map a c. datum to a 1-based transcript (cDNA) coordinate.
    # cds_type: :cds (c.N, N>=1), :utr5 (c.-N), :utr3 (c.*N)
    def cds_to_cdna(datum, cds_type)
      case cds_type
      when :cds  then @tx[:cds_start_i] + (datum - 1)
      when :utr5 then @tx[:cds_start_i] - datum
      when :utr3 then @tx[:cds_end_i] + datum
      else raise MapError, "unknown cds_type #{cds_type}"
      end
    end

    # 1-based cDNA coordinate -> 0-based genomic coordinate, walking the gap CIGAR.
    def cdna_to_genomic(cdna)
      exon = @exons.find { |e| cdna >= e[3] && cdna <= e[4] }
      raise MapError, "cDNA #{cdna} not within any exon (intronic/UTR out of range?)" unless exon

      g_start, g_end, _no, c_start, _c_end, gap = exon
      offset_in_exon = cdna - c_start # 0-based offset within exon, transcript 5'->3'

      genomic_consumed = walk_gap_to_genomic(gap, offset_in_exon)

      if plus?
        g_start + genomic_consumed
      else
        (g_end - 1) - genomic_consumed
      end
    end

    # Given a gap CIGAR and a transcript offset within the exon, return how many
    # genomic bases are consumed up to that transcript offset.
    #   M<n>: n aligned bases (consume transcript + genomic)
    #   I<n>: n bases inserted in transcript (consume transcript only)
    #   D<n>: n bases deleted from transcript (consume genomic only)
    def walk_gap_to_genomic(gap, target_tx_offset)
      return target_tx_offset if gap.nil? # 1:1 alignment

      tx = 0
      g = 0
      gap.split(/\s+/).each do |op|
        code = op[0]
        n = op[1..].to_i
        case code
        when "M"
          take = [n, target_tx_offset - tx].min
          take = 0 if take.negative?
          tx += take
          g += take
          return g if tx >= target_tx_offset
          tx += (n - take)
          g += (n - take)
        when "I" # transcript insertion: consumes transcript, no genomic
          take = [n, target_tx_offset - tx].min
          take = 0 if take.negative?
          tx += take
          return g if tx >= target_tx_offset
          tx += (n - take)
        when "D" # transcript deletion vs genome: consumes genomic only
          g += n
        else
          raise MapError, "unknown CIGAR op #{op.inspect}"
        end
      end
      g
    end
  end
end
