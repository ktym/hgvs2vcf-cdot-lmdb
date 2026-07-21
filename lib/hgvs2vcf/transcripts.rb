# frozen_string_literal: true

require "json"
require "zlib"

module Hgvs2vcf
  # Loads a cdot transcript JSON (or a prebuilt SQLite index) and answers the
  # three resolution questions the decoder needs:
  #   * gene symbol -> MANE Select transcript   (default when no accession given)
  #   * accession without version -> latest version available
  #   * ENST / NM accession -> transcript record
  class Transcripts
    ResolveError = Class.new(StandardError)

    attr_reader :build

    def initialize(build: "GRCh38")
      @build = build
      @by_id = {}                 # "NM_000603.4" => record
      @versions = Hash.new { |h, k| h[k] = [] } # "NM_000603" => [3,4,...]
      @mane_by_symbol = {}        # "NOS3" => "NM_000603.4"
      @canonical_by_symbol = {}   # fallback: Ensembl canonical / longest
    end

    # Load a cdot *.json or *.json.gz file (RefSeq or Ensembl build file).
    def load_cdot(path)
      raw = path.end_with?(".gz") ? Zlib::GzipReader.open(path, &:read) : File.read(path)
      data = JSON.parse(raw)
      data.fetch("transcripts").each do |tx_id, tx|
        next unless tx.dig("genome_builds", @build)

        @by_id[tx_id] = tx
        base, ver = tx_id.split(".", 2)
        @versions[base] << ver.to_i if ver

        symbol = tx["gene_name"]
        next unless symbol

        tags = Array(tx["genome_builds"][@build]["tag"]) | Array(tx["tag"])
        if tags.any? { |t| t.to_s.include?("MANE_Select") || t.to_s == "MANE Select" }
          @mane_by_symbol[symbol] = tx_id
        end
        @canonical_by_symbol[symbol] ||= tx_id if tags.any? { |t| t.to_s.include?("canonical") }
      end
      self
    end

    # Resolve any of the accepted reference forms to a concrete transcript id.
    def resolve_transcript_id(reference_type:, reference:, reference_base:, reference_version:)
      case reference_type
      when :gene_symbol
        resolve_symbol(reference)
      when :refseq_transcript, :ensembl_transcript
        reference_version ? reference : latest_version(reference_base)
      else
        raise ResolveError, "cannot resolve reference type #{reference_type} to a transcript"
      end
    end

    def resolve_symbol(symbol)
      @mane_by_symbol[symbol] ||
        @canonical_by_symbol[symbol] ||
        raise(ResolveError, "no MANE Select / canonical transcript for gene symbol #{symbol.inspect}")
    end

    def latest_version(base)
      vers = @versions[base]
      raise ResolveError, "unknown accession #{base.inspect}" if vers.empty?

      "#{base}.#{vers.max}"
    end

    def record(tx_id)
      @by_id[tx_id] || raise(ResolveError, "transcript not loaded: #{tx_id.inspect}")
    end

    # Convert a cdot record into the flat shape Mapper expects.
    def normalize(tx_id)
      tx = record(tx_id)
      b = tx["genome_builds"].fetch(@build)
      start_codon = tx["start_codon"] # transcript coord of CDS start (0-based in cdot)
      stop_codon  = tx["stop_codon"]
      {
        id: tx_id,
        gene: tx["gene_name"],
        protein: tx["protein"],
        contig: b["contig"],
        strand: b["strand"],
        exons: b["exons"],
        # c.1 == first CDS base. cdot start_codon is 0-based transcript offset,
        # exon cdna coords are 1-based, so first CDS base (c.1) is at cdna start_codon+1.
        cds_start_i: start_codon ? start_codon + 1 : nil,
        cds_end_i: stop_codon,   # last CDS base (transcript coord); validate per data
        coding: !start_codon.nil?
      }
    end
  end
end
