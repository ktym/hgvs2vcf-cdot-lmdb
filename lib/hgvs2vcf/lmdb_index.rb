# frozen_string_literal: true

require "json"
require_relative "packing"
require_relative "transcripts" # for Transcripts::ResolveError (shared error type)

module Hgvs2vcf
  # Read side of the LMDB index. Same interface as Transcripts (resolve + normalize)
  # so Decoder is storage-agnostic, but instead of holding the whole cdot JSON in
  # RAM it reads one small packed record per query straight from the mmap'd DBI.
  class LmdbIndex
    def initialize(kv:, txdb:, build: "GRCh38")
      @kv = kv
      @tx = txdb
      @build = build
      @contigs = nil
    end

    def resolve_transcript_id(reference_type:, reference:, reference_base:, reference_version:)
      case reference_type
      when :gene_symbol
        @kv["mane:#{reference}"] || @kv["canon:#{reference}"] ||
          raise(Transcripts::ResolveError, "no MANE Select / canonical transcript for #{reference.inspect}")
      when :refseq_transcript, :ensembl_transcript
        if reference_version
          reference
        else
          @kv["latest:#{reference_base}"] ||
            raise(Transcripts::ResolveError, "unknown accession #{reference_base.inspect}")
        end
      else
        raise Transcripts::ResolveError, "cannot resolve reference type #{reference_type}"
      end
    end

    def normalize(tx_id)
      bytes = @tx[tx_id] or raise Transcripts::ResolveError, "transcript not in index: #{tx_id.inspect}"
      r = Packing.unpack_transcript(bytes)
      {
        id: tx_id,
        gene: r[:gene],
        contig: contigs.fetch(r[:contig_id]),
        strand: r[:strand],
        exons: r[:exons],
        cds_start_i: r[:coding] ? r[:cds_start_i] : nil,
        cds_end_i: r[:coding] ? r[:cds_end_i] : nil,
        coding: r[:coding]
      }
    end

    private

    def contigs
      @contigs ||= begin
        raw = @kv["meta:contigs"] or raise Transcripts::ResolveError, "index missing meta:contigs (rebuild?)"
        JSON.parse(raw).each_with_object({}) { |(id, name), h| h[id.to_i] = name }
      end
    end
  end
end
