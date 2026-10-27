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
      # Serving one build's index against another build's FASTA raises no error
      # per query — it returns wrong coordinates — so refuse it up front.
      indexed = @kv["meta:build"]
      if indexed && indexed != build
        raise ArgumentError, "index was built for #{indexed} but BUILD is #{build}"
      end
    end

    def resolve_transcript_id(reference_type:, reference:, reference_base:, reference_version:)
      case reference_type
      when :gene_symbol
        @kv["mane:#{reference}"] || refuse_missing_mane(reference) || @kv["canon:#{reference}"] ||
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

    # A MANE gene whose MANE transcript is not in this build must not fall back
    # to the canonical tag: on GRCh37 that is Ensembl 82-87's choice, often a
    # different transcript, and the answer would be wrong without any warning.
    def refuse_missing_mane(symbol)
      tx_id = @kv["mane-missing:#{symbol}"] or return nil
      raise Transcripts::ResolveError,
            "MANE Select for #{symbol.inspect} is #{tx_id}, which is not in the #{@build} index; " \
            "give a transcript accession instead"
    end

    def contigs
      @contigs ||= begin
        raw = @kv["meta:contigs"] or raise Transcripts::ResolveError, "index missing meta:contigs (rebuild?)"
        JSON.parse(raw).each_with_object({}) { |(id, name), h| h[id.to_i] = name }
      end
    end
  end
end
