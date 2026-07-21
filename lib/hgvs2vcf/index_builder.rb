# frozen_string_literal: true

require "json"
require "zlib"
require_relative "packing"

module Hgvs2vcf
  # Streams a cdot transcript JSON into two key/value databases (any object that
  # responds to []/[]= — an LMDB named DBI in production, a Hash in tests).
  #
  #   kv  : resolution + metadata
  #     mane:<SYMBOL>   -> tx_id     (MANE Select)
  #     canon:<SYMBOL>  -> tx_id     (canonical fallback)
  #     latest:<BASE>   -> tx_id     (highest version of an accession)
  #     meta:contigs    -> JSON {id => contig_name}
  #     meta:build      -> genome build
  #   txdb: tx_id -> packed coordinate record (see Packing)
  #
  # Call #add_cdot for each cdot file, then #finalize! once.
  class IndexBuilder
    def initialize(kv:, txdb:, build: "GRCh38")
      @kv = kv
      @tx = txdb
      @build = build
      @contigs = {}       # name => id
      @contig_next = 0
      @latest = {}        # base => [max_version_int, tx_id]
    end

    def add_cdot(path)
      raw = path.end_with?(".gz") ? Zlib::GzipReader.open(path, &:read) : File.read(path)
      data = JSON.parse(raw)
      data.fetch("transcripts").each { |tx_id, t| add_transcript(tx_id, t) }
      self
    end

    def finalize!
      @latest.each { |base, (_ver, tx_id)| @kv["latest:#{base}"] = tx_id }
      @kv["meta:contigs"] = JSON.dump(@contigs.map { |name, id| [id.to_s, name] }.to_h)
      @kv["meta:build"] = @build
      self
    end

    private

    def add_transcript(tx_id, t)
      b = t.dig("genome_builds", @build)
      return unless b

      exons = b.fetch("exons").map { |e| [e[0], e[1], e[3], e[4], e[5]] } # g_start,g_end,cdna_start,cdna_end,gap
      start_codon = t["start_codon"]
      coding = !start_codon.nil?

      @tx[tx_id] = Packing.pack_transcript(
        strand: b.fetch("strand"),
        contig_id: contig_id(b.fetch("contig")),
        coding: coding,
        cds_start_i: coding ? start_codon + 1 : 0, # c.1 == cdna(start_codon + 1)
        cds_end_i: t["stop_codon"] || 0,
        gene: t["gene_name"],
        exons: exons
      )

      base, ver = tx_id.split(".", 2)
      if ver
        v = ver.to_i
        cur = @latest[base]
        @latest[base] = [v, tx_id] if cur.nil? || v > cur[0]
      end

      symbol = t["gene_name"]
      return unless symbol

      tags = Array(b["tag"]) | Array(t["tag"])
      if tags.any? { |x| x.to_s.include?("MANE_Select") || x.to_s == "MANE Select" }
        @kv["mane:#{symbol}"] = tx_id
      elsif tags.any? { |x| x.to_s.include?("canonical") }
        @kv["canon:#{symbol}"] = tx_id
      end
    end

    def contig_id(name)
      @contigs.fetch(name) do
        id = @contig_next
        @contig_next += 1
        @contigs[name] = id
        id
      end
    end
  end
end
