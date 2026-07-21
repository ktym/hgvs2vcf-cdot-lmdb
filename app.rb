# frozen_string_literal: true

# Local HGVS -> VCF decoding API, backed by the LMDB index (no JSON in RAM).
#   bundle exec ruby app.rb
# Env:
#   INDEX_DIR     LMDB index directory built by bin/build_index.rb
#   GENOME_FASTA  GRCh38 .fna with a samtools .fai
#   BUILD         genome build (default GRCh38)
require "sinatra"
require "json"
require "lmdb"
$LOAD_PATH.unshift File.expand_path("lib", __dir__)
require "hgvs2vcf/decoder"
require "hgvs2vcf/lmdb_index"

module App
  def self.decoder
    @decoder ||= begin
      # Read-only, shared mmap. Read txns are lock-free and zero-copy, so many
      # threads/processes can serve from the same index concurrently.
      # NOTE: we pass :nolock — LMDB opens lock.mdb for writing on env open,
      # which fails with EACCES when the index was built by another user (e.g.
      # root in docker) or lock.mdb is not writable. Skipping the lock is safe
      # for a static served index with no writer.
      # We intentionally omit :rdonly: opening named DBIs requires a write txn in
      # the lmdb gem, and a write txn on an :rdonly env raises Permission denied;
      # a read txn opens DBIs that become invalid after commit. :nolock is enough
      # for read-only serving — the app never writes.
      env = LMDB.new(ENV.fetch("INDEX_DIR"), maxdbs: 4, nolock: true)
      kv = txdb = nil
      env.transaction do
        kv = env.database("kv")
        txdb = env.database("tx")
      end
      idx = Hgvs2vcf::LmdbIndex.new(kv: kv, txdb: txdb, build: ENV.fetch("BUILD", "GRCh38"))
      seq = Hgvs2vcf::Sequence.new(ENV.fetch("GENOME_FASTA"))
      Hgvs2vcf::Decoder.new(transcripts: idx, sequence: seq)
    end
  end
end

configure do
  # Bind to all interfaces — inside a container the default 'localhost' is not
  # reachable through `-p host:container`, which shows up as a connection reset.
  set :bind, ENV.fetch("BIND", "0.0.0.0")
  set :port, ENV.fetch("PORT", 4567).to_i
  set :environment, ENV.fetch("APP_ENV", "production").to_sym
end

# GET /decode?hgvs=NM_000603.4:c.894T>G
get "/decode" do
  content_type :json
  hgvs = params["hgvs"] or halt(400, { error: "missing hgvs" }.to_json)
  begin
    App.decoder.decode(hgvs).to_json
  rescue Hgvs2vcf::Parser::ParseError, Hgvs2vcf::Transcripts::ResolveError, Hgvs2vcf::Mapper::MapError => e
    halt 422, { input: hgvs, error: e.message }.to_json
  end
end

# POST /decode  body: {"hgvs": ["...", "..."]}  (Content-Type: application/json)
post "/decode" do
  content_type :json
  request.body.rewind # Rack may have already read the input while parsing params
  raw = request.body.read.to_s
  halt(400, { error: "empty body; POST JSON like {\"hgvs\":[...]}" }.to_json) if raw.strip.empty?
  body = JSON.parse(raw)
  Array(body["hgvs"]).map do |h|
    App.decoder.decode(h)
  rescue StandardError => e
    { input: h, error: e.message }
  end.to_json
end

get("/health") { "ok" }
