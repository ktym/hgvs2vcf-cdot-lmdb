#!/usr/bin/env ruby
# frozen_string_literal: true

# Build the LMDB index from cdot JSON(s). Run offline, on each data refresh.
#   ruby bin/build_index.rb --out data/index --build GRCh38 \
#        data/cdot-0.2.33.refseq.grch38.json.gz data/cdot-0.2.33.ensembl.grch38.json.gz
require "optparse"
require "lmdb"
$LOAD_PATH.unshift File.expand_path("../lib", __dir__)
require "hgvs2vcf/index_builder"

opts = { out: "data/index", build: "GRCh38", mapsize: 8 * 1024**3 }
OptionParser.new do |o|
  o.on("--out DIR") { |v| opts[:out] = v }
  o.on("--build B") { |v| opts[:build] = v }
  o.on("--mapsize BYTES", Integer) { |v| opts[:mapsize] = v }
end.parse!
paths = ARGV
abort "usage: build_index.rb --out DIR [--build B] <cdot.json[.gz> ...]" if paths.empty?

require "fileutils"
FileUtils.mkdir_p(opts[:out])
env = LMDB.new(opts[:out], maxdbs: 4, mapsize: opts[:mapsize])
kv = env.database("kv", create: true)
txdb = env.database("tx", create: true)

builder = Hgvs2vcf::IndexBuilder.new(kv: kv, txdb: txdb, build: opts[:build])
# One big write transaction: LMDB commits once, atomically.
env.transaction do
  paths.each do |p|
    warn "loading #{p} ..."
    builder.add_cdot(p)
  end
  builder.finalize!
end
env.close
warn "done -> #{opts[:out]} (kv + tx DBIs)"
