# frozen_string_literal: true

# Score the decoder against bench/gold.jsonl (Ensembl Variant Recoder as oracle).
#
#   INDEX_DIR=... GENOME_FASTA=... ruby bench/run.rb [--verbose] [--category X]
#   ruby bench/run.rb --parse-only          # no index/genome needed
#
# --parse-only answers "which inputs does the parser even accept", which is the
# part of the gold set that can be checked without the 1 GB of reference data.
$LOAD_PATH.unshift File.expand_path("../lib", __dir__)
require "json"
require "hgvs2vcf/decoder"
require "hgvs2vcf/lmdb_index"

GOLD = File.expand_path("gold.jsonl", __dir__)

opts = {
  verbose: ARGV.include?("--verbose"),
  parse_only: ARGV.include?("--parse-only"),
  ensembl: ARGV.include?("--ensembl"),
  category: (i = ARGV.index("--category")) && ARGV[i + 1]
}

rows = File.readlines(GOLD).map { |l| JSON.parse(l) }
rows.select! { |r| r["category"] == opts[:category] } if opts[:category]

# --ensembl re-asks every question through the MANE Select *Ensembl* transcript.
# MANE guarantees the RefSeq and Ensembl members share sequence and CDS, so the
# same c./p. coordinates must land on the same genomic position — which makes the
# oracle's answer reusable as the expected result for the ENST path.
if opts[:ensembl]
  rows = rows.filter_map do |row|
    enst = row.dig("provenance", "mane_select", "ensembl_nuc") or next
    row = row.dup
    row["input"] = "#{enst}:#{row['input'].split(':', 2).last}"
    row["expected"] = row["expected"].merge("transcript" => enst)
    row
  end
end

def build_decoder
  require "lmdb"
  env = LMDB.new(ENV.fetch("INDEX_DIR"), maxdbs: 4, nolock: true)
  kv = txdb = nil
  env.transaction do
    kv = env.database("kv")
    txdb = env.database("tx")
  end
  Hgvs2vcf::Decoder.new(
    transcripts: Hgvs2vcf::LmdbIndex.new(kv: kv, txdb: txdb, build: ENV.fetch("BUILD", "GRCh38")),
    sequence: Hgvs2vcf::Sequence.new(ENV.fetch("GENOME_FASTA"))
  )
end

# The oracle returns rows in its own order; compare as sets.
def vcf_set(rows)
  rows.map { |v| [v["chrom"] || v[:chrom], v["pos"] || v[:pos], v["ref"] || v[:ref], v["alt"] || v[:alt]] }.sort
end

PAD = 60

# Two VCF rows describe the same variant when applying them to the reference
# yields the same sequence. This separates a real coordinate error from a mere
# normalization difference: VCF left-aligns, HGVS 3'-shifts, and the oracle does
# not always trim the anchor base, so equivalent rows routinely differ textually.
# The window must be identical for both sides, or two spellings of the same
# variant produce different strings purely because they start at different offsets.
def window(all_rows)
  [all_rows.map { |(_c, pos, _r, _a)| pos - 1 }.min - PAD,
   all_rows.map { |(_c, pos, ref, _a)| pos - 1 + ref.length }.max + PAD]
end

def haplotype(seq, (contig, pos, ref, alt), lo, hi)
  start0 = pos - 1
  seq.fetch(contig, lo, start0) + alt + seq.fetch(contig, start0 + ref.length, hi)
end

# Each row is applied on its own: an "ambiguous" result is a set of alternative
# variants, not a haplotype carrying all of them at once.
def haplotype_set(seq, rows, lo, hi)
  rows.map { |r| haplotype(seq, r, lo, hi) }.sort
end

def ref_matches_genome?(seq, rows)
  rows.all? { |(contig, pos, ref, _alt)| seq.fetch(contig, pos - 1, pos - 1 + ref.length) == ref }
end

def classify(row, got, seq)
  exp = row["expected"]
  want = vcf_set(exp["vcf"])
  have = vcf_set(got[:vcf])

  if want != have && want.map(&:first).uniq == have.map(&:first).uniq
    lo, hi = window(want + have)
    if haplotype_set(seq, want, lo, hi) == haplotype_set(seq, have, lo, hi)
      return [:equivalent, "#{have.inspect} (oracle #{want.inspect})"]
    end
  end

  if want == have
    # Identical coordinates via a different transcript id. Expected and harmless:
    # a bare gene symbol resolves to whichever member of the MANE pair the index
    # happens to hold, and the two are equivalent by definition.
    return [:pass, nil] if exp["transcript"].nil? || exp["transcript"] == got[:transcript]

    [:pass_other_transcript, "#{got[:transcript]} (oracle used #{exp['transcript']})"]
  elsif (want & have).any?
    [:partial, "#{(have - want).size} extra, #{(want - have).size} missing: #{have.inspect}"]
  else
    [:wrong_vcf, have.inspect]
  end
end

results = Hash.new { |h, k| h[k] = [] }

if opts[:parse_only]
  rows.each do |row|
    Hgvs2vcf::Parser.parse(row["input"])
    results[:parsed] << row
  rescue StandardError => e
    results[:unparsed] << [row, e.message]
  end
else
  dec = build_decoder
  seq = Hgvs2vcf::Sequence.new(ENV.fetch("GENOME_FASTA"))
  rows.each do |row|
    got = dec.decode(row["input"])
    got = { vcf: got[:vcf], transcript: got[:transcript], ambiguous: got[:ambiguous], warnings: got[:warnings] }
    status, detail = classify(row, got, seq)
    unless ref_matches_genome?(seq, vcf_set(got[:vcf]))
      status = :bad_ref_allele
      detail = vcf_set(got[:vcf]).inspect
    end
    results[status] << [row, detail, got]
  rescue StandardError => e
    results[:error] << [row, "#{e.class}: #{e.message}", nil]
  end
end

# --- report ---------------------------------------------------------------------
if opts[:parse_only]
  puts "parsed   #{results[:parsed].size}/#{rows.size}"
  puts "unparsed #{results[:unparsed].size}/#{rows.size}"
  results[:unparsed].group_by { |row, _| row["category"] }.sort.each do |cat, list|
    puts "\n  #{cat} (#{list.size})"
    list.each { |row, msg| puts "    #{row['input'].ljust(40)} #{msg}" }
  end
  exit 0
end

by_cat = rows.group_by { |r| r["category"] }
# Correct = same variant, whatever the textual normalization.
correct = %i[pass equivalent pass_other_transcript].flat_map { |s| results[s].map { |row, *| row } }
puts format("CORRECT %d/%d (%.0f%%)  — exact %d, equivalent-normalization %d, MANE-pair partner %d",
            correct.size, rows.size, 100.0 * correct.size / rows.size,
            results[:pass].size, results[:equivalent].size, results[:pass_other_transcript].size)
puts
puts "category                ok/total"
by_cat.sort.each do |cat, list|
  n = list.count { |r| correct.include?(r) }
  puts format("  %-22s %2d/%-2d %s", cat, n, list.size, n == list.size ? "" : "<-")
end

%i[error bad_ref_allele wrong_vcf partial].each do |status|
  next if results[status].empty?

  puts "\n#{status} (#{results[status].size})"
  results[status].each do |row, detail, got|
    puts "  #{row['input']}"
    puts "    expected #{vcf_set(row['expected']['vcf']).inspect} via #{row['expected']['transcript']}"
    puts "    got      #{detail}"
    puts "    warnings #{got[:warnings].inspect}" if opts[:verbose] && got
  end
end

exit(correct.size == rows.size ? 0 : 1)
