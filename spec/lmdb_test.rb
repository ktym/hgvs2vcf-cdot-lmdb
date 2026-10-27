# frozen_string_literal: true
# Integration test for the LMDB index path. Uses a Hash as a stand-in for an
# LMDB named DBI (both respond to []/[]=), so it runs without the lmdb gem.
$LOAD_PATH.unshift File.expand_path("../lib", __dir__)
require "json"
require "tempfile"
require "hgvs2vcf/index_builder"
require "hgvs2vcf/lmdb_index"
require "hgvs2vcf/mapper"
require "hgvs2vcf/decoder"
include Hgvs2vcf

$fails = 0
def check(desc, got, exp)
  ok = got == exp
  $fails += 1 unless ok
  puts "#{ok ? 'ok  ' : 'FAIL'} #{desc}#{ok ? '' : "  got=#{got.inspect} exp=#{exp.inspect}"}"
end

# ---- 1) real cdot test data: AOAH, minus strand, has a gap ---------------------
# Optional: point CDOT_TEST_JSON at cdot's tests/test_data/cdot.refseq.grch37.json.
cdot_path = ENV["CDOT_TEST_JSON"] || "/home/claude/cdot-0.2.26/tests/test_data/cdot.refseq.grch37.json"
if File.exist?(cdot_path)
src = JSON.parse(File.read(cdot_path))
kv = {}; txdb = {}
IndexBuilder.new(kv: kv, txdb: txdb, build: "GRCh37").add_cdot(cdot_path).finalize!

idx = LmdbIndex.new(kv: kv, txdb: txdb, build: "GRCh37")

# accession resolution: bare -> latest, versioned -> itself
check("latest NM_001637", idx.resolve_transcript_id(reference_type: :refseq_transcript,
      reference: "NM_001637", reference_base: "NM_001637", reference_version: nil), "NM_001637.3")
check("versioned passthrough", idx.resolve_transcript_id(reference_type: :refseq_transcript,
      reference: "NM_001637.3", reference_base: "NM_001637", reference_version: 3), "NM_001637.3")

tx = idx.normalize("NM_001637.3")
check("strand", tx[:strand], "-")
check("contig", tx[:contig], src["transcripts"]["NM_001637.3"]["genome_builds"]["GRCh37"]["contig"])
check("coding", tx[:coding], true)
check("exon count", tx[:exons].size, src["transcripts"]["NM_001637.3"]["genome_builds"]["GRCh37"]["exons"].size)

# End-to-end cross-check: c.1 must map to the genomic CDS start. On the minus
# strand that is the HIGHER genomic coordinate = cdot's cds_end (1-based).
mp = Mapper.new(tx)
pos0 = mp.c_to_genomic(datum: 1)[:pos0]
expected_cds_end = src["transcripts"]["NM_001637.3"]["genome_builds"]["GRCh37"]["cds_end"]
check("minus-strand c.1 == cds_end (1-based)", pos0 + 1, expected_cds_end)
else
  puts "skip  real cdot data (set CDOT_TEST_JSON to enable)"
end

# ---- 2) synthetic plus-strand transcript + tiny FASTA: full decode ------------
synthetic = {
  "cdot_version" => "test", "genome_builds" => ["GRCh38"], "genes" => {},
  "transcripts" => {
    "NM_999999.1" => {
      "id" => "NM_999999.1", "gene_name" => "TST",
      "start_codon" => 0, "stop_codon" => 30, "biotype" => ["protein_coding"],
      "tag" => ["MANE_Select"],
      "genome_builds" => { "GRCh38" => {
        "contig" => "NC_test", "strand" => "+", "cds_start" => 1000, "cds_end" => 1030,
        "exons" => [[1000, 1100, 0, 1, 100, nil]]
      } }
    }
  }
}
jf = Tempfile.new(["syn", ".json"]); jf.write(JSON.dump(synthetic)); jf.close
kv2 = {}; tx2 = {}
IndexBuilder.new(kv: kv2, txdb: tx2, build: "GRCh38").add_cdot(jf.path).finalize!
idx2 = LmdbIndex.new(kv: kv2, txdb: tx2, build: "GRCh38")

# FASTA: 1100 bases, 'G' at 0-based 1000, rest 'A'; single line for a simple .fai
seq = ("A" * 1000) + "G" + ("A" * 99)
fa = Tempfile.new(["syn", ".fna"]); fa.write(">NC_test\n#{seq}\n"); fa.close
File.write("#{fa.path}.fai", "NC_test\t1100\t9\t1100\t1101\n")
sequence = Sequence.new(fa.path)

dec = Decoder.new(transcripts: idx2, sequence: sequence)

r = dec.decode("NM_999999.1:c.1G>A")
check("decode via accession", r[:vcf], [{ chrom: "NC_test", pos: 1001, ref: "G", alt: "A" }])
check("no ref warning", r[:warnings], [])

r2 = dec.decode("TST:c.1G>A") # symbol -> MANE Select -> NM_999999.1
check("decode via symbol(MANE)", [r2[:transcript], r2[:vcf].first[:pos]], ["NM_999999.1", 1001])

# ---- 3) MANE from the NCBI summary (GRCh37 has no MANE tags in cdot) ---------
def tx37(gene, tag)
  { "gene_name" => gene, "start_codon" => 0, "stop_codon" => 30, "biotype" => ["protein_coding"],
    "genome_builds" => { "GRCh37" => {
      "contig" => "NC_test", "strand" => "+", "tag" => tag, "exons" => [[1000, 1100, 0, 1, 100, nil]]
    } } }
end
cdot37 = {
  "transcripts" => {
    "NM_100000.2" => tx37("OLDNAME", nil),                       # MANE, under an older symbol
    "ENST00000100000.1" => tx37("OLDNAME", "basic,Ensembl_canonical"),
    "NM_200000.1" => tx37("GONE", nil),                          # older version than MANE's .3
    "ENST00000200000.1" => tx37("GONE", "basic,Ensembl_canonical"),
    "ENST00000300000.1" => tx37("NOMANE", "basic,Ensembl_canonical")
  }
}
jf37 = Tempfile.new(["c37", ".json"]); jf37.write(JSON.dump(cdot37)); jf37.close
summary = Tempfile.new(["mane", ".txt"])
summary.write("#NCBI_GeneID\tsymbol\tRefSeq_nuc\tEnsembl_nuc\tMANE_status\n" \
              "1\tNEWNAME\tNM_100000.2\tENST00000999999.4\tMANE Select\n" \
              "2\tGONE\tNM_200000.3\tENST00000888888.2\tMANE Select\n" \
              "3\tNEWNAME\tNM_100001.1\tENST00000777777.1\tMANE Plus Clinical\n")
summary.close
kv3 = {}; tx3 = {}
b3 = IndexBuilder.new(kv: kv3, txdb: tx3, build: "GRCh37")
b3.add_cdot(jf37.path).add_mane(summary.path).finalize!
check("mane stats (Plus Clinical ignored)", b3.mane_stats, [1, 1])
idx3 = LmdbIndex.new(kv: kv3, txdb: tx3, build: "GRCh37")
sym = ->(s) { idx3.resolve_transcript_id(reference_type: :gene_symbol, reference: s, reference_base: s, reference_version: nil) }
check("summary symbol -> exact MANE version", sym.("NEWNAME"), "NM_100000.2")
gone = begin; sym.("GONE"); rescue Transcripts::ResolveError => e; e.message; end
check("MANE version absent -> 422, no canonical fallback", gone.include?("NM_200000.3"), true)
check("gene without MANE still uses canonical", sym.("NOMANE"), "ENST00000300000.1")

mismatch = begin
  LmdbIndex.new(kv: kv3, txdb: tx3, build: "GRCh38"); nil
rescue ArgumentError => e
  e.message
end
check("index/BUILD mismatch refused", mismatch, "index was built for GRCh37 but BUILD is GRCh38")

puts($fails.zero? ? "\nALL PASS" : "\n#{$fails} FAILED")
exit($fails.zero? ? 0 : 1)
