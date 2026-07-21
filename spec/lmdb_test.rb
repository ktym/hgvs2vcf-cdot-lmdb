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
cdot_path = "/home/claude/cdot-0.2.26/tests/test_data/cdot.refseq.grch37.json"
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

puts($fails.zero? ? "\nALL PASS" : "\n#{$fails} FAILED")
exit($fails.zero? ? 0 : 1)
