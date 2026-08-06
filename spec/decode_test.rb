# frozen_string_literal: true

# End-to-end decode tests over a synthetic transcript set and a tiny FASTA.
#   ruby spec/decode_test.rb
#
# The reference is a 16-base motif repeated to 300 bases, so any position is
# hand-checkable as MOTIF[pos0 % 16]. The motif has no repeated 1/2/3-mers, so
# small indels have a single left-aligned representation and the expected VCF
# records below are exact, not "whatever the normalizer happens to produce".
$LOAD_PATH.unshift File.expand_path("../lib", __dir__)
require "json"
require "tempfile"
require "hgvs2vcf/index_builder"
require "hgvs2vcf/lmdb_index"
require "hgvs2vcf/decoder"
include Hgvs2vcf

$fails = 0
def check(desc, got, exp)
  ok = got == exp
  $fails += 1 unless ok
  puts "#{ok ? 'ok  ' : 'FAIL'} #{desc}#{ok ? '' : "  got=#{got.inspect} exp=#{exp.inspect}"}"
end

def check_error(desc, pattern)
  yield
  $fails += 1
  puts "FAIL #{desc}  (no error raised)"
rescue StandardError => e
  ok = pattern.match?(e.message)
  $fails += 1 unless ok
  puts "#{ok ? 'ok  ' : 'FAIL'} #{desc}#{ok ? '' : "  message=#{e.message.inspect}"}"
end

MOTIF = "ACGTCAGTACTGCATG"
GENOME = (0...300).map { |i| MOTIF[i % 16] }.join

# Three transcripts on one contig:
#   TXP  plus strand, coding, two exons  (cDNA 1-30 @ 100-129, 31-60 @ 200-229)
#   TXM  minus strand, coding, two exons (cDNA 1-30 @ 229-200, 31-60 @ 129-100)
#   TXN  plus strand, non-coding, one exon
# CDS runs cDNA 6..50, so c.1 = cDNA 6, c.-1 = cDNA 5, c.*1 = cDNA 51.
CDOT = {
  "cdot_version" => "test", "genome_builds" => ["GRCh38"], "genes" => {},
  "transcripts" => {
    "NM_000001.1" => {
      "id" => "NM_000001.1", "gene_name" => "TXP", "tag" => ["MANE_Select"],
      "start_codon" => 5, "stop_codon" => 50,
      "genome_builds" => { "GRCh38" => {
        "contig" => "NC_test", "strand" => "+",
        "exons" => [[100, 130, 0, 1, 30, nil], [200, 230, 1, 31, 60, nil]]
      } }
    },
    "NM_000002.1" => {
      "id" => "NM_000002.1", "gene_name" => "TXM", "tag" => ["MANE_Select"],
      "start_codon" => 5, "stop_codon" => 50,
      "genome_builds" => { "GRCh38" => {
        "contig" => "NC_test", "strand" => "-",
        "exons" => [[200, 230, 0, 1, 30, nil], [100, 130, 1, 31, 60, nil]]
      } }
    },
    "NR_000003.1" => {
      "id" => "NR_000003.1", "gene_name" => "TXN", "tag" => ["canonical"],
      "genome_builds" => { "GRCh38" => {
        "contig" => "NC_test", "strand" => "+",
        "exons" => [[100, 130, 0, 1, 30, nil]]
      } }
    }
  }
}.freeze

jf = Tempfile.new(["cdot", ".json"]); jf.write(JSON.dump(CDOT)); jf.close
kv = {}; txdb = {}
IndexBuilder.new(kv: kv, txdb: txdb, build: "GRCh38").add_cdot(jf.path).finalize!

fa = Tempfile.new(["ref", ".fna"]); fa.write(">NC_test\n#{GENOME}\n"); fa.close
File.write("#{fa.path}.fai", "NC_test\t300\t9\t300\t301\n")

dec = Decoder.new(transcripts: LmdbIndex.new(kv: kv, txdb: txdb, build: "GRCh38"),
                  sequence: Sequence.new(fa.path))

def rec(pos, ref, alt)
  { chrom: "NC_test", pos: pos, ref: ref, alt: alt }
end

def one(dec, hgvs)
  r = dec.decode(hgvs)
  raise "expected exactly one VCF record, got #{r[:vcf].size}" unless r[:vcf].size == 1

  r[:vcf].first
end

# --- c. substitution, UTR and intron coordinates --------------------------------
# c.1 = cDNA 6 = 0-based 105 = MOTIF[9] = "C"
check("c. CDS SNV",      one(dec, "NM_000001.1:c.1C>T"), rec(106, "C", "T"))
# c.-5 = cDNA 1 = 100 = MOTIF[4] = "C"
check("c. 5'UTR SNV",    one(dec, "NM_000001.1:c.-5C>T"), rec(101, "C", "T"))
# c.*1 = cDNA 51 = 220 = MOTIF[12] = "C"
check("c. 3'UTR SNV",    one(dec, "NM_000001.1:c.*1C>A"), rec(221, "C", "A"))
# c.25 = cDNA 30 = 129 (last exonic base); +3 steps into the intron
check("c. intron +3",    one(dec, "NM_000001.1:c.25+3C>T"), rec(133, "C", "T"))
# c.26 = cDNA 31 = 200 (first base of exon 2); -3 steps back into the intron
check("c. intron -3",    one(dec, "NM_000001.1:c.26-3A>T"), rec(198, "A", "T"))
check("ref mismatch is warned",
      dec.decode("NM_000001.1:c.1A>T")[:warnings].first,
      "reference mismatch: HGVS implies A at NC_test:106 but genome has C")

# --- c. indels ------------------------------------------------------------------
# c.2_4 = cDNA 7..9 = 106..108 = "TGC". seq[105] == seq[108] == "C", so deleting
# "TGC" at 106 is the same event as deleting "CTG" at 105: left-aligned it anchors
# on seq[104] = "A". This is the HGVS-3'-shift vs VCF-left-align difference.
check("c. del range",    one(dec, "NM_000001.1:c.2_4del"), rec(105, "ACTG", "A"))
check("c. del stated seq", one(dec, "NM_000001.1:c.2_4delTGC"), rec(105, "ACTG", "A"))
check("c. del single",   one(dec, "NM_000001.1:c.2del"), rec(106, "CT", "C"))
check("c. dup",          one(dec, "NM_000001.1:c.2_4dup"), rec(105, "A", "ACTG"))
check("c. ins",          one(dec, "NM_000001.1:c.2_3insTT"), rec(106, "C", "CTT"))
check("c. delins",       one(dec, "NM_000001.1:c.2_4delinsGG"), rec(107, "TGC", "GG"))
check("c. inv",          one(dec, "NM_000001.1:c.2_4inv"), rec(107, "TGC", "GCA"))
check("stated del seq mismatch is warned",
      dec.decode("NM_000001.1:c.2_4delAAA")[:warnings].first,
      "reference mismatch: HGVS states AAA but the transcript has TGC at NC_test:107")
# `del3` counts the affected bases instead of spelling them: same assertion about
# the reference, so the same answer, and a disagreement is a warning not an error.
check("c. del stated count",  one(dec, "NM_000001.1:c.2_4del3"), rec(105, "ACTG", "A"))
check("c. dup stated count",  one(dec, "NM_000001.1:c.2_4dup3"), rec(105, "A", "ACTG"))
check("c. inv stated count",  one(dec, "NM_000001.1:c.2_4inv3"), rec(107, "TGC", "GCA"))
check("c. del single stated count", one(dec, "NM_000001.1:c.2del1"), rec(106, "CT", "C"))
check("stated del count mismatch is warned",
      dec.decode("NM_000001.1:c.2_4del9")[:warnings].first,
      "length mismatch: HGVS states 9 bases but the position range covers 3 at NC_test:107")
check("minus c. del stated count", one(dec, "NM_000002.1:c.2_4del3"), rec(221, "CATG", "C"))
# whole-intron deletion: c.25+1 = 130 .. c.26-1 = 199
check("c. intronic del",  one(dec, "NM_000001.1:c.25+1_26-1del")[:ref].length, 71) # 70 + anchor
# UTR ranges: c.-5_-3 = cDNA 1..3 = 100..102; c.*1_*3 = cDNA 51..53 = 220..222
check("c. 5'UTR del",    one(dec, "NM_000001.1:c.-5_-3del"), rec(100, "TCAG", "T"))
check("c. 3'UTR del",    one(dec, "NM_000001.1:c.*1_*3del"), rec(220, "GCAT", "G"))
check("c. del across the CDS start", one(dec, "NM_000001.1:c.-2_1del"), rec(103, "GTAC", "G"))

# --- minus strand ---------------------------------------------------------------
# c.1 = cDNA 6 = 224 = MOTIF[0] = "A"; transcript sees the complement, "T"
check("minus c. SNV",    one(dec, "NM_000002.1:c.1T>C"), rec(225, "A", "G"))
# c.2_4 = cDNA 7..9 = 223,222,221 -> genomic 221..223 = "ATG"
check("minus c. del",    one(dec, "NM_000002.1:c.2_4del"), rec(221, "CATG", "C"))
check("minus del stated seq (transcript orientation)",
      dec.decode("NM_000002.1:c.2_4delCAT")[:warnings], [])
# c.2 = 223, c.3 = 222; "TT" on the transcript is "AA" on the genome, inserted
# between 222 and 223 and anchored on seq[222] = "T".
check("minus c. ins is reverse-complemented",
      one(dec, "NM_000002.1:c.2_3insTT"), rec(223, "T", "TAA"))
# A duplicated genomic segment is strand-agnostic; only the anchor moves.
check("minus c. dup",    one(dec, "NM_000002.1:c.2_4dup"), rec(221, "C", "CATG"))
check("minus c. delins", one(dec, "NM_000002.1:c.2_4delinsAA"), rec(222, "ATG", "TT"))

# --- n. on a non-coding transcript ----------------------------------------------
check("n. SNV",          one(dec, "NR_000003.1:n.6C>T"), rec(106, "C", "T"))
check("n. del",          one(dec, "NR_000003.1:n.7_9del"), rec(105, "ACTG", "A")) # == c.2_4del
check_error("c. on a non-coding transcript is rejected", /non-coding/) do
  dec.decode("NR_000003.1:c.1C>T")
end

# --- coordinate validation ------------------------------------------------------
check_error("intron offset off an exon boundary", /last base of an exon/) do
  dec.decode("NM_000001.1:c.24+3G>A")
end
check_error("intron offset past the intron", /exceeds the 70 bp intron/) do
  dec.decode("NM_000001.1:c.25+71A>C")
end
check_error("non-adjacent insertion flanks", /not adjacent/) do
  dec.decode("NM_000001.1:c.2_4insTT")
end
check("splice-spanning del is warned",
      dec.decode("NM_000001.1:c.24_27del")[:warnings].first,
      "range spans a splice junction: the genomic interval also covers 70 intronic bases")

# --- p. -------------------------------------------------------------------------
# codon 1 = c.1_3 = 105..107 = "CTG" (Leu); codon 2 = c.4_6 = 108..110 = "CAT" (His)
check("p. substitution",  one(dec, "TXP:p.Leu1Met"), rec(106, "C", "A"))
check("p. del",           one(dec, "TXP:p.Leu1del"), rec(105, "ACTG", "A"))
check("p. dup",           one(dec, "TXP:p.His2dup"), rec(108, "G", "GCAT"))
check("p. delins (single codon)", one(dec, "TXP:p.Leu1delinsMet"), rec(106, "C", "A"))
ins = dec.decode("TXP:p.Leu1_His2insGly")
check("p. ins enumerates Gly codons", [ins[:ambiguous], ins[:vcf].size], [true, 4])
# Inserted between 107 and 108; seq[107] = "G" is the anchor, except for the GGG
# candidate, which rolls one base further left onto seq[106] = "T".
check("p. ins records", ins[:vcf].map { |v| [v[:pos], v[:alt]] }.sort,
      [[107, "TGGG"], [108, "GGGA"], [108, "GGGC"], [108, "GGGT"]])
check("p. reference mismatch is warned",
      dec.decode("TXP:p.Ala1del")[:warnings].first,
      "reference mismatch: codon 1 encodes L but HGVS states A")
check_error("p. frameshift is rejected with a reason", /frameshift/) do
  dec.decode("TXP:p.Leu1MetfsTer12")
end
# p.Leu1= — codon 1 is CTG; the minimal synonymous alternatives are CTT/CTC/CTA
# (c.3) and TTG (c.1). CTA needs two changes and is dropped.
silent = dec.decode("TXP:p.Leu1=")
check("p. silent enumerates synonymous codons",
      [silent[:ambiguous], silent[:vcf].map { |v| [v[:pos], v[:ref], v[:alt]] }.sort],
      [true, [[106, "C", "T"], [108, "G", "A"], [108, "G", "C"], [108, "G", "T"]]])
check_error("p. silent on a single-codon residue is rejected", /single codon/) do
  dec.decode("TXP:p.Trp1=")
end
check_error("p.? is rejected", /consequence unknown/) { dec.decode("TXP:p.?") }
check_error("p.Met1? is rejected", /initiation lost/) { dec.decode("TXP:p.Met1?") }
check_error("too many codon combinations", /cannot be pinned down/) do
  dec.decode("TXP:p.Leu1_His2insLeuLeuLeu") # 6^3 = 216 > 64
end

# --- accession forms ------------------------------------------------------------
check("gene symbol -> MANE Select", dec.decode("TXP:c.1C>T")[:transcript], "NM_000001.1")
check("accession without version -> latest", dec.decode("NM_000001:c.1C>T")[:transcript], "NM_000001.1")
check("parenthesised gene symbol is ignored",
      dec.decode("NM_000001.1(TXP):c.1C>T")[:transcript], "NM_000001.1")

puts($fails.zero? ? "\nALL PASS" : "\n#{$fails} FAILED")
exit($fails.zero? ? 0 : 1)
