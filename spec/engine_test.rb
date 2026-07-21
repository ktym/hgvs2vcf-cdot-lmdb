# frozen_string_literal: true

# Pure-logic regression tests — no cdot/genome data required.
#   ruby spec/engine_test.rb
$LOAD_PATH.unshift File.expand_path("../lib", __dir__)
require "hgvs2vcf/genetic_code"
require "hgvs2vcf/parser"
require "hgvs2vcf/protein"
require "hgvs2vcf/mapper"
require "hgvs2vcf/sequence"
include Hgvs2vcf

$fails = 0
def check(desc, got, exp)
  ok = got == exp
  $fails += 1 unless ok
  puts "#{ok ? 'ok  ' : 'FAIL'} #{desc}#{ok ? '' : "  got=#{got.inspect} exp=#{exp.inspect}"}"
end

# --- parser --------------------------------------------------------------------
p1 = Parser.parse("ALDH2:p.Glu504Lys")
check("parse HGVSp type", [p1[:reference_type], p1[:aa_ref], p1[:aa_pos], p1[:aa_alt]], [:gene_symbol, "E", 504, "K"])
p2 = Parser.parse("NM_000603.4:c.894T>G")
check("parse RefSeq version", [p2[:reference_base], p2[:reference_version], p2[:edit]], ["NM_000603", 4, :snv])
p3 = Parser.parse("ENST00000297494.3:c.894T>G")
check("parse ENST", p3[:reference_type], :ensembl_transcript)
p4 = Parser.parse("NM_000690:c.1510G>A")
check("parse bare accession (no version)", [p4[:reference_version], p4[:pos][:datum]], [nil, 1510])

# --- protein back-translation --------------------------------------------------
r = Protein.back_translate("GAG", "E", 504, "K") # ALDH2
edits = r[:candidates].map { |c| Protein.candidate_to_cdna_edit(c) }
check("E504K unique", [r[:ambiguous], edits], [false, [{ edit: :snv, c_pos: 1510, ref: "G", alt: "A" }]])

r = Protein.back_translate("GAT", "D", 298, "E") # NOS3 — ambiguous
edits = r[:candidates].map { |c| Protein.candidate_to_cdna_edit(c) }.sort_by { |e| e[:alt] }
check("D298E ambiguous", [r[:ambiguous], edits],
      [true, [{ edit: :snv, c_pos: 894, ref: "T", alt: "A" }, { edit: :snv, c_pos: 894, ref: "T", alt: "G" }]])

# --- mapper --------------------------------------------------------------------
def tx(strand:, exons:, cds_start_i: 1, cds_end_i: 999)
  { contig: "NC", strand: strand, exons: exons, cds_start_i: cds_start_i, cds_end_i: cds_end_i }
end
mp = Mapper.new(tx(strand: "+", exons: [[1000, 1100, 0, 1, 100, nil]]))
check("plus c.1",   mp.c_to_genomic(datum: 1)[:pos0], 1000)
check("plus c.100", mp.c_to_genomic(datum: 100)[:pos0], 1099)
mm = Mapper.new(tx(strand: "-", exons: [[1000, 1100, 0, 1, 100, nil]]))
check("minus c.1",  mm.c_to_genomic(datum: 1)[:pos0], 1099)
mtwo = Mapper.new(tx(strand: "+", exons: [[1000, 1050, 0, 1, 50, nil], [2000, 2050, 1, 51, 100, nil]]))
check("splice c.51", mtwo.c_to_genomic(datum: 51)[:pos0], 2000)
check("intron c.50+3", mtwo.c_to_genomic(datum: 50, offset: 3)[:pos0], 1052)

# --- normalization -------------------------------------------------------------
require "tempfile"
fa = Tempfile.new(["hp", ".fna"]); fa.write(">c\nGGAAAATCGG\n"); fa.close
File.write("#{fa.path}.fai", "c\t10\t3\t10\t11\n")
s = Sequence.new(fa.path)
check("left-align homopolymer del", Normalize.left_align(s, "c", 5, "A", ""), [1, "GA", "G"])

puts($fails.zero? ? "\nALL PASS" : "\n#{$fails} FAILED")
exit($fails.zero? ? 0 : 1)
