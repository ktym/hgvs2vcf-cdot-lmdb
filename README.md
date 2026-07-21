# hgvs2vcf

Offline HTTP API server that converts HGVS notation to VCF coordinates (`CHROM POS REF ALT`). At query time it does not call VEP, Mutalyzer, or NCBI; it reads pre-built transcript models (cdot) and the reference genome (GRCh38) from a local index.

## What it does

| Input | Example | Transcript resolution |
|-------|---------|----------------------|
| HGVSp (gene symbol) | `ALDH2:p.Glu504Lys` | symbol → **MANE Select** |
| HGVSc (gene symbol) | `ALDH2:c.1510G>A` | symbol → **MANE Select** |
| HGVSc (RefSeq) | `NM_000690:c.1510G>A` | version omitted → **latest** |
| HGVSc (RefSeq, versioned) | `NM_000603.4:c.894T>G` | exact |
| HGVSc (Ensembl) | `ENST00000297494.3:c.894T>G` | exact / latest |

`p.` / `c.` / `n.` substitutions are implemented. `del` / `dup` / `ins` / `delins` and intronic offsets are parsed and mapped; indels are left-aligned using `bcftools norm`-equivalent normalization.

Responses include `vcf` (coordinate rows), `ambiguous` (when one amino-acid change maps to several nucleotide changes), and `warnings` (e.g. reference-allele mismatch).

## Quick start

### 1. Prepare data

Source filenames carry versions that change between releases, so downloads are manual. Grab the current files from the sources below, place them in `data/`, and symlink to these fixed names:

```bash
mkdir -p data && cd data

# 1) cdot transcript models — latest from https://github.com/SACGF/cdot/releases
#    (e.g. data_v0.2.33: cdot-0.2.33.refseq.GRCh38.json.gz / .ensembl.GRCh38.json.gz)
ln -s cdot-0.2.33.refseq.GRCh38.json.gz   cdot.refseq.GRCh38.json.gz
ln -s cdot-0.2.33.ensembl.GRCh38.json.gz  cdot.ensembl.GRCh38.json.gz

# 2) MANE summary (auditing / provenance; not read by the running app)
#    latest from https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/
#    (e.g. MANE.GRCh38.v1.5.summary.txt.gz)
ln -s MANE.GRCh38.v1.5.summary.txt.gz     MANE.GRCh38.summary.txt.gz

# 3) reference genome — NCBI GCF_000001405.40 (RefSeq contig accessions == cdot 'contig')
#    download GCF_000001405.40_GRCh38.p14_genomic.fna.gz from:
#    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/
ln -s GCF_000001405.40_GRCh38.p14_genomic.fna.gz  GRCh38.fna.gz

cd ..
```

See [Data sources](#data-sources) for what each file is and how it is used.

### 2. Run with a container (Docker / Podman)

```bash
docker build -t hgvs2vcf .

# build the LMDB index (re-run on each data refresh)
docker run -it --rm --user $(id -u):$(id -g) -v $(pwd)/data:/data hgvs2vcf index

# serve (first run decompresses GRCh38.fna.gz and runs samtools faidx once, ~3GB)
docker run -it --rm --user $(id -u):$(id -g) -v $(pwd)/data:/data -p 4567:4567 hgvs2vcf
```

For Podman, substitute `podman` and add the `:z` SELinux relabel on the mount:

```bash
podman run -it --rm --userns=keep-id -v $(pwd)/data:/data:z -p 4567:4567 hgvs2vcf
```

`data/` holds the symlinks, the built `index/` (LMDB), and the prepared `GRCh38.fna` (plus `.fai`). The app binds to `0.0.0.0` inside the container (override with `BIND`).

### 3. Install without a container

Data preparation is the same as above. If you prefer to run on the host with Ruby and Bundler instead of a container, follow the steps below.

`bundle install` installs Ruby gems only. The `lmdb` gem is a **native extension** and needs the **LMDB C library on the OS** (development headers and link library) plus a C compiler. You also need `samtools` to index the reference FASTA. The container image (section 2) already includes these.

```bash
# e.g. Debian/Ubuntu
sudo apt install build-essential liblmdb-dev samtools

# e.g. macOS (Homebrew)
brew install lmdb samtools
```

```bash
bundle install
ruby bin/build_index.rb --out data/index --build GRCh38 \
     data/cdot.refseq.GRCh38.json.gz data/cdot.ensembl.GRCh38.json.gz
gunzip -kc data/GRCh38.fna.gz > data/GRCh38.fna && samtools faidx data/GRCh38.fna
INDEX_DIR=data/index GENOME_FASTA=data/GRCh38.fna bundle exec ruby app.rb
```

## API

| Method | Path | Description |
|--------|------|-------------|
| GET | `/health` | Liveness check (no data required) |
| GET | `/decode?hgvs=...` | Decode a single HGVS string |
| POST | `/decode` | Batch decode via JSON `{"hgvs":["...", ...]}` |

**Liveness**

```bash
curl localhost:4567/health
# -> ok
```

**Single decode** — `>` must be URL-encoded.

```bash
curl -G localhost:4567/decode --data-urlencode 'hgvs=NM_000603.4:c.894T>G'
```

```json
{"input":"NM_000603.4:c.894T>G","transcript":"NM_000603.4","gene":"NOS3",
 "vcf":[{"chrom":"NC_000007.14","pos":150999023,"ref":"T","alt":"G"}],
 "ambiguous":false,"warnings":[]}
```

**Batch decode** — set `Content-Type: application/json`.

```bash
curl -X POST localhost:4567/decode \
     -H 'Content-Type: application/json' \
     -d '{"hgvs":["ALDH2:p.Glu504Lys","ALDH2:c.1510G>A","NM_000690:c.1510G>A"]}'
```

Each item is decoded independently; per-item errors are returned as `{"input":...,"error":...}`.

## Background

### Architecture

```
                prepare data/  (manual download + symlink, see Quick start)
                        │
     ┌──────────────────┴───────────────────┐
     ▼                                        ▼
 cdot JSON (RefSeq+Ensembl, GRCh38)     GRCh38 FASTA + .fai
 transcript models: exon blocks,        (RefSeq accessions as
 CDS, strand, MANE/canonical tags,      contig ids: NC_0000..)
 transcript↔genome gap CIGAR
     │                                        │
     └──────────────┬─────────────────────────┘
                    ▼
              Decoder (lib/hgvs2vcf)
   Parser → Transcripts(resolve) → Mapper(c.→g.) → Sequence(ref/alt + normalize)
                    ▼
             Sinatra API  (GET/POST /decode)
```

### Why cdot as the transcript backbone

`biocommons/hgvs` + UTA is the reference stack but needs Postgres and aligns sequences itself. **[cdot](https://github.com/SACGF/cdot)** (`SACGF/cdot`, MIT) pre-converts the RefSeq GFF and Ensembl GTF into flat JSON with exactly what a projector needs, and crucially it carries the **transcript↔genome gap** (`"M196 I1 M61 …"`) so RefSeq transcripts that don't align perfectly to the reference map correctly — the single most common source of off-by-N errors in naive GFF-exon projectors. Latest data release at time of writing: **`data_v0.2.33`** (2026-06-26; Ensembl 116, RefSeq RS_2025_08 / GCF_000001405.40 GRCh38.p14), with MANE tags baked in.

cdot exon tuple: `[g_start(0), g_end(excl), exon_number, cdna_start(1), cdna_end(1), gap|null]`. Exons are listed in **genomic** order, so the mapper indexes them by cDNA interval, not array position (matters for the minus strand).

### Data sources

All sources are refreshed periodically and **never queried live** at decode time:

| Role | Source | Used at |
|------|--------|---------|
| Transcript models (RefSeq + Ensembl, GRCh38) | [cdot releases](https://github.com/SACGF/cdot/releases) | index build |
| Reference genome | [NCBI GCF_000001405.40_GRCh38.p14](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/) — `GCF_000001405.40_GRCh38.p14_genomic.fna.gz` (RefSeq contig accessions match cdot `contig`) | serve |
| MANE summary (audit / independent symbol map) | [NCBI MANE_human/current](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/) | auditing only |
| HGNC complete set (optional) | symbol aliases / withdrawn symbols | not used by default |

The running app reads only the LMDB index and the genome FASTA. cdot JSON is needed at index-build time; the MANE file is kept for provenance.

### Three decisions to be aware of

1. **HGVSp is intrinsically ambiguous.** An amino-acid substitution can map to several nucleotide variants. The reference codon (read from the transcript) disambiguates *which* reference codon, but not always *which* target codon:
   - `ALDH2:p.Glu504Lys`, ref codon `GAG` → unique `c.1510G>A` (GAG→AAG).
   - `NOS3:p.Asp298Glu`, ref codon `GAT` → **two** minimal SNVs: `c.894T>A` and `c.894T>G`. (The known rs1799983 is T>G, but HGVSp alone can't say.)
   The decoder returns **all** minimal candidates and sets `ambiguous: true`; candidates needing ≥2 nt changes are flagged `requires_mnv` and emitted as an MNV/delins. Decide your policy: return all, prefer single-nt, or refuse.

2. **MANE Select as the default for a bare symbol / bare accession.** Resolved from cdot tags; falls back to Ensembl-canonical when a gene has no MANE (non-coding, some patch genes). `ALDH2` → `NM_000690.4`, `NOS3` → `NM_000603.x`. Bare `NM_000690` → highest loaded version.

3. **ENST works the same way, and is actually easier.** Ensembl/GENCODE transcripts are derived from the genome, so `gap` is always null — the minus-strand/splice arithmetic is the whole story. Load the cdot *ensembl* file alongside refseq and ENST resolves identically.

### Gotchas (verify against your truth set)

- **CDS datum → cDNA offset.** cdot `start_codon`/`stop_codon` are transcript offsets; the code takes `c.1 = cdna(start_codon + 1)`. The exact ±1 and the UTR (`c.-N`, `c.*N`) and intron (`c.N+M`) conventions are the classic place to be off by one — validate before trusting.
- **Strand.** For the minus strand, `ref`/`alt` from the HGVS are transcript-oriented and get complemented into genome orientation before the ref-allele check and normalization.
- **Reference-allele check.** Every SNV verifies the implied ref base against the genome FASTA and emits a `warnings` entry on mismatch (wrong transcript version, liftover build mix-up, or a real transcript/genome discrepancy).
- **Normalization.** `Normalize.left_align` reproduces `bcftools norm` (trim + left-shift). HGVS 3′-shifts, VCF left-aligns, so indels *will* move.

### On-disk index (LMDB) — no JSON in RAM at serve time

The cdot JSON is parsed **once, offline**, into an LMDB environment (`bin/build_index.rb`). The serving app reads one small packed record per query straight from the mmap — the giant JSON never touches the app's heap.

Two named DBIs in one environment:

- `kv` — resolution + metadata: `mane:<SYMBOL>`, `canon:<SYMBOL>`, `latest:<BASE>`, `meta:contigs` (interned id→name), `meta:build`.
- `tx` — `tx_id` → packed coordinate record.

**Packed record layout** (little-endian; see `lib/hgvs2vcf/packing.rb`):

```
magic u8 | strand u8 | coding u8 | contig_id u16 |
cds_start_i u32 | cds_end_i u32 | gene_len u16 | gene[gene_len] |
n_exons u16 |
  n_exons × { g_start u32, g_end u32, cdna_start u32, cdna_end u32, gap_len u16, gap[gap_len] }
```

Query = one mmap'd record read + binary search over the exon cDNA intervals + arithmetic (gap CIGAR walked only for the rare transcripts that carry one). Build is a single LMDB write transaction (one atomic commit). At serve time the app opens the index with `nolock: true` (no writes, including to `lock.mdb`), so many threads/processes can share the same static index.

## Tests

```bash
ruby spec/engine_test.rb   # parser, HGVSp back-translation, mapper, normalization
ruby spec/lmdb_test.rb     # index build → LMDB → decode
```

## Validating data refreshes

Correctness depends on the data release. After updating cdot or the reference genome, spot-check against known truth:

- ClinVar `variant_summary.txt` or cdot's `clinvar_hgvs_*.tsv` samples
- VariantValidator / Mutalyzer for hard cases (minus-strand genes, gap transcripts, exon boundaries, `c.-N` / `c.*N`, etc.)

A sudden jump in mismatch rate usually means a transcript version or MANE reassignment, not a code bug.

## Environment variables

| Variable | Default | Description |
|----------|---------|-------------|
| `INDEX_DIR` | `data/index` (local) / `/data/index` (container) | LMDB index path |
| `GENOME_FASTA` | `data/GRCh38.fna` / `/data/GRCh38.fna` | Reference FASTA (requires `.fai`) |
| `BUILD` | `GRCh38` | Genome build name |
| `BIND` | `0.0.0.0` | Bind address |
| `PORT` | `4567` | Listen port |
