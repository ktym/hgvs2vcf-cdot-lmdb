#!/usr/bin/env bash
set -euo pipefail

: "${DATA:=/data}"
: "${INDEX_DIR:=$DATA/index}"
: "${GENOME_FASTA:=$DATA/GRCh38.fna}"
: "${BUILD:=GRCh38}"

# Ruby's plain-FASTA faidx reader needs an uncompressed .fna + .fai. This is a
# one-time bridge; the Rust rewrite will read GRCh38.fna.gz directly via htslib
# (bgzf) and drop this decompression step.
prep_fasta() {
  if [[ ! -f "$GENOME_FASTA" ]]; then
    echo "[entrypoint] decompressing $DATA/GRCh38.fna.gz -> $GENOME_FASTA (one-time, ~3GB)"
    gunzip -kc "$DATA/GRCh38.fna.gz" > "$GENOME_FASTA"
  fi
  if [[ ! -f "$GENOME_FASTA.fai" ]]; then
    echo "[entrypoint] samtools faidx $GENOME_FASTA (one-time)"
    samtools faidx "$GENOME_FASTA"
  fi
}

cmd="${1:-serve}"; shift || true
case "$cmd" in
  serve)
    [[ -d "$INDEX_DIR" ]] || { echo "no LMDB index at $INDEX_DIR — run 'build' first" >&2; exit 1; }
    prep_fasta
    exec bundle exec ruby app.rb
    ;;
  index)
    exec ruby bin/build_index.rb --out "$INDEX_DIR" --build "$BUILD" \
      "$DATA/cdot.refseq.GRCh38.json.gz" "$DATA/cdot.ensembl.GRCh38.json.gz"
    ;;
  prep-fasta)
    prep_fasta
    ;;
  *)
    exec "$cmd" "$@"
    ;;
esac
