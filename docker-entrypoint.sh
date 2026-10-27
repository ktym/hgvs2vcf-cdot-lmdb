#!/usr/bin/env bash
set -euo pipefail

: "${DATA:=/data}"
: "${BUILD:=GRCh38}"
# GRCh38 keeps its original paths; other builds get their own, so both can share /data.
if [[ "$BUILD" == "GRCh38" ]]; then
  : "${INDEX_DIR:=$DATA/index}"
else
  : "${INDEX_DIR:=$DATA/index-$BUILD}"
fi
: "${GENOME_FASTA:=$DATA/$BUILD.fna}"
: "${MANE_SUMMARY:=$DATA/MANE.GRCh38.summary.txt.gz}"
export INDEX_DIR GENOME_FASTA BUILD

# Ruby's plain-FASTA faidx reader needs an uncompressed .fna + .fai. This is a
# one-time bridge; the Rust rewrite will read the .fna.gz directly via htslib
# (bgzf) and drop this decompression step.
prep_fasta() {
  if [[ ! -f "$GENOME_FASTA" ]]; then
    echo "[entrypoint] decompressing $DATA/$BUILD.fna.gz -> $GENOME_FASTA (one-time, ~3GB)"
    gunzip -kc "$DATA/$BUILD.fna.gz" > "$GENOME_FASTA"
  fi
  if [[ ! -f "$GENOME_FASTA.fai" ]]; then
    echo "[entrypoint] samtools faidx $GENOME_FASTA (one-time)"
    samtools faidx "$GENOME_FASTA"
  fi
}

cmd="${1:-serve}"; shift || true
case "$cmd" in
  serve)
    [[ -d "$INDEX_DIR" ]] || { echo "no LMDB index at $INDEX_DIR — run 'index' first" >&2; exit 1; }
    prep_fasta
    exec bundle exec ruby app.rb
    ;;
  index)
    # Both files: only the ensembl one carries the canonical tags that back the
    # fallback for genes without a MANE Select.
    mane=()
    if [[ "$BUILD" != "GRCh38" ]]; then
      # cdot tags MANE only on GRCh38; elsewhere it must come from the summary.
      [[ -f "$MANE_SUMMARY" ]] || { echo "$BUILD needs the MANE summary at $MANE_SUMMARY" >&2; exit 1; }
      mane=(--mane "$MANE_SUMMARY")
    fi
    exec ruby bin/build_index.rb --out "$INDEX_DIR" --build "$BUILD" ${mane[@]+"${mane[@]}"} \
      "$DATA/cdot.refseq.$BUILD.json.gz" "$DATA/cdot.ensembl.$BUILD.json.gz"
    ;;
  prep-fasta)
    prep_fasta
    ;;
  *)
    exec "$cmd" "$@"
    ;;
esac
