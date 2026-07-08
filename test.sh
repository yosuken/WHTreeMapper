#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

TOOL="${WHTREE_MAPPER:-./wh_tree_mapper}"
NCPUS="${NCPUS:-1}"

run_case() {
  local mode="$1"
  local input="$2"
  local outdir="$3"

  echo "==> ${mode}: ${input} -> ${outdir}"
  rm -rf "$outdir"
  "$TOOL" "$mode" -i "$input" -o "$outdir" -n "$NCPUS"

  if [[ ! -s "${outdir}/detected.tsv" ]]; then
    echo "ERROR: ${outdir}/detected.tsv was not created or is empty" >&2
    exit 1
  fi

  local entries
  entries="$(tail -n +2 "${outdir}/detected.tsv" | wc -l)"
  if [[ "$entries" -lt 1 ]]; then
    echo "ERROR: no entries found in ${outdir}/detected.tsv" >&2
    exit 1
  fi

  echo "==> ${outdir}/detected.tsv: ${entries} entries"
}

run_case --prot test_seq/BFU66023.faa test/BFU66023
run_case --nucl test_seq/LC852818.fa test/LC852818

echo "==> test passed"
