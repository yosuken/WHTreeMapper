#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

TOOL="${WHTREE_MAPPER:-./wh_tree_mapper}"
NCPUS="${NCPUS:-1}"

run_case() {
  local mode="$1"
  local input="$2"
  local outdir="$3"
  local expected_query="$4"

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

  local query
  query="$(awk 'NR == 2 { print $1 }' "${outdir}/detected.tsv")"
  if [[ "$query" != "$expected_query" ]]; then
    echo "ERROR: expected query '${expected_query}', got '${query}' in ${outdir}/detected.tsv" >&2
    exit 1
  fi
}

run_case --prot test_seq/BFU66023.faa test/BFU66023 BFU66023.1
run_case --nucl test_seq/LC852818.fa test/LC852818 LC852818.1_1

echo "==> test passed"
