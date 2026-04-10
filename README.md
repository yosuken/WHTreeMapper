
# WHTreeMapper - WH Region Detection and Diamond BLASTP Analysis Tool

[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](/LICENSE)

## Overview

WHTreeMapper detects Walker Homology (WH) regions in query protein sequences using HMM-based prefiltering against 50 WH clade profiles, then performs diamond blastp against a WH reference database.

### Workflow

1. Prefilter query protein sequences by HMM similarity detection (hmmsearch) against 50 WH clade profiles
2. Extract detected WH regions from query sequences
3. Merge detected regions and run diamond blastp against WH reference database
4. Generate `detected.tsv` with clade annotation and diamond blastp results

## Install

### 1. Create conda environment

Using `environment.yaml`:

```bash
# Using micromamba (recommended)
micromamba create -n WHTreeMapper -f environment.yaml -y
micromamba activate WHTreeMapper

# Or using mamba
mamba env create -f environment.yaml
mamba activate WHTreeMapper

# Or using conda
conda env create -f environment.yaml
conda activate WHTreeMapper
```

Dependencies (installed via environment.yaml):
- ruby (>= 3.2)
- hmmer (3.3.2)
- diamond (2.1.9)
- GNU parallel (>= 20230822)

### 2. Set up bundled data

The following data must be placed in the WHTreeMapper directory:

- `refpkg/` - HMM profiles for 8 WH clades (A-H) and 42 subclades
- `db/wh.dmnd` - Diamond database of WH representative sequences
- `db/seq_info.tsv` - Sequence-to-clade mapping table

## Usage

```bash
micromamba activate WHTreeMapper
./WHTreeMapper [options] -q <query protein fasta(s)> -o <output dir>
```

### Options

```
[File/directory]
    -q, --query FILE(S)              Query sequence file(s) (protein fasta, can be gzipped) [required]
    -o, --outdir PATH                Output directory [required]
        --[no-]overwrite             Overwrite output directory (default: no overwrite)

[HMM prefilter]
    -e, --evalue NUM                 E-value threshold of hmmsearch (default: 1e-5)
        --evaluedom NUM              Domain E-value threshold of hmmsearch (default: 1e-2)
        --minhmmcov FLOAT            Minimum HMM coverage fraction in linked result (default: 0.8)
        --minhmmcovdom FLOAT         Minimum HMM coverage fraction in domain-level result (default: 0.2)

[Diamond BLASTP]
        --dmnd-evalue NUM            E-value threshold for diamond blastp (default: 1e-5)
        --dmnd-id NUM                Minimum identity % for diamond blastp (default: 40)
        --dmnd-subject-cover NUM     Minimum subject coverage % for diamond blastp (default: 80)
        --dmnd-max-target-seqs INT   Max target sequences per query for diamond blastp (default: 1)

[Computation]
    -n, --ncpus INT                  Number of CPUs to use (default: 1)

[General]
    -h, --help                       Show this help message
    -v, --version                    Show version
```

### Examples

```bash
# Basic usage
./WHTreeMapper -q query.faa -o results

# Multiple query files with 4 CPUs
./WHTreeMapper -q "queries/*.faa" -o results -n 4

# Overwrite existing output
./WHTreeMapper -q query.faa -o results --overwrite

# Custom HMM and diamond thresholds
./WHTreeMapper -q query.faa -o results -e 1e-3 --minhmmcov 0.5 --dmnd-id 30 --dmnd-subject-cover 60
```

## Output

```
<outdir>/
  detected.tsv                             -- main result
  diamond/dmnd.blastp.out                  -- raw diamond blastp output
  diamond/query/region.faa                 -- merged WH regions used as diamond query
  hmm_hits/<clade>/seq/region.fa           -- detected WH region sequences (per clade)
  hmm_hits/<clade>/seq/whole.fa            -- full-length sequences with WH hits (per clade)
```

### detected.tsv

Main output file. Each row corresponds to a query sequence with a WH region detected by hmmsearch and matched via diamond blastp to the WH reference database.

| # | Column | Description |
|---|--------|-------------|
| 0 | query | Query sequence ID (with region suffix) |
| 1 | detected_wh_region | Extracted WH region coordinates in the query (e.g., `10-109`) |
| 2 | clade | WH clade of the best diamond hit (A, B, C, ..., H) |
| 3 | subclade | WH subclade (A1, A2, B1, ..., H22). Empty if hit is at the clade level |
| 4 | index_in_tree | Index of the matched leaf in the WH reference tree |
| 5 | leaf_name | Leaf name in the WH reference tree |
| 6 | qseqid | Diamond: query sequence ID |
| 7 | sseqid | Diamond: subject sequence ID |
| 8 | pident | Diamond: percentage of identical matches |
| 9 | length | Diamond: alignment length |
| 10 | mismatch | Diamond: number of mismatches |
| 11 | gapopen | Diamond: number of gap openings |
| 12 | qstart | Diamond: start of alignment in query |
| 13 | qend | Diamond: end of alignment in query |
| 14 | sstart | Diamond: start of alignment in subject |
| 15 | send | Diamond: end of alignment in subject |
| 16 | evalue | Diamond: expect value |
| 17 | bitscore | Diamond: bit score |
| 18 | qlen | Diamond: query sequence length |
| 19 | slen | Diamond: subject sequence length |
| 20 | qcovhsp | Diamond: query coverage per HSP (%) |
| 21 | scovhsp | Diamond: subject coverage per HSP (%) |

### Default parameters

- HMM detection: `-e 1e-5 --evaluedom 1e-2 --minhmmcov 0.8 --minhmmcovdom 0.2`
- Diamond BLASTP: `--dmnd-evalue 1e-5 --dmnd-id 40 --dmnd-subject-cover 80 --dmnd-max-target-seqs 1`
- Diamond BLASTP (fixed): `--ultra-sensitive --dbsize 1e9`

## Test

```bash
micromamba activate WHTreeMapper
cd test/wh
./run_test.sh
```

Or using make directly:

```bash
cd test/wh
make test    # run test
make clean   # remove test output
```

## Citation
```
```
