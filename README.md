
# WHTreeMapper - Detection of WH RIPs and their classification based on the WH reference tree

[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](/LICENSE)

## Overview

WHTreeMapper detects winged-helix (WH) regions of WH RIPs (replication initiation proteins) in query protein sequences using HMM-based prefiltering against 50 WH clade profiles (8 clades and 42 subclades), then performs diamond blastp against a WH reference database. Nucleotide FASTA input can also be translated to proteins with Prodigal metagenome mode.

### Workflow

1. If nucleotide FASTA is given, predict protein sequences with Prodigal `-p meta`
2. Prefilter query protein sequences by HMM similarity detection (hmmsearch) against 50 WH clade profiles (8 clades and 42 subclades)
3. Extract detected WH regions from query sequences
4. Merge detected regions and run diamond blastp against WH reference database
5. Generate `detected.tsv` with clade annotation and diamond blastp results

## Install

### 1. Create conda environment and install command

Using `install.sh`:

```bash
bash install.sh
micromamba activate WHTreeMapper
wh_tree_mapper -h
```

`install.sh` uses `environment.yaml` with `micromamba`, `mamba`, or `conda`, then installs a `wh_tree_mapper` launcher into the environment's `bin/` directory.

Manual environment creation is also possible, but `wh_tree_mapper` will not be added to `PATH` automatically:

```bash
micromamba create -n WHTreeMapper -f environment.yaml -y
micromamba activate WHTreeMapper
./wh_tree_mapper -h
```

Dependencies (installed via environment.yaml):
- ruby (>= 3.2)
- hmmer (3.3.2)
- diamond (2.1.9)
- prodigal (>= 2.6.3)
- GNU parallel (>= 20230822)

### 2. Set up bundled data

The following data must be placed in the WHTreeMapper directory:

- `refpkg/` - HMM profiles for 8 WH clades (A-H) and 42 subclades
- `db/wh.dmnd` - Diamond database of WH representative sequences
- `db/seq_info.tsv` - Sequence-to-clade mapping table

## Usage

```bash
micromamba activate WHTreeMapper
wh_tree_mapper (--prot | --nucl) -i <input fasta(s)> -o <output dir>
```

Exactly one of `--prot` or `--nucl` is required.

### Options

```
[File/directory]
    -i, --input FILE(S)              Input sequence file(s) (FASTA, can be gzipped) [required]
    -o, --outdir PATH                Output directory [required]
        --[no-]overwrite             Overwrite output directory (default: no overwrite)

[Input]
        Exactly one of --prot or --nucl is required.
        --prot                       Treat input as protein FASTA
        --nucl                       Treat input as nucleotide FASTA and predict proteins with prodigal -p meta
        --codon-table INT            Translation table for prodigal -g when --nucl is used (default: 11)

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
wh_tree_mapper --prot -i query.faa -o results

# Multiple query files with 4 CPUs
wh_tree_mapper --prot -i "queries/*.faa" -o results -n 4

# Nucleotide FASTA input translated by Prodigal metagenome mode
wh_tree_mapper --nucl -i contigs.fna -o results --codon-table 11

# Overwrite existing output
wh_tree_mapper --prot -i query.faa -o results --overwrite

# Custom HMM and diamond thresholds
wh_tree_mapper --prot -i query.faa -o results -e 1e-3 --minhmmcov 0.5 --dmnd-id 30 --dmnd-subject-cover 60
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
- Prodigal for nucleotide input: `-p meta -g 11`
- Diamond BLASTP: `--dmnd-evalue 1e-5 --dmnd-id 40 --dmnd-subject-cover 80 --dmnd-max-target-seqs 1`
- Diamond BLASTP (fixed): `--ultra-sensitive --dbsize 1e9`

## Test

```bash
micromamba activate WHTreeMapper
cd test/wh
./run_test.sh
```

## Citation

Nishimura Y, Kaneko K, Kamijo T, Isogai N, Tokuda M, Xie H, Tsuda Y, Hirabayashi A, Moriuchi R, Dohra H, Kimbara K, Suzuki-Minakuchi C, Nojiri H, Suzuki H, Suzuki M, Shintani M. A replication-centered phylogeny illuminates the evolutionary landscape of bacterial plasmids. *bioRxiv* (2024). doi: [10.1101/2024.09.03.610885](https://doi.org/10.1101/2024.09.03.610885)
