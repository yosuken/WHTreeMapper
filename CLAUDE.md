# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

WHTreeMapper is a Ruby-based bioinformatics tool for detecting Walker Homology (WH) regions in protein sequences and performing diamond blastp analysis against a WH reference database. It uses HMM-based prefiltering against 50 WH clade profiles to extract homologous regions, then identifies best hits via diamond blastp. Nucleotide FASTA input can be translated to proteins with Prodigal metagenome mode.

This tool was derived from PiPP (Pipeline for Phylogenetic Placement) and specialized for WH region analysis.

## Running the Tool

```bash
# Basic usage
./wh_tree_mapper --prot -i <query_protein_fasta> -o <output_dir>

# Nucleotide FASTA input
./wh_tree_mapper --nucl -i <query_nucleotide_fasta> -o <output_dir> --codon-table 11

# With multiple CPUs
./wh_tree_mapper --prot -i "queries/*.faa" -o results -n 4

# Overwrite existing output
./wh_tree_mapper --prot -i query.faa -o results --overwrite -n 2
```

## Architecture Overview

### Core Components

1. **Main Script (`wh_tree_mapper`)**: Ruby CLI that validates dependencies, parses arguments, and invokes the Rake workflow

2. **Rake Workflow (`WHTreeMapper.rake`)**: Orchestrates the pipeline through 10 sequential tasks:
   - Query validation and preprocessing (01-1a-A, 01-1a-B)
   - Reference package validation (01-1b-A, 01-1b-B)
   - HMM-based prefiltering with hmmsearch (01-2a, 01-2b)
   - Detected sequence extraction (01-2c, 01-2d)
   - Region merging for diamond query (02-1)
   - Diamond BLASTP analysis (02-2)

3. **Ruby Processing Scripts (`script/`)**: Utilities for sequence validation, HMM result parsing, and refhmm validation

### Bundled Data

- `refhmm/`: 50 WH clade HMM profiles (A, A1, A2, B, B1, ..., H22), each containing `.fa`, `.hmm`, `.fasttree.newick`
- `db/wh.dmnd`: Diamond database of WH representative sequences

### Key Dependencies

- ruby (>= 2.0)
- hmmer (>= 3.0) for hmmsearch
- diamond (>= 2.1.9) for blastp
- prodigal for nucleotide FASTA input
- GNU parallel (when --ncpus > 1)

### Pipeline Workflow

```
Input: protein FASTA or nucleotide FASTA → Prodigal (--nucl only) → hmmsearch (50 HMM profiles) → parse hits → extract WH regions
→ merge regions → diamond blastp → output
```

### Output Structure

```
<outdir>/
  detected.tsv                             -- main result (clade annotation + diamond hits)
  diamond/dmnd.blastp.out                  -- raw diamond blastp output
  diamond/query/region.faa                 -- merged query for diamond
  hmm_hits/<clade>/seq/{whole,region}.fa   -- detected sequences per clade
  prefilter/                               -- intermediate hmmsearch results
  log/                                     -- execution logs
```

### Diamond BLASTP Parameters (default)

- `--ultra-sensitive --id 40 --subject-cover 80 --max-target-seqs 1 --dbsize 1e9 --evalue 1e-5`
- Output format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp scovhsp

### HMM Detection Parameters (default)

- `--evalue 1e-5 --evaluedom 1e-2 --minhmmcov 0.8 --minhmmcovdom 0.2`

## Development Notes

- Uses Ruby 2.0+ with Rake for task orchestration
- Batch job generation via WriteBatch/RunBatch lambdas for parallel execution
