
# PiPP - a Pipeline for Phylogenetic Placement

[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](/LICENSE)
[![size](https://img.shields.io/github/size/webcaetano/craft/build/phaser-craft.min.js.svg)]()

## currently PiPP is beta version. Any specification might be changed in a future version.
PiPP is developed as a tool for phylogenetic placement onto a clade or taxonomy defined phylogenetic tree through procedures below. Very large queries are acceptable.

## install (use conda environment)

### [1] make conda environment and install packages
```
$ conda create -n PiPP -y && conda activate PiPP
```

### [2] install packages
```
$ conda install -y -c conda-forge ruby=2.7.2 parallel=20210622
$ conda install -y -c bioconda hmmer=3.3.2 gappa=0.7.1 mafft=7.490
$ conda install -y -c bioconda pplacer=1.1.alpha19   ### <-- does not work for OSX environment
```

## usage 
```
### PiPP ver 0.1.0 (2020-05-11) ###

[description]
PiPP - Pipeline for phylogenetic placement.
PiPP is developed as a tool for phylogenetic placement onto a clade or taxonomy defined phylogenetic tree through procedures below. Very large queries are acceptable.

1. prefilter query sequences by similarity detection (hmmsearch). Queries and references should be protein sequences at this moment.
2. align query sequences to a given reference alignment ('mafft --add' [or 'hmmalign': not implemented yet])
3. perform phylogenetic placement ('pplacer' [or 'epa-ng': not implemented yet]) with efficient parallelization ('gappa prepare chunkify' and 'gappa prepare unchunkify')
4. analysis of placed sequences
  a. assign clade/taxonomy and generate statistics ('gappa examine assign')
  b. extract placed sequences and placement file for each  clade/taxonomy ('gappa prepare extract')
  c. extract placed sequences and placement file for each  clade/taxonomy ('gappa prepare extract')

[usage]
$ PiPP [options] -q <query fasta(s)> -r <refpkg dir(s)> -o <output dir>

[dependencies]
- hmmer (ver >= 3.0)
- mafft (tested by 7.453)
- gappa (tested by 0.6.0)
- pplacer (tested by 1.1.alpha19)
- ruby (ver >= 2.0)

[options]
  (general)
    -h, --help
    -v, --version

  (file/directory)
    -q, --query      [file(s)] (required)  -- query sequence file(s) (currently protein only)
                                              Multiple query file can be specfied with wildcard (e.g., 'dir/*.fa', quote required) or comma separated values (e.g., A.fa,B.fa).
                                              File names should be unique (used as output labels).
    -r, --refpkg     [dir(s)] (required)   -- reference package(s) made by taxtastic as 'taxit create -P <name> --tree-file <newick file> --tree-stats <file>', for example.
                                              The directory must contain additional 2 files:
                                              - aligned sequence (.fa, .mfa, or .fasta),
                                              - its hmm model (.hmm, made by 'hmmbuild')
                                              For full functionality, it is recommended that the directory contains additional 2 files:
                                              - taxon file ('taxon.tsv' -- 1. sequence_id, 2. clade or taxonomy) used by '--taxon-file <tsv>' of 'gappa examine assign'
                                              - alignment check file ('position.tsv' -- 1. label, 2. position(s) in alignment, which can be multiple per a label by using ',' as a separator)
                                                to check important positions in sequence alignment

                                              Multiple Reference packages can be specfied with wildcard (e.g., 'dir/*', quate required) or comma separated values (e.g., A,B).
                                              Directory names should be unique (used as output labels).
                                              Note that each query is placed onto at most one tree that is most significantly related among refpkgs (in terms of hmmsearch evalue).

                                              Further information
                                              - taxtastic: https://github.com/fhcrc/taxtastic
                                              - gappa examine assign: https://github.com/fhcrc/taxtasti://github.com/lczech/gappa/wiki/Subcommand:-assign
    -o, --outdir     [path] (required)     -- output directory (should not exist unless '--overwrite' is specified)
    --overwrite      [bool]                -- overwrite output directory

  (prefilter)
    -e, --evalue     [num] (default: 1e-5) -- evalue threshold of hmmsearch
    --minseqlen      [int] (default: 0)    -- minimum query aa length
    --minhmmlen      [int] (default: 0)    -- minimum hmm length in alignment of hmmsearch
    --trim-option    [merge or largest] (default: merge) -- take merged regions of hmmsearch hits (if multiple) or take largest hit

  (align)
    --mafft-method   [FFT-NS-2 or E-INS-i] (default: FFT-NS-2) -- mafft method used with '--add' and '--keeplength'.
                                                                  It can be either '--retree 2' (FFT-NS-2, faster) or '--genafpair' (E-INS-i, more accurate).
                                                                  For more details, see https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html.

  (extract)
    --extract-levels [int(s)] (default: 0) -- specify level of clade/taxonomy (written in the second column of TSV file in refpkg) to extract each clade/taxonomy placement.
                                              (e.g., 'Eukaryota; Amoebozoa; Myxogastria' has three levels. '--extract-levels 3' means 'Myxogastria' will be extracted.)
                                              Multiple levels can be specified with comma separeted values (e.g., '--extract-levels 2,3').
                                               '0' means all levels in the TSV file. '-1' means extraction is not performed.

  (computation)
    -c, --chunk-size [int] (default: 3000) -- chunk size of 'gappa prepare chunkify'. If the value is lower, smaller batch jobs (that can be parallelized) will be generated.
                                              Note that it is not efective to have very low value (e.g., less than 1,000)
    -n, --ncpus      [int] (default: 1)    -- num CPUs to use (for multicore computation, 'parallel' (GNU parallel) need to be available.)

[output files]
  result/<refpkg name>/{all,each}/{query,alignment,placement,assign,extract,...} -- result of placement and further analysis
```

## citation
```
```
