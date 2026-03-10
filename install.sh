v=PiPP_v0.3.0 && micromamba create -n $v -c conda-forge -c bioconda ruby=3.4.5 hmmer=3.4 parallel=20250822 gappa=0.9.0 pplacer=1.1.alpha19 mafft=7.520 fasttree=2.2.0 epa-ng=0.3.8 python=3.12 -y && micromamba activate $v && pip install apples taxtastic
# pip install weblogo
