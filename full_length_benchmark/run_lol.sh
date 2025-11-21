#!/bin/bash

foldseek easy-search \
    /queryDB \
    /targetDB \
    results/fs_lol_aln \
    tmp \
    --alignment-type 3 \
    --format-mode 4 \
    --max-seqs 2000 \
    --format-output "query,target,lddt,qtmscore,alntmscore,evalue,bits,fident,qstart,qend,tstart,tend,cigar" \
    --threads 128


