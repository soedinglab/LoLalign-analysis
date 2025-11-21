#!/bin/bash

foldseek easy-search \
    /queryDB \
    /targetDB \
    fs_results/fs_aln \
    tmp_fs \
    --alignment-type 0 \
    --format-mode 4 \
    --max-seqs 2000 \
    --format-output "query,target,lddt,qtmscore,alntmscore,evalue,bits,fident,qstart,qend,tstart,tend,cigar" \
    --threads 128


