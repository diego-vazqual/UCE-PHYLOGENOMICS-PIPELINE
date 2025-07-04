#!/bin/bash

for i in mafft--fastas-internal-no-trimmed-gblocks-clean/*.fasta; do
        file=$(basename $i)
        uce=${file%.fasta}
        sed -i "s.>.>${uce}_.g" $i      
done
