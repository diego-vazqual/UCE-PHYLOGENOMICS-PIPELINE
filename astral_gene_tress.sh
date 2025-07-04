#!/bin/bash

# Create output folders
mkdir -p raxml_out/{bestTree,info,log,result,parsimonyTree}

input_dir="mafft-phylip-internal-trimmed-gblocks-clean-50p"
output_dir="$(pwd)/astral_raxml_output"

TOTAL_THREADS=24
THREADS_PER_JOB=6
MAX_JOBS=$((TOTAL_THREADS / THREADS_PER_JOB))

RAxML_BIN="/home/intern/miniconda3/bin/raxmlHPC-PTHREADS"

export THREADS_PER_JOB RAxML_BIN output_dir

# Run RAxML in parallel
ls "$input_dir"/*.phylip | parallel -j "$MAX_JOBS" --eta --env THREADS_PER_JOB --env RAxML_BIN --env output_dir 'bash -c "
    aln=\"{}\"  
    base=\$(basename \"\$aln\" .phylip)
    seed=\$((1313 + RANDOM % 10000))

    \"\$RAxML_BIN\" -T \"\$THREADS_PER_JOB\" -N 10 -p \"\$seed\" \
        -s \"\$aln\" -m GTRGAMMA \
        -n \"\$base.genetree\" \
        -w \"\$output_dir/result\"

    mv \"\$output_dir\"/result/RAxML_bestTree.\$base.genetree \"\$output_dir\"/bestTree/
    mv \"\$output_dir\"/result/RAxML_info.\$base.genetree.RUN.* \"\$output_dir\"/info/ 2>/dev/null
    mv \"\$output_dir\"/result/RAxML_log.\$base.genetree.RUN.* \"\$output_dir\"/log/ 2>/dev/null
    mv \"\$output_dir\"/result/RAxML_parsimonyTree.\$base.genetree \"\$output_dir\"/parsimonyTree/ 2>/dev/null
    mv \"\$output_dir\"/result/RAxML_result.\$base.genetree.RUN.* \"\$output_dir\"/result/ 2>/dev/null

    echo \"Tree generated for \$base\"
"'
