#!/bin/bash

mkdir -p genes_tree

input_dir="mafft-phylip-nexus-internal-no-trimmed-gblocks-clean"
output_dir="$(pwd)/genes_tree"

# Total number of available threads
TOTAL_THREADS=24
# Threads per RAxML instance
THREADS_PER_JOB=6
# Maximum number of simultaneous processes
MAX_JOBS=$((TOTAL_THREADS / THREADS_PER_JOB))

# Path to the multithreaded RAxML executable
RAxML_BIN="/home/intern/miniconda3/bin/raxmlHPC-PTHREADS"

# Export variables for GNU parallel environment
export THREADS_PER_JOB output_dir RAxML_BIN

# Run RAxML in parallel
ls "$input_dir"/*.phylip | parallel -j "$MAX_JOBS" --eta --env THREADS_PER_JOB --env output_dir --env RAxML_BIN '
    base=$(basename {} .phylip)
    $RAxML_BIN -T $THREADS_PER_JOB -s {} -n $base -m GTRGAMMA -p 12345 -x 12345 -# 100 -f a -w $output_dir
    echo "Tree generated for $base using $THREADS_PER_JOB threads"
'
