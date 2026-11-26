#!/bin/bash
# Quality filtering and adapter trimming of paired-end reads with fastp.
# Usage: bash fastp.sh <input_dir> <output_dir>
#   <input_dir>  one subfolder per sample, each with its R1 and R2 files (*_R1*.fastq.gz, *_R2*.fastq.gz)
#   <output_dir> one subfolder per sample with <sample>-READ1.fastq, <sample>-READ2.fastq and the fastp reports
# Requires fastp in the PATH (e.g. activate the phyluce-1.7.2 environment before running the script).

input="$1"
output="$2"
threads=10

if [ -z "$input" ] || [ -z "$output" ]; then
    echo "Usage: bash fastp.sh <input_dir> <output_dir>"; exit 1
fi
mkdir -p "$output"

for sample_dir in "$input"/*/; do
    sample=$(basename "$sample_dir")
    echo "PROCESSING: $sample"

    r1=$(find "$sample_dir" -maxdepth 1 -type f -name "*_R1*.fastq.gz" | head -1)
    r2=$(find "$sample_dir" -maxdepth 1 -type f -name "*_R2*.fastq.gz" | head -1)

    if [ -n "$r1" ] && [ -n "$r2" ]; then
        out="$output/$sample"
        mkdir -p "$out"
        fastp -i "$r1" -I "$r2" \
              -o "$out/${sample}-READ1.fastq" -O "$out/${sample}-READ2.fastq" \
              --trim_poly_g --correction --overrepresentation_analysis \
              --html "$out/${sample}_fastp.html" --json "$out/${sample}_fastp.json" \
              --thread "$threads"
    else
        echo "[ERROR] R1 and/or R2 not found for $sample"
    fi
done
