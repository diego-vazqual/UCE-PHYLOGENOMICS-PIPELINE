#!/bin/bash
# Removal of duplicate read pairs with CD-HIT-DUP.
# Usage: bash cd-hit-dup.sh <input_dir> <output_dir>
#   <input_dir>  output of fastp.sh: one subfolder per sample with <sample>-READ1.fastq and <sample>-READ2.fastq
#   <output_dir> one subfolder per sample with the deduplicated reads, compressed (.fastq.gz)
# Requires cd-hit-dup in the PATH.

input="$1"
output="$2"

if [ -z "$input" ] || [ -z "$output" ]; then
    echo "Usage: bash cd-hit-dup.sh <input_dir> <output_dir>"; exit 1
fi
mkdir -p "$output"

for sample_dir in "$input"/*/; do
    sample=$(basename "$sample_dir")
    echo "PROCESSING: $sample"

    read1="$sample_dir/${sample}-READ1.fastq"
    read2="$sample_dir/${sample}-READ2.fastq"

    if [ -f "$read1" ] && [ -f "$read2" ]; then
        out="$output/$sample"
        mkdir -p "$out"
        # -u 30: compare only the first 30 bases of each read; -m false: do not require identical length
        cd-hit-dup -u 30 -m false \
                   -i "$read1" -i2 "$read2" \
                   -o "$out/${sample}-READ1.fastq" -o2 "$out/${sample}-READ2.fastq"
        gzip "$out/${sample}-READ1.fastq" "$out/${sample}-READ2.fastq"
        rm -f "$out"/*.clstr                     # cluster reports, not needed downstream
        echo "Duplicates removed from $sample"
    else
        echo "[ERROR] READ1 and/or READ2 not found for $sample"
    fi
done
