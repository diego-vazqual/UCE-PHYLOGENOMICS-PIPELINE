#!/usr/bin/env python3
 
# ============================================================
# remove_empty_sequences.py
#
# Remove from every NEXUS alignment the sequences made only of
# gaps or missing data (-, ?, N). These appear after Gblocks
# when a taxon only had sequence in the columns that were
# removed. RAxML refuses to build a tree from an alignment
# with such sequences, and phyluce counts them as present when
# building the occupancy matrices.
#
# Usage:
#   python3 remove_empty_sequences.py <input_dir> <output_dir>
#
# Arguments:
#   input_dir    Directory with the alignments in NEXUS format
#                (.nexus), e.g. the output of
#                phyluce_align_remove_locus_name_from_files
#   output_dir   Directory where the cleaned alignments will be
#                saved (same file names). Created if needed.
#
# Example:
#   python3 remove_empty_sequences.py \
#       mafft-nexus-internal-no-trimmed-gblocks-clean \
#       mafft-nexus-internal-no-trimmed-gblocks-clean-noempty
#
# Output (in the working directory):
#   empty_sequences_removed.tsv   loci and taxa removed
#
# Requirements: Biopython
# ============================================================
 
import os
import sys
from Bio.Nexus import Nexus
 
if len(sys.argv) != 3:
    raise SystemExit(
        "Usage: python3 remove_empty_sequences.py <input_dir> <output_dir>"
    )
 
INPUT_DIR = sys.argv[1]
OUTPUT_DIR = sys.argv[2]
LOG = "empty_sequences_removed.tsv"
 
EMPTY_CHARS = set("-?N")
 
os.makedirs(OUTPUT_DIR, exist_ok=True)
 
nexus_files = sorted(f for f in os.listdir(INPUT_DIR) if f.endswith(".nexus"))
 
if not nexus_files:
    raise SystemExit("ERROR: No .nexus files found in " + INPUT_DIR)
 
n_loci_changed = 0
n_removed = 0
 
with open(LOG, "w") as log:
    log.write("Locus\tTaxa_before\tTaxa_after\tRemoved\n")
 
    for nx_file in nexus_files:
        nx = Nexus.Nexus()
        nx.read(os.path.join(INPUT_DIR, nx_file))
 
        empty = [t for t, seq in nx.matrix.items()
                 if not set(str(seq).upper()) - EMPTY_CHARS]
 
        # Same file name as the input: phyluce uses it as the locus name
        nx.write_nexus_data(filename=os.path.join(OUTPUT_DIR, nx_file),
                            delete=empty)
 
        if empty:
            n_loci_changed += 1
            n_removed += len(empty)
            locus = nx_file.replace(".nexus", "")
            n = len(nx.matrix)
            log.write(f"{locus}\t{n}\t{n - len(empty)}\t{','.join(empty)}\n")
 
print(f"Alignments processed: {len(nexus_files)}")
print(f"Alignments with empty sequences: {n_loci_changed}")
print(f"Empty sequences removed: {n_removed}")
print(f"Details: {LOG}")
print(f"Output: {OUTPUT_DIR}")
