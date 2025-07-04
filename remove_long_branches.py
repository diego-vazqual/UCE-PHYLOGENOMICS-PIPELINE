import os
import numpy as np
from collections import defaultdict, Counter
from Bio import Phylo
from Bio import SeqIO
from Bio.Nexus import Nexus

# ===============================
# USER CONFIGURATION
# ===============================
TREE_DIR = "genes_tree/RAxML_bipartitions"
NEXUS_DIR = "mafft-nexus-internal-trimmed-gblocks-clean"
OUTPUT_TREES = "clean_genes_trees"
OUTPUT_NEXUS = "clean-mafft-nexus-internal-trimmed-gblocks-clean"
TAXON_FILE = "spms_info.txt"
DISCARDS_LOG = "long_branch_discards.tsv"
SUMMARY_LOG = "summary_report.txt"
NEXUS_LOG = "filtered_nexus_log.tsv"

# ===============================
# FILTERING PARAMETERS
# ===============================
TAXON_SSD_THRESHOLD = 2.7
TAXON_MSD_THRESHOLD = 2.7
LOCUS_SD_THRESHOLD = 2.7

# ===============================
# LOAD TAXA
# ===============================
taxon_to_group = {}
with open(TAXON_FILE) as f:
    next(f)
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) != 2:
            continue
        group, taxon = parts
        taxon_to_group[taxon] = taxon_to_group.get(taxon, group)

all_taxa = list(taxon_to_group.keys())

# ===============================
# SUPPORT FUNCTIONS
# ===============================
def get_branch_lengths(tree, taxa):
    lengths = {}
    for terminal in tree.get_terminals():
        name = terminal.name
        if name in taxa:
            path = tree.get_path(terminal)
            if path:
                lengths[name] = path[-1].branch_length or 0.0
    return lengths

# ===============================
# TREE PROCESSING
# ===============================
for d in [OUTPUT_TREES, OUTPUT_NEXUS]:
    if not os.path.exists(d):
        os.makedirs(d)

print("Reading trees and detecting long branches based on global and group statistics...")
tree_files = [f for f in os.listdir(TREE_DIR) if f.startswith("RAxML_bipartitions.")]

locus_lengths = {}
taxon_branch_lengths_s = defaultdict(list)
taxon_branch_lengths_m = defaultdict(list)
group_counts_global = Counter(taxon_to_group.values())

discarded_by_locus = defaultdict(list)
all_discards = []

for tree_file in tree_files:
    locus = tree_file.split(".")[1]
    path = os.path.join(TREE_DIR, tree_file)
    tree = Phylo.read(path, "newick")
    brlens = get_branch_lengths(tree, all_taxa)

    locus_lengths[locus] = []
    present_taxa = list(brlens.keys())
    present_groups = [taxon_to_group[t] for t in present_taxa]
    group_counts_local = Counter(present_groups)

    for taxon, brlen in brlens.items():
        group = taxon_to_group[taxon]
        is_singleton = group_counts_global[group] == 1
        if is_singleton:
            taxon_branch_lengths_s[taxon].append(brlen)
        else:
            taxon_branch_lengths_m[taxon].append(brlen)
        locus_lengths[locus].append((taxon, brlen))

# Calculate means and SDs per taxon (global)
taxon_stats_s = {t: (np.mean(v), np.std(v)) for t, v in taxon_branch_lengths_s.items() if len(v) > 0}
taxon_stats_m = {t: (np.mean(v), np.std(v)) for t, v in taxon_branch_lengths_m.items() if len(v) > 0}

# Detect long branches per locus
for tree_file in tree_files:
    locus = tree_file.split(".")[1]
    path = os.path.join(TREE_DIR, tree_file)
    tree = Phylo.read(path, "newick")
    brlens = get_branch_lengths(tree, all_taxa)

    values = list(brlens.values())
    locus_mean = np.mean(values)
    locus_sd = np.std(values)

    for taxon, brlen in brlens.items():
        group = taxon_to_group[taxon]
        is_singleton = group_counts_global[group] == 1

        if is_singleton and taxon in taxon_stats_s:
            mean_s, sd_s = taxon_stats_s[taxon]
            if brlen > locus_mean + LOCUS_SD_THRESHOLD * locus_sd and brlen > mean_s + TAXON_SSD_THRESHOLD * sd_s:
                discarded_by_locus[locus].append(taxon)
                all_discards.append(taxon)
        elif taxon in taxon_stats_m:
            mean_m, sd_m = taxon_stats_m[taxon]
            if brlen > locus_mean + LOCUS_SD_THRESHOLD * locus_sd and brlen > mean_m + TAXON_MSD_THRESHOLD * sd_m:
                discarded_by_locus[locus].append(taxon)
                all_discards.append(taxon)

# ===============================
# WRITE DISCARDS OUTPUT
# ===============================
print("Writing discard file...")
with open(DISCARDS_LOG, "w") as out:
    out.write("Locus\tNum_Discarded\tTaxa\n")
    for locus, taxa in discarded_by_locus.items():
        out.write(f"{locus}\t{len(taxa)}\t{','.join(taxa)}\n")

# ===============================
# FILTER NEXUS ALIGNMENTS
# ===============================
print("Filtering NEXUS alignments...")
nexus_files = [f for f in os.listdir(NEXUS_DIR) if f.endswith(".nexus")]

with open(NEXUS_LOG, "w") as log:
    log.write("Locus\tOriginal\tKept\tRemoved\n")

    for nx_file in nexus_files:
        locus = nx_file.replace(".nexus", "")
        to_discard = set(discarded_by_locus.get(locus, []))

        nexus_path = os.path.join(NEXUS_DIR, nx_file)
        output_path = os.path.join(OUTPUT_NEXUS, f"{locus}_filtered.nexus")

        nx = Nexus.Nexus()
        nx.read(nexus_path)

        all_taxa_in_file = list(nx.matrix.keys())
        removed = 0
        for taxon in list(nx.matrix):
            if taxon in to_discard:
                del nx.matrix[taxon]
                removed += 1

        kept = len(nx.matrix)
        log.write(f"{locus}\t{len(all_taxa_in_file)}\t{kept}\t{removed}\n")

        with open(output_path, "w") as out_nx:
            out_nx.write(nx.write_nexus_data())

# ===============================
# FILTER TREES
# ===============================
print("Filtering trees...")
for tree_file in tree_files:
    locus = tree_file.split(".")[1]
    path = os.path.join(TREE_DIR, tree_file)
    tree = Phylo.read(path, "newick")
    to_discard = set(discarded_by_locus.get(locus, []))
    tree.prune_targets = [t for t in tree.get_terminals() if t.name in to_discard]
    for node in tree.prune_targets:
        tree.prune(node)
    out_path = os.path.join(OUTPUT_TREES, f"{tree_file}_filtered")
    Phylo.write(tree, out_path, "newick")

# ===============================
# SUMMARY
# ===============================
print("Generating summary...")
total_loci = len(discarded_by_locus)
total_discards = len(all_discards)
discard_counts = Counter(all_discards)

with open(SUMMARY_LOG, "w") as f:
    f.write(f"Total loci with discards: {total_loci}\n")
    f.write(f"Total individual discards: {total_discards}\n\n")
    f.write("Taxa discarded in total:\n")
    f.write("Taxon\t# discards\n")
    for taxon, count in discard_counts.most_common():
        f.write(f"{taxon}\t{count}\n")

print("Process completed.")
