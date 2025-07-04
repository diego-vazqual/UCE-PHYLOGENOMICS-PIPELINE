import os
from Bio.Nexus import Nexus

# Paths
input_dir = "mafft-nexus-internal-no-trimmed-gblocks-clean"
output_dir_phylip = "mafft-phylip-nexus-internal-no-trimmed-gblocks-clean"
mapping_file = "taxa_map.txt"
log_file = "rename_log.txt"

# Create output directory if it doesn't exist
os.makedirs(output_dir_phylip, exist_ok=True)

# Read the mapping file
species_to_taxa = {}
with open(mapping_file, 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) == 2:
            species, taxa = parts
            species_to_taxa[species.replace("-", "_")] = taxa

# Initialize error log
errors = []

# Process NEXUS files
for filename in os.listdir(input_dir):
    if filename.endswith(".nexus"):
        input_path = os.path.join(input_dir, filename)
        output_path_phylip = os.path.join(output_dir_phylip, filename.replace(".nexus", ".phylip"))

        nexus = Nexus.Nexus()
        nexus.read(input_path)

        new_matrix = {}
        for taxon, sequence in nexus.matrix.items():
            taxon_clean = taxon.replace("-", "_")
            new_name = species_to_taxa.get(taxon_clean)
            if new_name:
                new_matrix[new_name] = sequence
            else:
                errors.append(f"{filename}: '{taxon}' not found in the mapping file")
                new_matrix[taxon] = sequence  # Use original name if not found in mapping

        # Write PHYLIP file
        with open(output_path_phylip, "w") as out_f:
            ntax = len(new_matrix)
            nchar = len(next(iter(new_matrix.values())))
            out_f.write(f"{ntax} {nchar}\n")
            for taxon, seq in new_matrix.items():
                short_name = taxon[:10]  # Strict PHYLIP allows max. 10 characters for taxon name
                out_f.write(f"{short_name:<10} {seq}\n")

        print(f"Processed: {filename}")

# Write error log
if errors:
    with open(log_file, 'w') as log:
        for err in errors:
            log.write(err + "\n")
    print(f"\nFound {len(errors)} errors. See '{log_file}' for details.")
else:
    print("\nAll files were processed without errors.")
