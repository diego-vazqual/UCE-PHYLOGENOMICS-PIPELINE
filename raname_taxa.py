import os
import re

# Step 1: Read the taxon map
taxa_map = {}
with open("taxa_map.txt", "r") as f:
    for line in f:
        if line.strip():
            parts = line.strip().split()
            if len(parts) == 2:
                name, taxa = parts
                taxa_map[taxa] = name

# Step 2: Process each file of type RAxML_bipartitions*
folder = "raxml_out/bestTree"

for filename in os.listdir(folder):
    if filename.startswith("RAxML_bestTree."):
        filepath = os.path.join(folder, filename)

        # Read the tree
        with open(filepath, "r") as f:
            tree = f.read()

        # Step 3: Replace each TaxaX with the real name
        for taxa_code, real_name in taxa_map.items():
            # Ensures that only whole taxa are replaced
            tree = re.sub(r'\b{}\b'.format(re.escape(taxa_code)), real_name, tree)

        # Step 4: Overwrite the file with the real names
        with open(filepath, "w") as f:
            f.write(tree)

print("✅ All trees have been updated with the real names.")
