# Guide for UCE Phylogenomic Analysis


## INDEX

1. [Requirements](#1-requirements)  
2. [Clean reads](#2-clean-reads)  
3. [Contig assembly](#3-contig-assembly)  
4. [Target enrichment](#4-target-enrichment)  
5. [Alignment](#5-alignment)  
   - [5.1 Aligning UCE loci](#51-aligning-uce-loci)  
   - [5.2 Gblocks](#52-gblocks)  
   - [5.3 Alignment cleaning](#53-alignment-cleaning)
   - [5.4 Remove empty sequences](#54-remove-empty-sequences)  
   - [5.5 Remove long terminal branches](#54-remove-long-terminal-branches)  
   - [5.6 Preparing data for phylogenomic analysis](#55-preparing-data-for-phylogenomic-analysis)  
   - [5.7 Final data matrices](#56-final-data-matrices)  
6. [Phylogenomic Analysis](#6-phylogenomic-analysis)  
   - [6.1 IQ-TREE2](#61-iq-tree2)  
   - [6.2 ExaBayes](#62-exabayes)  
   - [6.3 ASTRAL-III](#63-astral-iii)  
7. [References](#7-references)
---

## 1. REQUIREMENTS
* ASTRAL-III ([https://sites.google.com/view/cd-hit](https://github.com/smirarab/ASTRAL))
* CD-HIT-DUP (https://sites.google.com/view/cd-hit)
* ExaBayes v1.5 (https://cme.h-its.org/exelixis/web/software/exabayes/)
* Fastp (https://github.com/OpenGene/fastp)
* IQ-TREE v2.2.5 (https://iqtree.github.io/)
* OpenMPI 4.1.0-10 (https://www.open-mpi.org)
* Phyluce v1.7.2 (https://phyluce.readthedocs.io/en/latest/)
---

## 2. Clean reads

### 2.1 Counting the read data
Before cleaning, check that the R1 and R2 files of each sample contain the same number of reads. The command below expects one folder per sample containing its R1 and R2 files, and prints the sample name followed by the number of reads in R1 and in R2.
```
for r1 in */*_R1*.fastq.gz; do
    r2=${r1/_R1/_R2}
    echo -e "$(dirname $r1)\t$(( $(zcat $r1 | wc -l) / 4 ))\t$(( $(zcat $r2 | wc -l) / 4 ))"
done
```
⚠️ Important: If R1 and R2 have a different number of reads, the files are probably truncated or corrupted. Download or copy them again before continuing.

### 2.2 FASTP
Quality filtering (removal of low-quality bases and reads), adapter trimming, trimming of poly-G tails and base correction in overlapping read pairs are performed with fastp, using the script `fastp.sh`. The script takes as input a folder with one subfolder per sample, each containing its R1 and R2 files (`*_R1*.fastq.gz`, `*_R2*.fastq.gz`), and writes the cleaned reads (`<sample>-READ1.fastq`, `<sample>-READ2.fastq`) and the fastp quality reports (HTML and JSON) to one subfolder per sample in the output directory.

Activate the environment where fastp is installed before running the script:
```
conda activate phyluce-1.7.2
```
```
bash fastp.sh data/fastq data/clean_reads_fastq
```
☁️ Check the HTML report of each sample to confirm that adapters have been removed and that read quality is acceptable after filtering.

### 2.3 CD-HIT-DUP

Duplicate read pairs are removed with CD-HIT-DUP before assembling the contigs with SPAdes. The script cd-hit-dup.sh takes as input the output folder of fastp.sh and writes the deduplicated reads, compressed (`<sample>-READ1.fastq`, `<sample>-READ2.fastq`), to one subfolder per sample in the output directory. Two read pairs are considered duplicates when their first 30 bases are identical (-u 30).
```
bash cd-hit-dup.sh data/clean_reads_fastq data/clean_reads_cdhitdup
```

## 3. Contig assembly

### 3.1 SPAdes
Clean reads are assembled into contigs with SPAdes, as implemented in Phyluce. Phyluce requires a configuration file (`assembly.conf`) with a `[samples]` header followed by one line per sample with the sample name and the absolute path to the folder containing its clean reads, separated by a colon.

From the working directory, create the configuration file from the output folder of section 2.3. First, write the header:
```:
```
echo "[samples]" > assembly.conf
```
Then, add one line per sample with its name and the path to its clean reads:
```
for i in data/clean_reads_cdhitdup/*/; do
    s=$(basename "$i")
    echo "$s:$(pwd)/data/clean_reads_cdhitdup/$s"
done >> assembly.conf
```
✏️ Example `assembly.conf`:
```
[samples]
Acteon_sp:/home/user/work_directory/data/clean_reads_cdhitdup/Acteon_sp
Acteon_tornatilis:/home/user/work_directory/data/clean_reads_cdhitdup/Acteon_tornatilis
Akera_bullata:/home/user/work_directory/data/clean_reads_cdhitdup/Akera_bullata
```
⚠️ Important:
- Phyluce looks for read files ending in `-READ1.fastq.gz` and `-READ2.fastq.gz`. The scripts in section 2 already name the files this way.
- Sample names should not contain spaces or hyphens (`-`). These names are carried to the contigs, alignments and trees, and can cause errors in later steps.

Activate the environment where Phyluce is installed and run the assembly from the working directory (the folder containing `assembly.conf`):
```
conda activate phyluce-1.7.2
```
```
phyluce_assembly_assemblo_spades \
    --conf assembly.conf \
    --output data/spades_assemblies \
    --cores 24
```
📌 Adjust `--cores` to the available resources and the number of samples.

The assembled contigs of each sample are written to `data/spades_assemblies/<sample>.contigs.fasta`. This `contigs/` folder is the input for the next step.

☁️ To check the quality of the assemblies (number of contigs, total length, mean, median and maximum contig length per sample), first create the header of the statistics file:
```
echo "Sample,contigs,total bp,mean length,95 CI length,min length,max length,median length,contigs >1kb" > assembly_stats.csv
```
And calculate the statistics for each sample:
```
for i in data/spades_assemblies/*.contigs.fasta; do
    phyluce_assembly_get_fasta_lengths --input "$i" --csv >> assembly_stats.csv
done
```

## 4. Target enrichment

### 4.1 Finding UCE loci
After assembling the contigs, the next step is to identify which contigs correspond to UCE loci and to exclude non-target sequences. To do this, the contigs of each sample are compared with the probe set used for target enrichment (in this example, `Probeset-70nt.fasta`).

From this point on, all commands are run from the working directory, which contains the probe file and the `data/` folder. Keeping a clean and consistent folder structure is especially important from here on to avoid errors.

```
phyluce_assembly_match_contigs_to_probes \
    --contigs data/spades_assemblies \
    --probes Probeset-70nt.fasta \
    --output uce-search-results \
    --keep-duplicates duplicates.txt \
    --csv uce_search_results.csv
```

By default, `phyluce_assembly_match_contigs_to_probes` removes two types of duplicates from each sample:
- **Probes hitting multiple contigs:** the probes of a single UCE locus match more than one contig.
- **Contigs hitting multiple probes:** a single contig is matched by probes designed for different UCE loci.

Both cases suggest that the locus is problematic, either because of assembly errors or because of biological events such as gene duplication (potential paralogs), so the affected loci and contigs are discarded for that sample. The option `--keep-duplicates` saves the list of duplicates detected in each sample to `duplicates.txt`, which is used in section 4.2.2 to recover part of these loci.

☁️ The search results (number of contigs, number of UCE loci recovered and number of duplicates removed per sample) are printed to the terminal and, with the option `--csv`, also saved to `uce_search_results.csv`.

### 4.2 Extracting UCE loci
Next, define the set of taxa to be included in the analysis in a configuration file (`taxon-set.conf`). The first line contains the name of the taxon set in square brackets (e.g. `[taxon_set1]`), followed by the names of the taxa, one per line.

⚠️ Important: The taxon names must exactly match the sample names used in the assembly. To include all assembled samples and avoid typing errors, create the file directly from the contig files. First, write the name of the taxon set:
```
echo "[taxon_set1]" > taxon-set.conf
```
Then, add the names of all assembled samples:
```
for f in data/spades_assemblies/*.contigs.fasta; do basename "$f" .contigs.fasta; done >> taxon-set.conf
```
To exclude some samples, delete their lines from the file (e.g. with `nano`).

✏️ Example `taxon-set.conf`:
```
[taxon_set1]
Acteon_sp
Acteon_tornatilis
Akera_bullata
```
💡 Several taxon sets can be defined in the same file, each under its own header (e.g. `[taxon_set2]`), and selected with the option `--taxon-group`.

Then, generate the list of UCE loci recovered for each taxon of the set:
```
phyluce_assembly_get_match_counts \
    --locus-db uce-search-results/probe.matches.sqlite \
    --taxon-list-config taxon-set.conf \
    --taxon-group 'taxon_set1' \
    --incomplete-matrix \
    --output taxon_set1-taxa-incomplete.conf
```
And extract the corresponding sequences in FASTA format:
```
phyluce_assembly_get_fastas_from_match_counts \
    --contigs data/spades_assemblies \
    --locus-db uce-search-results/probe.matches.sqlite \
    --match-count-output taxon_set1-taxa-incomplete.conf \
    --output taxon_set1-taxa-incomplete.fasta \
    --incomplete-matrix taxon_set1-taxa-incomplete.incomplete \
    --log-path log
```
The option `--incomplete-matrix` keeps loci that are missing in some taxa. Loci are filtered by occupancy later, when building the final data matrices (section 5.7).

#### 4.2.1 Statistics per taxon
☁️ To obtain the number and length of the UCE loci recovered for each taxon, first split `taxon_set1-taxa-incomplete.fasta` into one FASTA file per taxon:
```
phyluce_assembly_explode_get_fastas_file \
    --input taxon_set1-taxa-incomplete.fasta \
    --output exploded-fastas \
    --by-taxon
```
Create the header of the statistics file:
```
echo "Sample,UCE loci,total bp,mean length,95 CI length,min length,max length,median length,contigs >1kb" > exploded_fastas.csv
```
And calculate the statistics for each taxon:
```
for i in exploded-fastas/*.fasta; do
    phyluce_assembly_get_fasta_lengths --input "$i" --csv >> exploded_fastas.csv
done
```
✏️ Example `exploded_fastas.csv`:
```
Sample,UCE loci,total bp,mean length,95 CI length,min length,max length,median length,contigs >1kb
Acteon-sp.unaligned.fasta,1395,756820,542.5,3.8,208,1245,560.0,12
Acteon-tornatilis.unaligned.fasta,1455,776310,533.5,4.6,201,1189,548.0,9
```

#### 4.2.2 Recovering duplicate UCE loci
The loci removed as duplicates in section 4.1 are lost for the affected taxa. To recover part of them, we use the script `phyluce_assembly_parse_duplicates_file.py`, written by Brant C. Faircloth for PHYLUCE (Faircloth, 2016). The script extracts from the assemblies the contigs of the loci listed under **probes hitting multiple contigs** in `duplicates.txt`. Contigs hitting multiple probes are not recovered.

First, create a folder for the recovery of duplicates:
```
mkdir -p duplicates
```
Move the list of duplicates generated in section 4.1 into this folder:
```
mv duplicates.txt duplicates/
```
Run the script `phyluce_assembly_parse_duplicates_file.py` within the Phyluce environment to extract the duplicate contigs in FASTA format:
```
python3 phyluce_assembly_parse_duplicates_file.py \
    --contigs data/spades_assemblies \
    --duplicates-file duplicates/duplicates.txt \
    --output duplicates/duplicates.fasta \
    --log-path log
```
For each duplicated locus and taxon, the script labels the contigs as `_DUPE1`, `_DUPE2`, etc., following their order in the assembly file. Because SPAdes writes the contigs sorted by decreasing length, `DUPE1` is always the longest contig of each duplicated locus. We keep only this contig, as it is the most complete copy of the locus: it covers the UCE core region targeted by the probes and the largest part of the more variable flanking regions, which contain most of the phylogenetically informative sites.

Keep only the `DUPE1` contigs and remove the `_DUPE1` label, so that the headers have the same format as the rest of the loci:
```
awk '/^>/ {keep = ($1 ~ /_DUPE1$/)} keep' duplicates/duplicates.fasta | sed 's/_DUPE1//' > duplicates/duplicates_DUPE1.fasta
```
☁️ Count the number of recovered loci (all taxa):
```
grep -c '^>' duplicates/duplicates_DUPE1.fasta
```
⚠️ Important:
- `duplicates.txt` includes all the samples searched in section 4.1. If `taxon-set.conf` does not include all of them, remove the sequences of the excluded taxa from `duplicates_DUPE1.fasta` before the next step.
- A recovered locus may still be a paralog. Sequences that behave as paralogs usually produce abnormally long terminal branches in the gene trees, and are removed by the long-branch filter (section 5.5).

Finally, add the recovered loci to the rest of the UCE loci:
```
cat taxon_set1-taxa-incomplete.fasta duplicates/duplicates_DUPE1.fasta > final-taxon_set1-taxa-incomplete.fasta
```

☁️ To check how many loci have been recovered per taxon, repeat the statistics of section 4.2.1 on the final file and compare them with `exploded_fastas.csv`. First, split the final file by taxon:
```
phyluce_assembly_explode_get_fastas_file \
    --input final-taxon_set1-taxa-incomplete.fasta \
    --output exploded-fastas-final \
    --by-taxon
```
Create the header of the statistics file:
```
echo "Sample,UCE loci,total bp,mean length,95 CI length,min length,max length,median length,contigs >1kb" > exploded_fastas_final.csv
```
And calculate the statistics for each taxon:
```
for i in exploded-fastas-final/*.fasta; do
    phyluce_assembly_get_fasta_lengths --input "$i" --csv >> exploded_fastas_final.csv
done
```
## 5. Alignment
### 5.1 Aligning UCE loci
The recovered UCE loci were aligned using the MAFFT algorithm. For the alignment, we specified the options `--no-trim`, which disables the automatic trimming of unaligned regions and retains all positions in the alignments, and `--taxa`, where we must indicate the number of taxa in our dataset.

```
phyluce_align_seqcap_align \
    --input final-taxon_set1-taxa-incomplete.fasta \
    --output mafft-nexus-no-internal-trimmed \
    --taxa 117 \
    --aligner mafft \
    --cores 24 \
    --incomplete-matrix \
    --output-format fasta \
    --no-trim \
    --log-path log
```

### 5.2 Gblocks
To mask regions of sequences that are highly variable, ambiguous, or error-prone, Gblocks is used as a bioinformatics masking tool. Gblocks identifies and removes unreliable or poorly aligned sections in DNA, RNA, or protein sequences, retaining only the regions that are well and consistently aligned. Therefore, this tool improves the quality of multiple sequence alignments in phylogenetic studies.

Phyluce includes a specific command to run Gblocks directly:
```
phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
        --alignments mafft-nexus-no-internal-trimmed  \
        --output mafft-nexus-internal-no-trimmed-gblocks \
        --b1 0.5 \
        --b2 0.85 \
        --b3 4 \
        --b4 8 \
        --cores 24 \
        --log-path log
```
There are three possible Gblocks configurations, which differ in the values of the parameters `--b1`, `--b2`, `--b3`, and `--b4`. The choice of configuration depends on the nature of the data and the goals of the analysis:

**Configuration 1**: --b1 0.5 --b2 0.85 --b3 4 --b4 8  # very restrictive

**Configuration 2**: --b1 0.5 --b2 0.5 --b3 6 --b4 6  # intermediate

**Configuration 3**: --b1 0.5 --b2 0.5 --b3 10 --b4 4  # very conservative

### 5.3 Alignment cleaning
Once the UCE loci have been aligned and the multiple alignments masked, each file contains both the taxon name and the locus name. For downstream analyses, it is necessary to remove the locus name from the files so that they include only the taxon name.

To remove the locus name, use the following command:
```
phyluce_align_remove_locus_name_from_files \
    --alignments mafft-nexus-internal-no-trimmed-gblocks \
    --output mafft-nexus-internal-no-trimmed-gblocks-clean \
    --cores 24 \
    --log-path log
```

### 5.4 Remove empty sequences
After Gblocks, some alignments may contain sequences made only of gaps or missing data (`-`, `?`, `N`). This happens when a taxon only had sequence in the columns removed by Gblocks. These empty sequences cause two problems:

• RAxML refuses to build a gene tree from an alignment that contains them (`ERROR: Sequence ... consists entirely of undetermined values`). As a result, those loci have no gene tree for the long-branch filter (section 5.5) or for ASTRAL (section 7.3), and pass the long-branch filter without being checked.

• Phyluce counts them as present taxa when building the occupancy matrices (section 5.7), so the occupancy of those loci is overestimated.

To remove them, we use the script `remove_empty_sequences.py`, which takes as input the folder with the cleaned alignments in nexus format and writes the alignments without empty sequences, with the same file names, to a new folder:
```
python3 remove_empty_sequences.py \
    mafft-nexus-internal-no-trimmed-gblocks-clean \
    mafft-nexus-internal-no-trimmed-gblocks-clean-noempty
```
☁️ The script reports how many alignments contained empty sequences and how many sequences were removed, and saves the details (locus, number of taxa before and after, and removed taxa) in `empty_sequences_removed.tsv`.

⚠️ Important: From this step on, use the folder `mafft-nexus-internal-no-trimmed-gblocks-clean-noempty` as input for the long-branch filter.

### 5.5 Remove long terminal branches
Even after the initial cleaning steps, alignments may still contain atypical sequences, which could represent chimeras, contaminants, or non-homologous alignments. These sequences can introduce artifacts in subsequent phylogenetic inference, particularly by generating abnormally long terminal branches, which may lead to the phenomenon known as long branch attraction (LBA). This occurs when highly divergent sequences are incorrectly grouped due to elevated substitution rates or high sequence dissimilarity.

To mitigate this issue, we identified and removed abnormally long terminal branches using the following steps.

To identify and remove abnormally long terminal branches, we first build gene trees from each alignment. For this step, alignments must be in  `.phylip` format.

⚠️ Important: The PHYLIP format requires taxon names to be no longer than 10 characters. If names exceed this limit, programs may throw read errors, truncate names, or produce incorrect results.

✏️To avoid these issues, taxon names must be renamed using a mapping file `taxa_map.txt` like the example below:
```
Acteon_sp    Taxa1
Acteon_tornatilis    Taxa2
Akera_bullata    Taxa3
Ammonicera_sp    Taxa4
```
In this file `taxa_map.txt`, the first column contains the original taxon names (as they appear in the alignments), and the second column, separated by a tab, contains the new shortened names. In this example, we use generic identifiers like Taxa1, Taxa2, etc.

The script `alignment_phylip_format.py` renames taxon names and converts alignments from nexus format to phylip format. 
This script takes as input the folder containing the alignments in nexus format (`mafft-nexus-internal-no-trimmed-gblocks-clean-noempty`) and the mapping file with the original and new taxon names (`taxa_map.txt`). The output is saved in a directory called `mafft-phylip-nexus-internal-no-trimmed-gblocks-clean`.

An example of running the script would be:
```
python3 alignment_phylip_format.py
```

Once the alignments have been converted to phylip format, we use the script `single_gene_trees.sh` to reconstruct the individual gene trees. This is done using RAxML under the GTRGAMMA model with 100 bootstrap replicates.
The script takes as input the alignments in phylip format `mafft-phylip-nexus-internal-no-trimmed-gblocks-clean` and generates the output files in a directory named `gene_trees`.

An example of running the script would be:
```
bash single_gene_trees.sh
```

Once the script has finished, we organize the files generated inside the `gene_trees` directory using the following commands:
```
cd gene_trees
```
```
mkdir RAxML_bestTree
mv RAxML_bestTree.uce* RAxML_bestTree
```
```
mkdir RAxML_bipartitions
mv RAxML_bipartitions.uce* RAxML_bipartitions
```
```
mkdir RAxML_bipartitionsBranchLabels
mv RAxML_bipartitionsBranchLabels.uce* RAxML_bipartitionsBranchLabels
```
```
mkdir RAxML_bootstrap
mv RAxML_bootstrap.uce* RAxML_bootstrap
```
```
mkdir RAxML_info
mv RAxML_info.uce* RAxML_info
```

Once everything is organized, we will rename the tips of the individual gene trees using their original names, but only those located in the `genes_tree/RAxML_bipartitions` directory, as these are the ones that will be used for long branch filtering. To do this, we run the `rename_taxa.py` script, which takes as input the gene trees for each locus from the `genes_tree/RAxML_bipartitions` folder and a mapping file called `taxa_map.txt`.

An example of how to run the script from the `work_directory/` is:
```
python3 rename_taxa.py
```
To detect and remove abnormally long terminal branches, we use the `remove_long_branches.py` script, which applies a dual filtering approach:

• Branch length distribution per locus: branch lengths are analyzed within each individual gene tree.

• Branch length distribution per taxon: branch lengths are analyzed across all gene trees where the taxon is present.

Each taxon is first assigned to a group (e.g., at the family level) using a predefined file `spms_info.txt`. Based on this classification, taxa are labeled in each tree as singletons (the only representative of their group in that tree) or non-singletons (present along with at least one other member of the same group).

Branch length distributions are calculated separately for both categories. A terminal branch is removed if its length exceeds the mean + 2.7 standard deviations of both distributions: at the locus level and the taxon level.

The threshold value used here is 2.7 (Fedosov et al., 2024), but it can be adjusted depending on the nature of the data and the desired level of stringency (higher values are more stringent; lower values are less so). If you wish to adjust these thresholds, you must modify the corresponding values in the **FILTERING PARAMETERS** section of the `remove_long_branches.py` script.

✏️ Example `spms_info.txt`:
```
Bigtaxon    Taxon
Acteonidae    Acteon_sp    
Acteonidae    Acteon_tornatilis
Aplysiida    Akera_bullata
Ammonicera    Ammonicera_sp 
```
In the `spms_info.txt` file, the first column contains the assigned taxonomic group, headed by `Bigtaxon`, and the second column, headed by `Taxon` and separated by a tab, includes the original names of the taxa or terminals. For the script to run correctly, this file must follow that structure.

The script takes as input the individual gene trees from the `genes_tree/RAxML_bipartitions` directory and the alignments in nexus format (`mafft-nexus-internal-no-trimmed-gblocks-clean-noempty`). As output, it generates the filtered individual gene trees in the `clean_genes_trees` directory and **directly modifies the original alignments**.

⚠️ Important: If you wish to preserve the original alignments before filtering atypical sequences, make a copy, as the script directly overwrites the original files.

☁️ Additionally, it generates detailed statistics, including: how many terminals were removed per taxon (`summary_report.txt`), in which loci those removals occurred (`long_branch_discards.tsv`), and how many terminals remained per alignment in each locus (`filtered_nexus_log.tsv`).

An example of how to run the script from the `work_directory/` is:
```
python3 remove_long_branches.py
```

☁️ To determine how many loci have been removed per taxon, we can recalculate the statistics from steps `4.2.1` and `4.2.2` and compare how many UCE loci have been recovered.
 
 First we convert the alignments to fasta format:
  ```
phyluce_align_convert_one_align_to_another \
--alignments mafft-nexus-internal-trimmed-gblocks-clean \
--output mafft-fastas-internal-trimmed-gblocks-clean \
--input-format nexus \
--output-format fasta \
--cores 24 \
--log-path log

```
Next, we will tag each file with the name of the corresponding UCE using the `add_tag.sh` script.
 ```
bash add_tag.sh
```
Once all the UCEs are labeled, we need to concatenate all the files into a single monolithic file, similar to the process done previously.
```
cd mafft-fastas-internal-trimmed-gblocks-clean
cat * >>  all-fastas-clean
 ```
 ```
phyluce_assembly_explode_get_fastas_file \
--input all-fastas-clean \
--output exploded-fastas-clean \
--by-taxon
```
```
echo "Sample ID,UCE loci,total bp,mean length,95 CI length,min length,max length, median legnth, contigs >1kb" > fasta_lengths_clean.csv
for i in exploded-fastas-clean/*.fasta; do  
    phyluce_assembly_get_fasta_lengths --input "$i" --csv >> fasta_lengths_clean.csv  
done
```

### 5.6 Final data matrices
At this stage, we are interested in minimizing noise in the data and increasing confidence in phylogenetic inferences by removing loci or genes that may be less informative or more prone to error, as well as eliminating missing data. To achieve this, we generate occupancy matrices with different thresholds.

☁️ Using the Phyluce function `phyluce_align_get_align_summary_data`, we can first obtain information about the number of loci in matrices with varying levels of completeness, the total number of available loci, among other useful metrics.

To obtain this summary, use the following command:
```
phyluce_align_get_align_summary_data \
    --alignments mafft-nexus-internal-no-trimmed-gblocks-clean \
    --cores 24 \
    --log-path log 
```

A 50% occupancy matrix means that each column of the matrix (i.e., each locus) must be present in at least 50% of the taxa to be included. This helps reduce noise and improve data quality for phylogenetic analysis.

For example, to generate a 50% occupancy matrix, use the following command:
```
phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-nexus-internal-no-trimmed-gblocks-clean \
    --taxa 117 \
    --percent 0.50 \
    --output mafft-clean-nexus-internal-trimmed-gblocks-clean-50p \
    --cores 24 \
    --log-path log
```
In the `--taxa` option, you must specify the number of taxa present, and in the `--percent` option, the desired occupancy threshold for the matrix (0.5 = 50%; 0.75 = 75%...).

☁️ To obtain statistics on how many ECUs have been conserved for each species, we can recalculate the statistics in a similar way as above.

 ```
phyluce_align_convert_one_align_to_another \
--alignments mafft-nexus-internal-trimmed-gblocks-clean-50p \
--output mafft-fastas-internal-trimmed-gblocks-clean-50p \
--input-format nexus \
--output-format fasta \
--cores 24 \
--log-path log
```
 ```
bash add_tag.sh
```
 ```
cd mafft-clean-fastas-internal-trimmed-gblocks-clean-50p
cat * >>  all-fastas-50p
```
 ```
phyluce_assembly_explode_get_fastas_file \
--input all-fastas-50p \
--output exploded-fastas-50p \
--by-taxon
```
 ```
echo "Sample ID,UCE loci,total bp,mean length,95 CI length,min length,max length, median legnth, contigs >1kb" > fasta_lengths_50p.csv
for i in exploded-fastas-50p/*.fasta; do  
    phyluce_assembly_get_fasta_lengths --input "$i" --csv >> fasta_lengths_50p.csv  
done
```
### 5.7 Preparing data for phylogenomic analysis
To ensure that IQ-TREE2, ExaBayes, and ASTRAL can recognize our files, we need to concatenate all our UCEs from the matrix into a single file in phylip format.

To do this, run:
 ```
phyluce_align_concatenate_alignments \
    --alignments mafft-clean-nexus-internal-trimmed-gblocks-clean-50p \
    --output mafft-clean-nexus-internal-trimmed-gblocks-clean-50p-IQTree \
    --phylip \
    --log-path log
```

## 6. Phylogenomic Analysis
### 6.1 IQ-TREE2
To infer the maximum likelihood tree, we will use IQ-TREE2. In our case, we will use the GHOST model (GTR+FO*H4) and estimate node support values using the Ultrafast Bootstrap approximation with 1,500 replicates.
These parameters can be adjusted to suit your data by modifying the `-m` (model) and `-B` (number of replicates) options.
 ```
iqtree2 --seqtype DNA --ninit 10 -B 1500 -s mafft-clean-nexus-internal-trimmed-gblocks-clean-50p-IQTree/mafft-clean-nexus-internal-trimmed-gblocks-clean-50p-IQTree.phylip --prefix iqtree-GHOST-50p -m GTR+FO*H4 -T 24 --rcluster 10 --mrate G,R,E
```
### 6.2 Exabayes
To infer the tree using Bayesian statistics, we will use the ExaBayes program. 

📌 This software requires a considerable amount of memory and is a slow process. Therefore, a compiler is used to run parallel processes (MPI), and in our case, the analysis was performed in a high-performance computing (HPC) environment managed with Slurm. We also recommend using tmux or another background process manager, as this is a time-consuming task.

To automate the process, we created a script called `exabayes.sh`, which runs four successive executions. To run this analysis using the script, you will need a configuration file `config.nex` that specifies the parameters and variables for the Bayesian analysis. Be sure to review and adjust this file according to the requirements of your study.

Once the execution is complete, we will check the parameters to confirm that it was carried out correctly.

• Check if parameters converged (ESS should be >100, PSRF should be <1.1)
```
postProcParam -f ExaBayes_parameters.run-0.run1 ExaBayes_parameters.run-0.run2 ExaBayes_parameters.run-0.run3 ExaBayes_parameters.run-0.run4 -n combinedParams
```
• Calculate standard deviation of split frequencies (aka convergence, ASDSF should be <1%)
```
sdsf -f ExaBayes_topologies.run-0.run1 ExaBayes_topologies.run-0.run2 ExaBayes_topologies.run-0.run3 ExaBayes_topologies.run-0.run4
```
• Verify that the chains have explored multiple topologies and that there is consistency among replicates (Number of topologies should be > 1,000).
```
grep -v "^#" ExaBayes_topologies.run-0.run1 | sort | uniq | wc -l
grep -v "^#" ExaBayes_topologies.run-0.run2 | sort | uniq | wc -l
grep -v "^#" ExaBayes_topologies.run-0.run3 | sort | uniq | wc -l
grep -v "^#" ExaBayes_topologies.run-0.run4 | sort | uniq | wc -l
```

Finally, the consensus tree must be generated:
```
consense -f ExaBayes_topologies.run-0.run1 ExaBayes_topologies.run-0.run2 ExaBayes_topologies.run-0.run3 ExaBayes_topologies.run-0.run4
```

### 6.3 ASTRAL-III
Finally, to infer the coalescent trees, we will utilize ASTRAL-III, implemented in Java.

To initiate the process, we first need to employ RAxML to generate a single gene trees. Just like in section `5.4`, we first need to rename the taxon names and convert the alignments from Nexus format to Phylip format, but in this case for the alignments of the 50% occupancy matrix `mafft-fastas-internal-trimmed-gblocks-clean-50p`.

To do this, we run the script `alignment_phylip_format.py` on the alignment `mafft-fastas-internal-trimmed-gblocks-clean-50p`.
Make sure the `taxa_map.txt` file is available in order to truncate the taxon names.

Next, we reconstruct the individual gene trees using RAxML (-N 10 -m GTRGAMMA) with the script `astral_gene_trees.sh`. The script takes as input the phylip format alignments from the 50% occupancy matrix and generates an output containing the individual gene trees in a directory called `astral_raxml_output`.
An example of how to run the script from the `work_directory/` is:
```
bash astral_gene_trees.sh 
```

Then, we rename the individual trees with their original taxon names using the `rename_taxa.py` script, just as in section `5.4`, but this time applied to the trees located in the `astral_raxml_output/bestTree` directory.

Finally, we run the following command to obtain the final ASTRAL-III tree:
```
java -jar /home/intern/Desktop/apps/ASTRAL/astral.5.7.8.jar -i astral_raxml_output/bestTree/*.genetree -o astral_sptree.treefile
```

## 7. References
• Aberer, A. J., Kobert, K., & Stamatakis, A. (2014). ExaBayes: massively parallel Bayesian tree inference for the whole-genome era. Molecular Biology and Evolution, 31(10), 2553–2556. https://doi.org/10.1093/molbev/msu236

• Bankevich, A., Nurk, S., Antipov, D., Gurevich, A. A., Dvorkin, M., Kulikov, A. S., Lesin, V. M., Nikolenko, S. I., Pham, S., Prjibelski, A. D., Pyshkin, A. V., Sirotkin, A. V., Vyahhi, N., Tesler, G., Alekseyev, M. A., & Pevzner, P. A. (2012). SPAdes: A new genome assembly algorithm and its applications to single-cell sequencing. Journal of Computational Biology, 19(5), 455–477. https://doi.org/10.1089/cmb.2012.0021

• Castresana, J. (2000). Selection of conserved blocks from multiple alignments for their use in phylogenetic analysis. Molecular Biology and Evolution, 17(4), 540–552. https://academic.oup.com/mbe/article/17/4/540/1127654

• Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884-i890. https://doi.org/10.1093/bioinformatics/bty560

• Crotty, S. M., Minh, B. Q., Bean, N. G., Holland, B. R., Tuke, J., Jermiin, L. S., & Haeseler, A. V. (2020). GHOST: recovering historical signal from heterotachously evolved sequence alignments. Systematic Biology, 69(2), 249–264. https://doi.org/10.1093/sysbio/syz051

• Cunha, T. J., & Giribet, G. (2019). A congruent topology for deep gastropod relationships. Proceedings of the Royal Society B: Biological Sciences, 286(1898), 20182776. https://doi.org/10.1098/rspb.2018.2776

• Faircloth, B. C. (2016). PHYLUCE is a software package for the analysis of conserved genomic loci. Bioinformatics, 32(5), 786–788. https://doi.org/10.1093/bioinformatics/btv646

• Fedosov, A. E., Zaharias, P., Lemarcis, T., Modica, M. V., Holford, M., Oliverio, M., Kantor, Y. I. & Pullandre, N. (2024). Phylogenomics of Neogastropoda: the backbone hidden in the bush. Systematic Biology, 73(3), 521–531. https://doi.org/10.1093/sysbio/syae010

• Fu, L., Niu, B., Zhu, Z., Wu, S., & Li, W. (2012). CD-HIT: Accelerated for clustering the next-generation sequencing data. Bioinformatics, 28(23), 3150–3152. https://doi.org/10.1093/bioinformatics/bts565

• Hoang, D. T., Chernomor, O., Von Haeseler, A., Minh, B. Q., & Vinh, L. S. (2018). UFBoot2: improving the ultrafast bootstrap approximation. Molecular Biology and Evolution, 35(2), 518–522. https://doi.org/10.1093/molbev/msx281

• Katoh, K., & Standley, D. M. (2013). MAFFT Multiple Sequence Alignment software version 7: Improvements in performance and usability. Molecular Biology and Evolution, 30(4), 772–780. https://doi.org/10.1093/molbev/mst010

• Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., Von Haeseler, A., & Lanfear, R. (2020). IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. Molecular Biology and Evolution, 37(5), 1530–1534. https://doi.org/10.1093/molbev/msaa015

• Mirarab, S., Reaz, R., Bayzid, M. S., Zimmermann, T., Swenson, M. S., & Warnow, T. (2014). ASTRAL: genome-scale coalescent-based species tree estimation. Bioinformatics, 30(17), i541–i548. https://doi.org/10.1093/bioinformatics/btu462

• Sayyari, E., & Mirarab, S. (2016). Fast coalescent-based computation of local branch support from quartet frequencies. Molecular Biology and Evolution, 33(7), 1654–1668. https://doi.org/10.1093/molbev/msw079

• Stamatakis, A. (2006). RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models. Bioinformatics, 22(21), 2688–2690. https://doi.org/10.1093/bioinformatics/btl446

• Stamatakis, A. (2014). RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics, 30(9), 1312-–1313. https://doi.org/10.1093/bioinformatics/btu033


• For this work we have used the seqs2occupancy.py (https://github.com/tauanajc/phylo_scripts/blob/master/seqs2occupancy.py) script created by Dr. Tauana Cunha and Dr. Bruno Medeiros.

• We have also used the scripts and commands available in the PHYLUCE tutorial (https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-1.html).

