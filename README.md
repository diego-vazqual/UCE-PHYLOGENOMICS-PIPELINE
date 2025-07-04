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
   - [5.4 Remove long terminal branches](#54-remove-long-terminal-branches)  
   - [5.5 Preparing data for phylogenomic analysis](#55-preparing-data-for-phylogenomic-analysis)  
   - [5.6 Final data matrices](#56-final-data-matrices)  
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
This command counts the number of reads in files R1 or R2 (in the example R1). It is useful to verify that R1 and R2 have the same number of reads.
The command divides the total number of lines by 4, since each read in FASTQ format occupies 4 lines.

```
for i in *_R1_*.fastq.gz; do echo $i; gunzip -c $i | wc -l | awk '{print $1/4}'; done
```
⚠️ Important: If R1 and R2 have a different number of readings, it will be necessary to match them before continuing.

### 2.2 FASTP
To perform quality filtering (removal of low-quality bases) and cleaning of the FASTQ files (removal of adapter sequences and duplicate reads), we use a script called `fastp.sh`. This script takes as input a folder containing the raw FASTQ files and generates cleaned sequences in an output directory (e.g., clean-fastq).

An example of running the script would be:
```
bash fastp.sh [input path:~/Desktop/Diego/work_directory/data/fastq] [output path:~/Desktop/Diego/work_directory/data/clean_reads_fastq]
```

### 2.3 CD-HIT-DUP

Before running assemblies with SPAdes, it is recommended to remove potential contamination and duplicate reads present in the raw data. This is done using CD-HIT-DUP, a tool designed to detect duplicate sequences in FASTQ files.

The process is automated through the `cd-hit-dup.sh` script, which takes the raw data as input and generates cleaned output files ready for assembly.
```
bash cd-hit-dup.sh [input path:~/Desktop/Diego/work_directory/data/clean_reads_fastq] [output path:~/Desktop/Diego/work_directory/data/clean_reads_cdhitdup]
```

## 3. Contig assembly

### 3.1 SPADES
The SPAdes program integrated into Phyluce requires a configuration file (assembly.conf) with a [samples] header and a list of each sample name followed by the path to its cleaned sequences.

From within the folder containing your cleaned sample directories:
```
echo "[samples]" > ../assembly.conf
for i in *; do echo "$i:/home/intern/Desktop/Diego/work_directory/data/clean_reads_cdhitdup/$i/"; done >> ../assembly.conf
```
✏️ Example assembly.conf:
```
[samples]
Acteon_sp:/home/intern/Desktop/Diego/work_directory/data/clean_reads_cdhitdup/Acteon_sp
Acteon_tornatilis:/home/intern/Desktop/Diego/data/work_directory/clean_reads_cdhitdup/Acteon_tornatilis
Akera_bullata:/home/intern/Desktop/Diego/work_directory/data/clean_reads_cdhitdup/Akera_bullata
Ammonicera_sp:/home/intern/Desktop/Diego/work_directory/data/clean_reads_cdhitdup/Ammonicera_sp
```
Run SPAdes assembly (from within the folder containing config file):

⚠️ Important: Before running `phyluce_assembly_assemblo_spades`, make sure you have activated the environment where Phyluce is installed:
```
conda activate phyluce-1.7.2
```
```
phyluce_assembly_assemblo_spades \
    --conf assembly.conf \
    --output data/spades_assemblies \
    --cores 30
```
📌 CPU threads depending on dataset size. This example used 117 samples and 2225 captured UCEs.

## 4. Target enrichment

### 4.1 Finding UCE loci
After assembling contigs from the clean reads, the next step is to identify which contigs correspond to UCE loci (Ultra Conserved Elements) and exclude non-target sequences.

From this point forward, maintaining a clean and consistent folder structure becomes especially important to avoid errors and streamline the workflow.

To identify UCE loci, the following command is used, which compares the contigs with a set of probes. In this example, the command is executed from within the `Diego/` directory, where the probe file `Probeset-70nt.fasta` is located:
```
phyluce_assembly_match_contigs_to_probes \
    --contigs data/spades_assemblies \
    --probes Probeset-70nt.fasta \
    --output uce-search-results \
    --keep-duplicates duplicates.txt
    --csv uce_serach_results.csv
```

By default, the Phyluce function `phyluce_assembly_match_contigs_to_probes` filters out UCE loci and contigs identified as duplicates in the dataset. These are identified as duplicates when probes designed for different UCE loci retrieve the same contig, or when multiple contigs, supposedly representing different genomic regions, are matched with probes targeting a single UCE locus. 

To recover these duplicates, the `--keep-duplicates` option is used. This option allows the recovery, in the output file `duplicates.txt`, of duplicates detected per taxon at each locus. In subsequent processing steps, we will recover part of these loci.

☁️ By default, the search results are displayed in the terminal. These results include the number of UCE loci captured within the total contigs for each taxon, as well as the UCE loci and contigs identified as duplicates that were removed. However, by using the --csv option, these results can be saved to a CSV file.

### 4.2 Extracting UCE loci

Once the UCE loci have been identified, the next step is to define the taxa to be included in the analysis. To do this, it is necessary to compile a list with the names of these taxa and create a configuration file indicating which UCE loci are present in each one.

If you don't have a list of the taxa you want to use, but you do have a folder containing the corresponding files for those taxa (e.g., spades_assemblies), you can retrieve their names using the following command:
```
ls -1 [Name of the directory containing the files of the taxa you want to use]
```
⚠️ Important: Make sure that the taxon names in your list exactly match the names of the folders or files in the assembly directory. For this reason, if no list has been prepared in advance, it's useful to extract the names directly from the assembly folder to avoid mistakes or missing taxa. 

Based on the list of taxa we want to use, we will create a file using the `nano` editor, which we will name ` taxon-set.conf`.  On the first line of the file, you should write the name you want to assign to the taxon set, for example, `[taxon_set1]`. Then, paste the names of the taxa, one per line.

✏️ Example taxon-set.conf:
```
[taxon_set1]
Acteon_sp
Acteon_tornatilis
Akera_bullata
Ammonicera_sp
```
💡 If you later want to change the taxon set, you can edit this configuration file (` taxon-set.conf`) directly using the `nano` editor and add a new list with a different header, such as `[taxon_set2]`.

Having configured the taxon list, we can run the following command to generate the initial list of UCE loci enriched in each taxon:
```
phyluce_assembly_get_match_counts \
    --locus-db uce-search-results/probe.matches.sqlite \
    --taxon-list-config taxon-set.conf \
    --taxon-group 'taxon_set1' \
    --incomplete-matrix \
    --output taxon_set1-taxa-incomplete.conf
```
The above command generates an output called `all-taxa-incomplete.conf`. We need this output in FASTA format for the next step. To do this, we run the following command:
```
phyluce_assembly_get_fastas_from_match_counts \
    --contigs data/spades_assemblies \
    --locus-db uce-search-results/probe.matches.sqlite \
    --match-count-output taxon_set1-taxa-incomplete.conf \
    --output taxon_set1-taxa-incomplete.fasta \
    --incomplete-matrix taxon_set1-taxa-incomplete.incomplete \
    --log-path log
```
#### 4.2.1 Analyse the statistics for each taxon
☁️ If we want to analyze the statistics for each taxon, we need to separate the all-taxa-incomplete.fasta file into individual FASTA files for each taxon, each containing its respective captured UCEs.

To do this, you must first split the `taxon_set1-taxa-incomplete.fasta` file by taxon using the following command:
```
phyluce_assembly_explode_get_fastas_file \
    --input taxon_set1-taxa-incomplete.fasta \
    --output exploded-fastas \
    --by-taxon
```
Then, extract the statistics using the following command:
```
echo "Sample ID,UCE loci,total bp,mean length,95 CI length,min length,max length, median legnth, contigs >1kb" > exploded_fastas.csv
for i in exploded-fastas/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv exploded_fastas.csv;
done
```
The results will be printed to the terminal and saved to a csv file, `exploded_fastas.csv` in this example.

Example of statistics output:
```
## samples id,UCE loci,total bp,mean length,95 CI length,min length,max length,median legnth,contigs >1kb
Acteon_sp.unaligned.fasta,1395,756820,477.8030888030888,3.765000046699546,308,577,520.0,12
Acteon_tornatilis.unaligned.fasta,1455,776310,605.4563888030888,4.657000046699546,208,420,402.0,12
```

#### 4.2.2 Recovering duplicate UCE loci
Recall that in step 4.1, all UCE loci identified as duplicates for a given taxon are discarded—this information is stored in the `duplicates.txt` file. This means that any UCE locus for which more than one contig was found is completely excluded from the dataset for that taxon, resulting in a loss of information.

To recover part of this information, we process the `duplicates.txt` file to select, for each duplicated UCE locus, the best representative contig.

First, we create a directory where we will carry out the recovery of the duplicates:
```
mkdir -p duplicates
```
Then, we move the list, which indicates which duplicates have been identified for each UCE locus and taxon into, this directory:
```
mv duplicates.txt duplicates/
```
Next, we run the Python script `phyluce_assembly_parse_duplicates_file.py` to recover these duplicates in FASTA format:
```
python ./phyluce_assembly_parse_duplicates_file.py --contigs data/spades_assemblies --duplicates-file duplicates/duplicates.txt --output duplicates/duplicates.fasta
```
And inside ` duplicates/`:
```
cd duplicates/
```
```
phyluce_assembly_explode_get_fastas_file \
    --input duplicates.fasta \
    --output exploded-duplicates-fastas \
    --by-taxon
```

To select the most representative contig, we rely on two criteria: the longest contig and the highest percentage identity relative to the corresponding locus in other taxa. By default, the script `phyluce_assembly_parse_duplicates_file.py` labels this contig as DUPE1. To recover these DUPE1 contigs, we run the following commands:
```
cd exploded-duplicates-fastas
```
```
cat *-DUPE1.unaligned.fasta >> duplicates_DUPE1.fasta
```
```
sed -i 's/_DUPE1//g' duplicates_DUPE1.fasta
```

Now, the selected contigs stored in the FASTA file `duplicates_DUPE1.fasta` are concatenated with the rest of the contigs from the FASTA file `taxon_set1-taxa-incomplete.fasta`.

First, we make a copy of the `duplicates_DUPE1.fasta` file into the main working directory, in this case `work_directory/`:
```
cp duplicates_DUPE1.fasta ../../ 
```
And finally, we concatenate both files:

```
cat taxon_set1-taxa-incomplete.fasta duplicates_DUPE1.fasta >> final-taxon_set1-taxa-incomplete.fasta
```

☁️ To find out how many duplicates we have recovered per taxon, we can recalculate the statistics from step `4.2.1 ` and compare how many UCE loci have been recovered.

To do this, we repeat the steps from section 4.2.1, but this time using the file `final-taxon_set1-taxa-incomplete.fasta`:
```
phyluce_assembly_explode_get_fastas_file \
    --input final-taxon_set1-taxa-incomplete.fasta \
    --output exploded-fastas-final-set1 \
    --by-taxon
```
```
echo "Sample ID,UCE loci,total bp,mean length,95 CI length,min length,max length, median legnth, contigs >1kb" > exploded_fastas_final.csv
for i in exploded-fastas-final-set1/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv exploded_fastas_final.csv;
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

### 5.4 Remove long terminal branches
Even after the initial cleaning steps, alignments may still contain atypical sequences, which could represent chimeras, contaminants, or non-homologous alignments. These sequences can introduce artifacts in subsequent phylogenetic inference, particularly by generating abnormally long terminal branches, which may lead to the phenomenon known as long branch attraction (LBA). This occurs when highly divergent sequences are incorrectly grouped due to elevated substitution rates or high sequence dissimilarity.

To mitigate this issue, we identified and removed abnormally long terminal branches using the following steps.

To identify and remove abnormally long terminal branches, we first build gene trees from each alignment. For this step, alignments must be in  `.phylip` format.

⚠️ Important: The PHYLIP format requires taxon names to be no longer than 10 characters. If names exceed this limit, programs may throw read errors, truncate names, or produce incorrect results.

✏️To avoid these issues, taxon names must be renamed using a mapping file `taxa_map.txt` like the example below:
```
Acteon_sp    Taxa1
Acteon_tornatilis    Taxa2
Akera_bullata    Taxa3
Ammonicera_sp    Taxa 4
```
In this file `taxa_map.txt`, the first column contains the original taxon names (as they appear in the alignments), and the second column, separated by a tab, contains the new shortened names. In this example, we use generic identifiers like Taxa1, Taxa2, etc.

The script `alignment_phylip_format.py` renames taxon names and converts alignments from nexus format to phylip format. 
This script takes as input the folder containing the alignments in nexus format (`mafft-phylip-nexus-internal-no-trimmed-gblocks-clean`) and the mapping file with the original and new taxon names (`taxa_map.txt`). The output is saved in a directory called `taxa_map.txt`.

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

The script takes as input the individual gene trees from the `genes_tree/RAxML_bipartitions` directory and the alignments in nexus format (`mafft-nexus-internal-no-trimmed-gblocks-clean`). As output, it generates the filtered individual gene trees in the `clean_genes_trees` directory and **directly modifies the original alignments**.

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
```
```
echo "Sample ID,UCE loci,total bp,mean length,95 CI length,min length,max length, median legnth, contigs >1kb" > fasta_lengths_clean.csv
for i in exploded-fastas-clean/*.fasta; do  
    phyluce_assembly_get_fasta_lengths --input "$i" --csv >> fasta_lengths_clean.csv  
done
```

### 5.5 Final data matrices
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
```
 ```
echo "Sample ID,UCE loci,total bp,mean length,95 CI length,min length,max length, median legnth, contigs >1kb" > fasta_lengths_50p.csv
for i in exploded-fastas-50p/*.fasta; do  
    phyluce_assembly_get_fasta_lengths --input "$i" --csv >> fasta_lengths_50p.csv  
done
```
### 5.6 Preparing data for phylogenomic analysis
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

• Crotty, S.M., Minh, B.Q., Bean, N.G., Holland, B.R., Tuke, J., Jermiin, L.S., Haeseler, A.V., 2019. GHOST: Recovering historical signal from heterotachously-evolved sequence alignments. Syst. Biol. 69, 249–264. https://doi.org/10.1093/sysbio/syz051

• Wang, H.C., Minh, B.Q., Susko, S., Roger, A.J., 2018. Modeling site heterogeneity with posterior mean site frequency profiles accelerates accurate phylogenomic estimation. Syst. Biol., 67:216-235. https://doi.org/10.1093/sysbio/syx068

• Kalyaanamoorthy, S., Minh, B.Q., Wong, T.K.F., Haeseler, A.V., Jermiin, L.S., 2017. ModelFinder: Fast Model Selection for Accurate Phyloge- netic Estimates, Nature Methods, 14:587–589. https://doi.org/10.1038/nmeth.4285

• Faircloth, B.C., 2016. PHYLUCE is a software package for the analysis of conserved genomic loci. Bioinformatics 32, 786–788. https://doi.org/10.1093/bioinformatics/btv646.

• Andre J. Aberer, Kassian Kobert, Alexandros Stamatakis, ExaBayes: Massively Parallel Bayesian Tree Inference for the Whole-Genome Era, Molecular Biology and Evolution, Volume 31, Issue 10, October 2014, Pages 2553–2556, https://doi.org/10.1093/molbev/msu236

• Lam-Tung Nguyen, Heiko A. Schmidt, Arndt von Haeseler, Bui Quang Minh, IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies, Molecular Biology and Evolution, Volume 32, Issue 1, January 2015, Pages 268–274, https://doi.org/10.1093/molbev/msu300

• Fu, L., Niu, B., Zhu, Z., Wu, S., & Li, W. (2012). CD-HIT: accelerated for clustering the next-generation sequencing data. Bioinformatics (Oxford, England), 28(23), 3150–3152. https://doi.org/10.1093/bioinformatics/bts565

• For this work we have used the fastqCombinePairedEnd.py (https://github.com/enormandeau/Scripts) script created by Dr. Enric Normandeau, and the seqs2occupancy.py (https://github.com/tauanajc/phylo_scripts/blob/master/seqs2occupancy.py) script created by Dr. Tauana Cunha and Dr. Bruno Medeiros.

• We have also used the scripts and commands available in the PHYLUCE tutorial (https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-1.html), and the script zorro.py created by Dr. Tauana Cunha that has been modified for the present study.

