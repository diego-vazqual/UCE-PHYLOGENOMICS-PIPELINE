# Guide for UCE Phylogenomic Analysis
![logo](https://github.com/diego-vazqual/UCE-PHYLOGENOMICS-PIPELINE/blob/main/ultraconserved-header.png)


## INDEX

1. [Requirements](#1-requirements)
2. [Clean reads](#2-clean-reads)
3. [Contig assembly](#3-contig-assembly)
4. [Target enrichment](#4-target-enrichment)
5. [Alignment](#5-Alignment)
6. [Phylogenomic Analysis](#6-phylogenomic-analysis)
7. [References](7-references)

---

## 1. REQUIREMENTS

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
⚠️ If R1 and R2 have a different number of readings, it will be necessary to match them before continuing.

### 2.2 FASTP
To perform quality filtering (removal of low-quality bases) and cleaning of the FASTQ files (removal of adapter sequences and duplicate reads), we use a script called fastp.sh. This script takes as input a folder containing the raw FASTQ files and generates cleaned sequences in an output directory (e.g., clean-fastq).

An example of running the script would be:
```
bash fastp.sh [input path:~/Desktop/Diego/work_directory/data/fastq] [output path:~/Desktop/Diego/work_directory/data/clean_reads_fastq]
```

### 2.3 CD-HIT-DUP

Before running assemblies with SPAdes, it is recommended to remove potential contamination and duplicate reads present in the raw data. This is done using CD-HIT-DUP, a tool designed to detect duplicate sequences in FASTQ files.

The process is automated through the cd-hit-dup.sh script, which takes the raw data as input and generates cleaned output files ready for assembly.
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

⚠️ Before running `phyluce_assembly_assemblo_spades`, make sure you have activated the environment where Phyluce is installed:
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
⚠️ Make sure that the taxon names in your list exactly match the names of the folders or files in the assembly directory. For this reason, if no list has been prepared in advance, it's useful to extract the names directly from the assembly folder to avoid mistakes or missing taxa. 

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

☁️ To find out how many duplicates we have recovered per taxon, we can recalculate the statistics from step 4.2.1 and compare how many UCE loci have been recovered.

To do this, we repeat the steps from section 4.2.1, but this time using the file `final-taxon_set1-taxa-incomplete.fasta`:
```
phyluce_assembly_explode_get_fastas_file \
    --input final-taxon_set1-taxa-incomplete.fasta \
    --output exploded-fastas-final-set1 \
    --by-taxon
```
```
for i in exploded-fastas-final-set1/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv exploded_fastas.csv;
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

To avoid these issues, taxon names must be renamed using a mapping file `taxa_map.txt` like the example below:
```
Acteon_sp    Taxa1
Acteon_tornatilis    Taxa2
Akera_bullata    Taxa3
Ammonicera_sp    Taxa 4
```
In this file `taxa_map.txt`, the first column contains the original taxon names (as they appear in the alignments), and the second column, separated by a tab, contains the new shortened names. In this example, we use generic identifiers like Taxa1, Taxa2, etc.



Oaprimero hay  que cambiar los nombres de estos 
fueron reconstruidos utilizando RAxML v8.2.12 (Stamatakis, 2006) bajo el modelo GTRGAMMA con 100 réplicas de bootstrap (script single_gene_trees.sh). Los árboles resultantes se usaron para identificar y eliminar ramas terminales anormalmente largas,




List of cited tools, publications, and external scripts used in the workflow. Be sure to include citations for:

* PHYLUCE
* IQ-TREE
* ExaBayes
* Gblocks
* Zorro
* Any external Python or shell scripts used

---

> 🧠 *This README is meant to provide clear step-by-step instructions for setting up and running a full UCE phylogenomic pipeline, based on publicly available tools and custom scripts.*

