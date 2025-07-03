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
bash fastp.sh [input path:~/Desktop/Diego/data/fastq] [output path:~/Desktop/Diego/data/clean_reads_fastq]
```

### 2.3 CD-HIT-DUP

Before running assemblies with SPAdes, it is recommended to remove potential contamination and duplicate reads present in the raw data. This is done using CD-HIT-DUP, a tool designed to detect duplicate sequences in FASTQ files.

The process is automated through the cd-hit-dup.sh script, which takes the raw data as input and generates cleaned output files ready for assembly.
```
bash cd-hit-dup.sh [input path:~/Desktop/Diego/data/clean_reads_fastq] [output path:~/Desktop/Diego/data/clean_reads_cdhitdup]
```

## 3. Contig assembly

### 3.1 SPADES
The SPAdes program integrated into Phyluce requires a configuration file (assembly.conf) with a [samples] header and a list of each sample name followed by the path to its cleaned sequences.

From within the folder containing your cleaned sample directories:
```
echo "[samples]" > ../assembly.conf
for i in *; do echo "$i:/home/intern/Desktop/Diego/data/clean_reads_cdhitdup/$i/"; done >> ../assembly.conf
```
✏️ Example assembly.conf:
```
[samples]
Acteon_sp:/home/intern/Desktop/Diego/data/clean_reads_cdhitdup/Acteon_sp
Acteon_tornatilis:/home/intern/Desktop/Diego/data/clean_reads_cdhitdup/Acteon_tornatilis
Akera_bullata:/home/intern/Desktop/Diego/data/clean_reads_cdhitdup/Akera_bullata
Ammonicera_sp:/home/intern/Desktop/Diego/data/clean_reads_cdhitdup/Ammonicera_sp
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
```
By default, the Phyluce function `phyluce_assembly_match_contigs_to_probes` filters out UCE loci and contigs identified as duplicates in the dataset. These are identified as duplicates when probes designed for different UCE loci retrieve the same contig, or when multiple contigs, supposedly representing different genomic regions, are matched with probes targeting a single UCE locus. 

To recover these duplicates, the `--keep-duplicates` option is used. This option allows the recovery, in the output file `duplicates.txt`, of duplicates detected per taxon at each locus. In subsequent processing steps, we select among these duplicates, for each UCE locus, the contig with the greatest length and highest percentage of identity relative to the corresponding locus in other taxa.

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
### 4.3 ASTRAL

Generate gene trees with IQ-TREE, then run:

```bash
java -jar astral.jar -i *.treefile -o output.treefile
```

---

## 5. REFERENCES

List of cited tools, publications, and external scripts used in the workflow. Be sure to include citations for:

* PHYLUCE
* IQ-TREE
* ExaBayes
* Gblocks
* Zorro
* Any external Python or shell scripts used

---

> 🧠 *This README is meant to provide clear step-by-step instructions for setting up and running a full UCE phylogenomic pipeline, based on publicly available tools and custom scripts.*

