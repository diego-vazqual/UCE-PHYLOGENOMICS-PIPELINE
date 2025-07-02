# UCE-PHYLOPIPE
![logo](https://github.com/diego-vazqual/UCE-PHYLOGENOMICS-PIPELINE/blob/main/ultraconserved-header.png)
Guide for UCE Phylogenomic Analysis

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
Example assembly.conf:
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

### 4.2 ExaBayes

Run with `exabayes.sh`, analyze convergence with `postProcParam`, `sdsf`, `credibleSet`, etc.

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

