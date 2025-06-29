# UCE-PHYLOPIPE

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

* CD-HIT-DUP [Link](https://sites.google.com/view/cd-hit)
* ExaBayes v1.5 [Link](https://cme.h-its.org/exelixis/web/software/exabayes/)
* Fastp [GitHub](https://github.com/OpenGene/fastp)
* IQ-TREE v2.2.5 [Link](https://iqtree.github.io/)
* OpenMPI 4.1.0-10 [Link](https://www.open-mpi.org)
* Phyluce v1.7.2 [Link](https://phyluce.readthedocs.io/en/latest/)
* 
---

## 2. Clean reads

### 2.1 COUNTING THE READ DATA
This command counts the number of reads in files R1 or R2 (in the example R1). It is useful to verify that R1 and R2 have the same number of reads.
The command divides the total number of lines by 4, since each read in FASTQ format occupies 4 lines.

```bash
for i in *_R1_*.fastq.gz; do echo $i; gunzip -c $i | wc -l | awk '{print $1/4}'; done
```
⚠️ If R1 and R2 have a different number of readings, it will be necessary to match them before continuing.

### 2.2 FASTP
Para ejecutar el filtrado de calidad (eliminación de bases de baja calidad) y limpieza de los archivos FASTQ (eliminación adaptadores y lecturas duplicadas), utilizamos un script llamado fastp.sh. Este script toma como entrada una carpeta que contiene los archivos FASTQ brutos y genera las secuencias limpias en una carpeta de salida (por ejemplo, clean-fastq).

An example of running the script would be:
```bash
bash fastp.sh ~/Desktop/Diego/fastq ~/Desktop/Diego/clean-fastq
```

### 2.3 ExaBayes Installation

Install via source with GCC 10 and MPI. See ExaBayes manual.

### 2.4 CD-HIT Installation

```bash
conda install -c bioconda cd-hit
conda install -c bioconda cd-hit-auxtools
```

### 2.5 ASTRAL Installation

Download JAR file and add path to environment or use absolute path.

### 2.6 Fastp Installation

```bash
conda install -c bioconda fastp
```

---

## 3. PHYLOGENOMIC ANALYSIS

* **Counting Reads**
* **Pre-processing** (fastp, cd-hit-dup, adapter trimming)
* **Assembly** using SPAdes via `phyluce_assembly_assemblo_spades`
* **UCE Loci Identification**
* **Loci Extraction and Filtering**
* **Alignment with MAFFT**
* **Masking with Gblocks / Zorro**
* **Building Final Matrices** (50% occupancy, filtered loci)

---

## 4. DOWNSTREAM ANALYSIS

### 4.1 IQ-TREE

Maximum Likelihood trees with GHOST models:

```bash
iqtree2 -s input.phylip -m GTR+FO*H4 -B 1500 -T AUTO
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

