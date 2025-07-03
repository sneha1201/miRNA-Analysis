# miRNA Processing Pipeline

This repository contains two modular Python scripts that together form a complete miRNA-seq data processing pipeline. The pipeline performs preprocessing, alignment, miRNA prediction, annotation, and count matrix generation for both known and novel miRNAs.

---

## 🧩 Pipeline Overview

| Step | Script | Description |
|------|--------|-------------|
| 1–7  | `mirna_alignment.py` | Preprocessing raw FASTQ files: quality control, adapter trimming, alignment to reference genome, BAM conversion, and read collapsing |
| 8–12 | `mirna_after_bam_count.py` | Mapping collapsed reads, predicting known & novel miRNAs using miRDeep2, generating GFFs, and counting miRNA expressions |

---

## 🔧 Requirements

Ensure the following tools are installed and accessible in your environment:

- [fastqc]
- [fastp]
- [seqkit]
- [bowtie]
- [samtools]
- [fastx_toolkit]
- [miRDeep2]
- [featureCounts (from Subread)]

Python version: ≥3.6

---

## 🚀 Usage

### 1. Preprocessing and Alignment

```bash
python mirna_alignment.py \
    --raw /path/to/raw_fastq \
    --output /path/to/output_dir \
    --adapter /path/to/adapter.fa \
    --reference /path/to/bowtie_index_prefix \
    --threads 8


This script performs:

Raw FASTQ quality check (FastQC, SeqKit)

Adapter trimming (fastp)

Post-trimming quality check

Bowtie alignment (SAM/BAM/Index)

Alignment stats report

Merges and collapses reads to prepare for miRDeep2

2. miRNA Prediction and Counting

python mirna_after_bam_count.py \
    --collapsed /path/to/output_dir/5_Collapsed \
    --output /path/to/output_dir \
    --reference /path/to/bowtie_index \
    --ref_fasta_cleaned /path/to/genome_clean.fa \
    --mature /path/to/mature.fa \
    --hairpin /path/to/hairpin.fa \
    --bam_folder /path/to/output_dir/4_Alignment \
    --threads 8


This script performs:

Mapping collapsed reads using mapper.pl

Known miRNA prediction (miRDeep2)

GFF and counts table generation using featureCounts

Novel miRNA discovery and count generation


Output Directory Structure

output_dir/
├── 1_QC_raw/                 # FastQC + seqkit raw
├── 2_Clean_data_miRNA/      # Trimmed reads (.fq.gz)
├── 3_QC_clean/              # FastQC + seqkit clean
├── 4_Alignment/             # SAM, BAM, BAM index, stats
├── 5_Collapsed/             # Collapsed FASTA
├── 6_Mapper/                # ARF file for miRDeep2
├── 7_miRDeep2/              # Known miRNA prediction output
├── 8_Known_miRNA/           # BED, GFF, counts (known)
├── 9_Novel_miRNA/           # BED, GFF, counts (novel)
├── pipeline_step1_7.log     # Log from preprocessing script
└── pipeline_step8_12.log    # Log from miRNA prediction script


🧬 Input File Notes
Adapter file: A FASTA file of known adapters (for use with fastp)

Reference index: Prefix of a Bowtie genome index (use bowtie-build)

Collapsed reads: Output from fastx_collapser in the first script

Genome FASTA (cleaned): Required by miRDeep2.pl, no line breaks in sequences

mature.fa / hairpin.fa: miRBase or other curated mature/hairpin sequences


To run the pipeline on a new dataset, first execute mirna_alignment.py, then pass the output to mirna_after_bam_count.py.


