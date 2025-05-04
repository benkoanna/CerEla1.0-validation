# 🦌 CerEla1.0 Genome Validation Scripts

This repository contains the set of bioinformatics scripts used for validating the **CerEla1.0 genome assembly**, which represents the *Cervus elaphus* (red deer) reference genome.

The scripts support tasks such as:
- Merging and filtering alignment files (BAM)
- Selecting the best alignments per SNP
- Quality control and basic statistics
- Preparing genome references for BWA alignment

## 📁 Repository Structure

raw_files/
Initial input files used during the processing.

scripts/
Numbered scripts following the validation steps in the correct order.



## 🧬 Main Steps

1. **Reference Preparation**  
   Set up the CerEla1.0 genome for alignment using BWA.

2. **Alignment & Merging**  
   Align reads to individual chromosomes and merge the resulting BAM files.

3. **Filtering**  
   Apply MAPQ-based filtering to retain high-quality alignments.

4. **SNP Validation**  
   Choose the best alignment per SNP and summarize the allele distribution.

5. **Statistics**  
   Report summary metrics such as total SNP count, average MAPQ, and allele frequencies.

## ⚙️ Requirements

- samtools >= 1.10
- bwa
- seqret (from EMBOSS suite)
- awk, grep, cut, sort (standard Unix tools)

Make sure these tools are installed and accessible in your PATH.

## 🚀 Running the Pipeline

Scripts are written in bash.
