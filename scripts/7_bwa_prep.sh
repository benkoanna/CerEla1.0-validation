#!/bin/bash
# Prepare reference genome for BWA alignment

# Create directory structure
mkdir -p bwa_analysis/{genome,scripts,output,reads}

# Process reference genome
cd bwa_analysis/genome || exit
cp /path/to/Celaphus1.0.fasta .

# Standardize chromosome names
sed -i 's/chr_/chr/g' Celaphus1.0.fasta

# Split genome into chromosomes
seqretsplit Celaphus1.0.fasta

# Remove scaffold files
rm -f mkhe*.fasta

# Index each chromosome
for chr in chr*.fasta; do
    bwa index "$chr"
done

echo "Reference genome preparation complete"