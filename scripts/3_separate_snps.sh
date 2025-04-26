#!/bin/bash
# Separate SNPs into reference and alternative alleles

input_fasta="DeerSequences.fa"
ref_output="DeerSequences_1SNP.fa"
alt_output="DeerSequences_2SNP.fa"

# Process SNPs
sed "s/\[//g" "$input_fasta" | sed "s/\]//g" | sed "s/\///g" | cut -c 1-151,153-302 > "$ref_output"
sed "s/\[//g" "$input_fasta" | sed "s/\]//g" | sed "s/\///g" | cut -c 1-150,152-302 > "$alt_output"

# Convert line endings
sed -i 's/\r$//' "$ref_output"
sed -i 's/\r$//' "$alt_output"

# Convert to single-line format
awk '/^>/ {if (seq) print seq; seq=$0" "; next} {seq=seq""$0} END {print seq}' "$ref_output" > "DeerSequences_1SNP_.fa"
awk '/^>/ {if (seq) print seq; seq=$0" "; next} {seq=seq""$0} END {print seq}' "$alt_output" > "DeerSequences_2SNP_.fa"

echo "SNP separation complete"