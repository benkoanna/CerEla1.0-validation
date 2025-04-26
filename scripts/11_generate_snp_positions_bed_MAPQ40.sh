#!/bin/bash
# Generate BED file from alignment results

INPUT_BAM="merged"
OUTPUT_BED="$INPUT_BAM/final_snps_detailed.bed"

echo "Generating detailed BED file..."
samtools view -@ 8 "$INPUT_BAM/final_snps_filtered.bam" | \
awk '
BEGIN { 
    OFS = "\t"
    READ_LENGTH = 301
    SNP_OFFSET = 151
}
{
    read_start = $4
    
    snp_pos = read_start + SNP_OFFSET - 1
    
    read_end = read_start + READ_LENGTH - 1
    
    split($1, parts, "_")
    
    print $3, snp_pos-1, snp_pos, $1, $5, $2, 
          parts[4]"_"substr(parts[5],1,1), read_start, read_end
}' > "$OUTPUT_BED"

echo "Detailed BED file generated: $OUTPUT_BED"
echo "Total SNPs: $(wc -l < "$OUTPUT_BED")"
echo "First 5 entries:"
column -t "$OUTPUT_BED" | head -n 5