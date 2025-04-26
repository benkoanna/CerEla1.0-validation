#!/bin/bash
# Select best allele with MAPQ≥40 filtering

INPUT_DIR="output"
OUTPUT_DIR="merged"
MIN_MAPQ=30  # Quality threshold

mkdir -p "$OUTPUT_DIR"

# Merge and keep original file
samtools merge -f -@ 8 "$OUTPUT_DIR/merged_all.bam" "$INPUT_DIR"/chr*_sorted.bam || exit 1

# Filtering - new file with _filtered suffix
samtools view -@ 8 -h "$OUTPUT_DIR/merged_all.bam" | \
awk -v minq=$MIN_MAPQ '
BEGIN { FS=OFS="\t" }
/^@/ { print; next }
{
    if ($5 < minq) next
    snp_id = $1
    sub(/_(ref_[ACGT]|alt_[ACGT])$/, "", snp_id)
    
    if (!best[snp_id] || $5 > best_score[snp_id]) {
        best[snp_id] = $0
        best_score[snp_id] = $5
    }
}
END {
    for (id in best) print best[id]
}' | \
samtools view -@ 8 -b -o "$OUTPUT_DIR/best_alleles_filtered.bam" - || exit 1

# Final output based on filtered data
samtools sort -@ 8 -o "$OUTPUT_DIR/final_snps_filtered.bam" "$OUTPUT_DIR/best_alleles_filtered.bam"
samtools index "$OUTPUT_DIR/final_snps_filtered.bam"

# Statistics
total_snps=$(samtools view -c "$OUTPUT_DIR/final_snps_filtered.bam")
total_expected=38083
avg_mapq=$(samtools view "$OUTPUT_DIR/final_snps_filtered.bam" | awk '{sum+=$5} END{printf "%.1f", sum/NR}')
ref_count=$(samtools view "$OUTPUT_DIR/final_snps_filtered.bam" | grep -c "_ref_")
alt_count=$((total_snps - ref_count))

echo "Total SNP: $total_snps ($(awk -v t=$total_snps -v e=$total_expected 'BEGIN{printf "%.1f", (t/e)*100}')% of expected)"
echo "Expected SNPs: $total_expected"
echo "Avg. MAPQ: $avg_mapq"
echo "   - $(awk -vq=$avg_mapq 'BEGIN{
    if (q >= 60) print "Excellent quality (MAPQ ≥ 60)";
    else if (q >= 40) print "Good quality (MAPQ ≥ 40)";
    else print "Acceptable quality (MAPQ < 40)"}')"
echo ""
echo "Allele type distribution:"
echo "   Reference (major) alleles: $ref_count ($(awk -v r=$ref_count -v t=$total_snps 'BEGIN{printf "%.1f", (r/t)*100}')%)"
echo "   Alternative (minor) alleles: $alt_count ($(awk -v a=$alt_count -v t=$total_snps 'BEGIN{printf "%.1f", (a/t)*100}')%)"
echo ""
samtools view "$OUTPUT_DIR/final_snps_filtered.bam" | cut -f1 | awk -F'_' '{print $NF}' | sort | uniq -c | awk '{print "   "$1, $2}'

# Index original files - if needed
if [ ! -f "$OUTPUT_DIR/merged_all.bam.bai" ]; then
    samtools index -@ 8 "$OUTPUT_DIR/merged_all.bam"
fi