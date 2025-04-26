#!/bin/bash
# Select best allele for each SNP

INPUT_DIR="output"
OUTPUT_DIR="merged"


mkdir -p "$OUTPUT_DIR"

# Merge all chromosome BAMs
echo "1. Merging BAM files..."
samtools merge -f -@ 8 "$OUTPUT_DIR/merged_all.bam" "$INPUT_DIR"/chr*_sorted.bam || exit 1

# Select best alignment per SNP
echo "2. Selecting best allele for each SNP..."
samtools view -@ 8 -h "$OUTPUT_DIR/merged_all.bam" | \
awk '
BEGIN { FS=OFS="\t" }
/^@/ { print; next }
{
    snp_id = $1
    sub(/_(ref_[ACGT]|alt_[ACGT])$/, "", snp_id)

    if (!best[snp_id] || $5 > best_score[snp_id]) {
        best[snp_id] = $0
        best_score[snp_id] = $5
        best_type[snp_id] = ($1 ~ /_ref_/) ? "ref" : "alt"
    }
}
END {
    for (id in best) print best[id]
}' | \
samtools view -@ 8 -b -o "$OUTPUT_DIR/best_alleles.bam" - || exit 1

# Sort and index final BAM
echo "3. Creating final file..."
samtools sort -@ 8 -o "$OUTPUT_DIR/final_snps.bam" "$OUTPUT_DIR/best_alleles.bam"
samtools index "$OUTPUT_DIR/final_snps.bam"

# Generate statistics
echo "4. Final statistics:"
total_snps=$(samtools view -c "$OUTPUT_DIR/final_snps.bam")
avg_mapq=$(samtools view "$OUTPUT_DIR/final_snps.bam" | awk '{sum+=$5} END{printf "%.1f", sum/NR}')
ref_count=$(samtools view "$OUTPUT_DIR/final_snps.bam" | grep -c "_ref_")
alt_count=$((total_snps - ref_count))

echo "Total SNP: $total_snps (100.0%)"
echo "Expected SNPs: 38083"
echo "Avg. MAPQ: $avg_mapq"
echo "   - $(awk -vq=$avg_mapq 'BEGIN{
    if (q >= 60) print "Excellent quality (MAPQ ≥ 60)";
    else if (q >= 40) print "Good quality (40 ≤ MAPQ < 60)";
    else print "Acceptable quality (MAPQ < 40)"}')"
echo "Allele type distribution:"
echo "   Reference (major) alleles: $ref_count ($(awk -v r=$ref_count -v t=$total_snps 'BEGIN{printf "%.1f", (r/t)*100}')%)"
echo "   Alternative (minor) alleles: $alt_count ($(awk -v a=$alt_count -v t=$total_snps 'BEGIN{printf "%.1f", (a/t)*100}')%)"
samtools view "$OUTPUT_DIR/final_snps.bam" | cut -f1 | awk -F'_' '{print $NF}' | sort | uniq -c | awk '{print "   "$1, $2}'