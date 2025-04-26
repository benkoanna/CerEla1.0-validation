#!/bin/bash
# Merge and process SNP data with TableS1 information

# Process BED file
awk 'BEGIN {OFS="\t"} {
    pos_mb = sprintf("%.2f", $3/1000000)
    gsub(",", ".", pos_mb)
    
    strand = ($6 == 0 ? "+" : ($6 == 16 ? "-" : $6))

    split($4, parts, "_")
    snp_id = parts[1] "_" parts[2] "_" parts[3] "_" parts[4]
    allele_type = (parts[5] == "ref" ? "ref" : "alt")
    allele_base = parts[6]
    
    print $1, pos_mb, snp_id "\t" allele_type "\t" allele_base, strand
}' final_snps_detailed.bed > final_snps_detailed_mod.txt

# Process TableS1
awk -F'\t' 'BEGIN {OFS="\t"} NR>1 {
    est_pos = sprintf("%.2f", $11/1000000)
    gsub(",", ".", est_pos)
    
    chr = $4
    sub(/^34$/, "X", chr)
    
    print $5, $1, "chr" chr, $7, $8, est_pos
}' TableS1_CervusElaphus_Final_Linkage_Map.txt > Table.txt

# Merge files
awk -F'\t' 'BEGIN {OFS = FS}
NR==FNR && FNR>1 {
    table_snp[$2] = $0
    chr_table[$2] = $3
    next
}
{
    snp_id = $3
    chr_snp = $1
    if (snp_id in table_snp && chr_snp == chr_table[snp_id]) {
        print $0, table_snp[snp_id]
    }
}' Table.txt final_snps_detailed_mod.txt > merged_SNPs.txt

# Add headers
echo -e "Chr\tPos(Mb)\tSNP_ID\tAllele_Type\tBase\tStrand\tCEL.order\tSNP.Name\tChromosome\tcMPosition.Female\tcMPosition.Male\tEstimated.Mb.Position" > header.txt
cat header.txt merged_SNPs.txt > final_merged_SNPs.txt
rm header.txt merged_SNPs.txt
echo "Merged file created: final_merged_SNPs.txt"