#!/bin/bash
# Create chromosome-specific CSV files

mkdir -p chr_CSV

# Process merged file
awk -F'\t' 'BEGIN {OFS=","} NR>1 {
    chr = $1
    sub(/^chr/, "", chr)
    
    # Create output filename
    output_file = "chr_CSV/chr" chr ".csv"
    
    # Print header if file doesn't exist
    if (!header_printed[chr]) {
        print "Genome,Est_Mbp,f_cM,m_cM" > output_file
        header_printed[chr] = 1
    }
    
    print $2, $12, $10, $11 >> output_file
    
    data[chr][$7] = $2 OFS $12 OFS $10 OFS $11
} END {
    # Process stored data
    for (chr in data) {
        output_file = "chr_CSV/chr" chr ".csv"
        print "Genome,Est_Mbp,f_cM,m_cM" > output_file
        n = asorti(data[chr], indices)
        for (i = 1; i <= n; i++) {
            print data[chr][indices[i]] >> output_file
        }
    }
}' final_merged_SNPs.txt

echo "Chromosome CSV files created in chr_CSV directory"