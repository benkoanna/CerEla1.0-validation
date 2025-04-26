#!/bin/bash
# Convert single-line FASTA to multi-line format

sed 's/ /\n/g' DeerSequences_1SNP_out.fa > DeerSequences_1SNP_out_.fa
sed 's/ /\n/g' DeerSequences_2SNP_out.fa > DeerSequences_2SNP_out_.fa

# Add allele information to FASTA headers
awk 'BEGIN {FS="\n"; OFS="\n"} 
     /^>/ { 
        if (seq) print header "_ref_" ref "\n" seq; 
        header=$0; 
        seq=""; 
        next 
     } 
     !/^>/ { 
        seq=$0; 
        ref=substr(seq, 151, 1); 
     } 
     END { 
        if (seq) print header "_ref_" ref "\n" seq; 
}' DeerSequences_1SNP_out_.fa > DeerSequences_1SNP_out_final.fa

a következő:
awk 'BEGIN {FS="\n"; OFS="\n"} 
     /^>/ { 
        if (seq) { 
            alt = (length(seq) >= 151) ? substr(seq, 151, 1) : "N"; 
            print header "_alt_" alt "\n" seq; 
        } 
        header=$0; 
        seq=""; 
        next 
     } 
     { 
        seq = seq $0
     } 
     END { 
        if (seq) { 
            alt = (length(seq) >= 151) ? substr(seq, 151, 1) : "N"; 
            print header "_alt_" alt "\n" seq; 
        } 
}' 3_DeerSequences_2SNP_out_.fa > 3_DeerSequences_2SNP_out_final.fa


# Combine reference and alternative allele files
cat DeerSequences_1SNP_out_final.fa DeerSequences_2SNP_out_final.fa > DeerSequences_SNP_final.fa