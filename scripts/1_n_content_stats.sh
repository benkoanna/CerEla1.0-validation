#!/bin/bash
# Calculate N-content statistics in a FASTA file

input_fasta="DeerSequences.fa"

awk '/^>/ { 
       if (seqlen) {
          ratio = ncount / seqlen * 100;
          if (ratio < 10) under_10++;
          else if (ratio <= 25) under_25++;
          else if (ratio <= 50) under_50++;
          else over_50++;
          total++;
       }
       seqlen=0; ncount=0; next
     } 
     { 
       seqlen += length($0); 
       ncount += gsub(/[Nn]/, "") 
     } 
     END { 
       if (seqlen) {
          ratio = ncount / seqlen * 100;
          if (ratio < 10) under_10++;
          else if (ratio <= 25) under_25++;
          else if (ratio <= 50) under_50++;
          else over_50++;
          total++;
       }
       print "N-content statistics:";
       print "----------------------------------";
       print "Sequences with <10% N-content: " under_10 " (" under_10/total*100 "%)";
       print "Sequences with 10-25% N-content: " under_25 " (" under_25/total*100 "%)";
       print "Sequences with 25-50% N-content: " under_50 " (" under_50/total*100 "%)";
       print "Sequences with >50% N-content: " over_50 " (" over_50/total*100 "%)";
       print "----------------------------------";
       print "Total sequences: " total;
     }' "$input_fasta"