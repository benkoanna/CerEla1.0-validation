#!/bin/bash
# Create pattern file from TableS1
grep -v "^SNP.Name" TableS1_CervusElaphus_Final_Linkage_Map.txt | awk '{print $1}' > patterns.txt