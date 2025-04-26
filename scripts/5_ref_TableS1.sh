#!/bin/bash
# Filter FASTA sequences based on TableS1 patterns

patterns_file="patterns.txt"
input_file="DeerSequences_2SNP_.fa"
output_file="DeerSequences_2SNP_out.fa"
log_file="process2.log"
error_file="error2.log"

# Initialize output files
> "$output_file"
> "$log_file"
> "$error_file"

echo "=== Filtering process started: $(date) ===" >> "$log_file"

while IFS= read -r pattern; do
    echo "Searching for pattern: $pattern" >> "$log_file"
    
    if grep -wi "$pattern" "$input_file" >> "$output_file"; then
        echo "Match found: $pattern" >> "$log_file"
    else
        echo "No match: $pattern" >> "$error_file"
        echo "ERROR: '$pattern' not found in input file." >> "$log_file"
    fi
done < "$patterns_file"

echo "=== Filtering process completed: $(date) ===" >> "$log_file"

echo "Results saved to '$output_file'"
echo "Errors logged to '$error_file'"
echo "Full process log at '$log_file'"