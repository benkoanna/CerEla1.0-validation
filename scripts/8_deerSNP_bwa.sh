#!/bin/bash
# Perform BWA alignment of SNP sequences

# Set directories
READS_DIR="reads"
REF_GENOME_DIR="genome"
OUTPUT_DIR="output"
LOG_DIR="log"
INPUT_FASTA="DeerSequences_SNP_final.fa"

# Check input file
if [[ ! -f "$READS_DIR/$INPUT_FASTA" ]]; then
    echo "Error: Input FASTA file not found" | tee -a "$LOG_DIR/errors.log"
    exit 1
fi

# Create output directories
mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# Check and create indexes if needed
for CHR in "$REF_GENOME_DIR"/chr*.fasta; do
    if [[ ! -f "${CHR}.bwt" ]]; then
        echo "Index file missing: $CHR. Indexing..." | tee -a "$LOG_DIR/indexing.log"
        bwa index "$CHR" 2>> "$LOG_DIR/indexing.log"
    fi
done

# Main processing loop
for CHR in "$REF_GENOME_DIR"/chr*.fasta; do
    CHR_NAME=$(basename "$CHR" .fasta)
    CURRENT_TIME=$(date "+%Y-%m-%d %H:%M:%S")
    
    echo "[$CURRENT_TIME] Processing: $CHR_NAME" | tee -a "$LOG_DIR/processing.log"
    
    # Run BWA-MEM with specified parameters for FASTA input
    bwa mem -k 15 -w 100 -T 20 -t 8 "$CHR" "$READS_DIR/$INPUT_FASTA" 2>> "$LOG_DIR/bwa_${CHR_NAME}.log" | \
    samtools view -@ 4 -b - 2>> "$LOG_DIR/samtools_${CHR_NAME}.log" | \
    samtools sort -@ 4 -o "$OUTPUT_DIR/${CHR_NAME}_sorted.bam" - 2>> "$LOG_DIR/samtools_${CHR_NAME}.log"
    
    # Create index
    samtools index "$OUTPUT_DIR/${CHR_NAME}_sorted.bam" 2>> "$LOG_DIR/indexing.log"
    
    # Basic statistics
    TOTAL_READS=$(samtools view -c "$OUTPUT_DIR/${CHR_NAME}_sorted.bam")
    MAPPED_READS=$(samtools view -c -F 4 "$OUTPUT_DIR/${CHR_NAME}_sorted.bam")
    
    echo "[$CURRENT_TIME] $CHR_NAME done. Total reads: $TOTAL_READS, Mapped: $MAPPED_READS" | tee -a "$LOG_DIR/processing.log"
    echo "Mapping rate: $((100 * MAPPED_READS / TOTAL_READS))%" | tee -a "$LOG_DIR/processing.log"
done

# Final summary
echo "Processing completed!" | tee -a "$LOG_DIR/success.log"
echo "BAM files: $OUTPUT_DIR" | tee -a "$LOG_DIR/success.log"
echo "Log files: $LOG_DIR" | tee -a "$LOG_DIR/success.log"