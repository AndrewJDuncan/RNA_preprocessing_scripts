#!/bin/bash

# Activate environment
source ~/miniforge3/bin/activate rna-tools

# Set directories
RAW_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/rawdata"
PREPROC_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/preproc"
REF_GENOME="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/references/hg38.fa"
INTER_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/intermediary_files"

mkdir -p "$INTER_DIR" "$PREPROC_DIR"

# Optional: skip SAM generation
SKIP_SAM=true

# Loop over samples
for R1_FILE in "$RAW_DIR"/*_R1_001.fastq.gz; do
    SAMPLE=$(basename "$R1_FILE" | sed 's/_R1_001.fastq.gz//')
    R2_FILE="${RAW_DIR}/${SAMPLE}_R2_001.fastq.gz"

    echo "=========================================="
    echo "Processing sample: $SAMPLE"
    echo "=========================================="

    # Extract UMIs from start of R1 (first 8bp) into RX tag
    umi_tools extract \
        --bc-pattern=NNNNNNNN \
        --stdin="$R1_FILE" \
        --stdout="$INTER_DIR/${SAMPLE}_extracted_R1.fastq.gz" \
        --read2-in="$R2_FILE" \
        --read2-out="$INTER_DIR/${SAMPLE}_extracted_R2.fastq.gz"

    # Align reads to reference genome
    hisat2 -x "$REF_GENOME" \
        -1 "$INTER_DIR/${SAMPLE}_extracted_R1.fastq.gz" \
        -2 "$INTER_DIR/${SAMPLE}_extracted_R2.fastq.gz" \
        -S "$INTER_DIR/${SAMPLE}.sam"

    # Convert SAM to sorted BAM
    samtools view -bS "$INTER_DIR/${SAMPLE}.sam" | samtools sort -o "$INTER_DIR/${SAMPLE}.sorted.bam"

    if [ "$SKIP_SAM" = true ]; then
        rm "$INTER_DIR/${SAMPLE}.sam"
    fi

    # Deduplicate using umi_tools
    umi_tools dedup \
        -I "$INTER_DIR/${SAMPLE}.sorted.bam" \
        -S "$PREPROC_DIR/${SAMPLE}.dedup.bam"

    # BAM flagstat for stats
    samtools flagstat "$PREPROC_DIR/${SAMPLE}.dedup.bam" > "$PREPROC_DIR/${SAMPLE}_dedup_stats.txt"

done

# Run output checker
bash check_pipeline_output.sh

echo "✅ Pipeline run complete!"
