#!/bin/bash

# Activate conda environment for tools
source /raid/VIDRL-USERS/HOME/aduncan/miniforge3/etc/profile.d/conda.sh
conda activate rna-tools

# Set directories
RAW_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/rawdata"
INTER_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/intermediary_files"
PREPROC_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/preproc"
REF_GENOME="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/references/hisat2_index/genome"

# Make sure intermediary and output directories exist
mkdir -p "$INTER_DIR"
mkdir -p "$PREPROC_DIR"

# Loop through R1 files in raw data directory
for R1_FILE in "$RAW_DIR"/*_R1_001.fastq.gz; do
    # Derive sample name
    BASENAME=$(basename "$R1_FILE")
    SAMPLE=$(echo "$BASENAME" | sed 's/_R1_001.fastq.gz//')

    R2_FILE="$RAW_DIR/${SAMPLE}_R2_001.fastq.gz"

    echo "============================"
    echo "Processing sample: $SAMPLE"
    echo "============================"

    # Step 1: Extract UMIs into read names
    umi_tools extract \
        --extract-method=string \
        --bc-pattern=NNNNNNNN \
        --stdin="$R1_FILE" \
        --stdout="$INTER_DIR/${SAMPLE}_R1_extracted.fastq.gz" \
        --read2-in="$R2_FILE" \
        --read2-out="$INTER_DIR/${SAMPLE}_R2_extracted.fastq.gz" \
        --log="$INTER_DIR/${SAMPLE}_extract.log"

    # Step 2: Align extracted reads with hisat2
    hisat2 -x "$REF_GENOME" \
        -1 "$INTER_DIR/${SAMPLE}_R1_extracted.fastq.gz" \
        -2 "$INTER_DIR/${SAMPLE}_R2_extracted.fastq.gz" \
        -S "$INTER_DIR/${SAMPLE}.sam"

    # Step 3: Convert SAM to sorted BAM
    samtools view -bS "$INTER_DIR/${SAMPLE}.sam" | samtools sort -o "$INTER_DIR/${SAMPLE}.sorted.bam"

    # Remove SAM file to save space
    rm "$INTER_DIR/${SAMPLE}.sam"

    # Step 4: Deduplicate BAM using umi_tools (UMIs in read names)
    umi_tools dedup \
        --extract-umi-method=read_id \
        -I "$INTER_DIR/${SAMPLE}.sorted.bam" \
        -S "$PREPROC_DIR/${SAMPLE}.dedup.bam" \
        --log="$PREPROC_DIR/${SAMPLE}_dedup.log"

    # Step 5: Get BAM stats
    samtools flagstat "$PREPROC_DIR/${SAMPLE}.dedup.bam" > "$PREPROC_DIR/${SAMPLE}_dedup_stats.txt"

    echo "Sample $SAMPLE complete."
    echo ""
done

echo "============================"
echo "Pipeline run completed!"
echo "============================"
