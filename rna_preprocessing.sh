#!/bin/bash

# ==============================
# RNA-Seq Preprocessing Pipeline
# Clean, updated, nohup-compatible
# ==============================

# Set directories
RAW_DIR="../rawdata"
PREPROC_DIR="../preproc"
INTERMED_DIR="../intermediary_files"
REF_GENOME="../references/hg38.fa"

# Create directories if they don't exist
mkdir -p "$PREPROC_DIR" "$INTERMED_DIR"

# Process each sample
for R1_FILE in $RAW_DIR/*_R1_001.fastq.gz; do
    BASENAME=$(basename "$R1_FILE" | sed 's/_R1_001.fastq.gz//')
    R2_FILE="$RAW_DIR/${BASENAME}_R2_001.fastq.gz"

    echo -e "\n=============================="
    echo "Processing: $BASENAME"
    echo -e "=============================="

    # Align to reference genome (output BAM directly, no .sam)
    bwa mem -t 8 "$REF_GENOME" "$R1_FILE" "$R2_FILE" | \
      samtools view -bS - > "$INTERMED_DIR/${BASENAME}.bam"

    # Sort BAM
    samtools sort -@ 8 -o "$INTERMED_DIR/${BASENAME}.sorted.bam" "$INTERMED_DIR/${BASENAME}.bam"
    rm "$INTERMED_DIR/${BASENAME}.bam"

    # Index BAM
    samtools index "$INTERMED_DIR/${BASENAME}.sorted.bam"

    # Deduplicate using umi_tools (UMIs are first 8bp of R1)
    umi_tools dedup \
        --extract-umi-method=read_id \
        --umi-separator="_" \
        -I "$INTERMED_DIR/${BASENAME}.sorted.bam" \
        -S "$PREPROC_DIR/${BASENAME}.dedup.bam" \
        --log="$PREPROC_DIR/${BASENAME}.dedup_log.txt"

    # BAM stats
    samtools flagstat "$PREPROC_DIR/${BASENAME}.dedup.bam" > "$PREPROC_DIR/${BASENAME}_stats.txt"

done

# ==============================
# Run pipeline output check
# ==============================

bash check_pipeline_output.sh "$PREPROC_DIR"

# ==============================
echo "Pipeline run complete."
