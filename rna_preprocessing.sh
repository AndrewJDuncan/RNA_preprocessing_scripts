#!/bin/bash

# ===============================
# RNA Preprocessing Pipeline
# Recommended: Run with:
# nohup bash rna_preprocessing.sh > pipeline_run.log 2>&1 &
# or within a screen session.
# ===============================

# Load conda environment
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rna-tools

# Define directories
BASE_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data"
RAW_DIR="$BASE_DIR/rawdata"
PREPROC_DIR="$BASE_DIR/preproc"
INTER_DIR="$BASE_DIR/intermediary_files"
REF_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/references"
LOG_DIR="$PREPROC_DIR/logs"

# Reference files
GENOME_REF="$REF_DIR/hg38.fa"

# Create output directories if missing
mkdir -p "$PREPROC_DIR" "$INTER_DIR" "$LOG_DIR"

# Process each pair of FASTQ files
for R1_FILE in "$RAW_DIR"/*_R1_001.fastq.gz; do
    SAMPLE=$(basename "$R1_FILE" | sed 's/_R1_001.fastq.gz//')
    R2_FILE="$RAW_DIR/${SAMPLE}_R2_001.fastq.gz"

    echo "========================="
    echo "Processing sample: $SAMPLE"
    echo "R1: $R1_FILE"
    echo "R2: $R2_FILE"
    echo "========================="

    # 1. Align reads to reference and generate SAM
    echo "Aligning reads for $SAMPLE"
    /raid/VIDRL-USERS/HOME/aduncan/bbmap/bbmap.sh ref="$GENOME_REF" in1="$R1_FILE" in2="$R2_FILE" \
        out="$INTER_DIR/$SAMPLE.sam" nodisk overwrite=t

    # 2. Convert SAM to sorted BAM, index BAM for umi_tools dedup
    echo "Converting to BAM for $SAMPLE"
    samtools sort -o "$INTER_DIR/$SAMPLE.sorted.bam" "$INTER_DIR/$SAMPLE.sam"
    samtools index "$INTERMED_DIR/${SAMPLE}.sorted.bam"
    rm "$INTER_DIR/$SAMPLE.sam"

    # 3. Deduplicate using UMI-tools
    echo "Deduplicating for $SAMPLE"
    umi_tools dedup -I "$INTER_DIR/$SAMPLE.sorted.bam" -S "$PREPROC_DIR/$SAMPLE.dedup.bam" --extract-umi-method=read_id --umi-separator=":" --log="$PREPROC_DIR/$SAMPLE.dedup_log.txt"

    # 4. Collect stats
    echo "Generating BAM stats for $SAMPLE"
    samtools flagstat "$PREPROC_DIR/$SAMPLE.dedup.bam" > "$PREPROC_DIR/$SAMPLE.stats.txt"

    # 5. Separator line
    echo "-------------------------"

done

echo "========================="
echo "Pipeline run completed!"
echo "========================="
echo "    "


# Run output check

echo "========================="
echo " Checking pipeline output "
echo "========================="

PREPROC_DIR="$BASE_DIR/preproc"

for BAM_FILE in "$PREPROC_DIR"/*.dedup.bam; do
    [ -e "$BAM_FILE" ] || continue  # skip if none exist

    SAMPLE=$(basename "$BAM_FILE" | sed 's/.dedup.bam//')

    echo "Sample: $SAMPLE"

    # Check BAM file size
    if [ ! -s "$BAM_FILE" ]; then
        echo " ❌ dedup.bam missing or empty!"
    else
        echo " ✅ dedup.bam present"
    fi

    # Check stats file
    if [ ! -s "$PREPROC_DIR/$SAMPLE.stats.txt" ]; then
        echo " ❌ stats.txt missing or empty!"
    else
        echo " ✅ stats.txt present"
    fi

    # Check dedup log
    if [ ! -s "$PREPROC_DIR/$SAMPLE.dedup_log.txt" ]; then
        echo " ⚠️  dedup_log.txt missing or empty"
    else
        echo " ✅ dedup_log.txt present"
    fi

    echo "-------------------------"
done

echo "========================="
echo " Output check complete! "
echo "========================="


