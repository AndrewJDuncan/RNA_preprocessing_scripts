#!/bin/bash

# Activate environment
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rna-tools

# Set directories
PROJECT_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data"
PREPROC_DIR="${PROJECT_DIR}/preproc"
INTERMED_DIR="${PROJECT_DIR}/intermediary_files"

echo "=============================="
echo "  Pipeline Output Check & Summary"
echo "=============================="
echo

# Header
printf "%-35s %-15s %-15s %-15s %-20s\n" "Sample" "Pre-dedup Reads" "Post-dedup Reads" "Dedup Stats" "File Check"

# Loop through dedup BAM files
for DEDUP_BAM in ${PREPROC_DIR}/*.dedup.bam; do
    SAMPLE=$(basename "$DEDUP_BAM" .dedup.bam)
    ALIGN_BAM="${INTERMED_DIR}/${SAMPLE}.sorted.bam"
    STATS_FILE="${PREPROC_DIR}/${SAMPLE}_dedup_stats.txt"

    FILE_STATUS="✔️"
    # Check files exist and non-zero
    for FILE in "$DEDUP_BAM" "$ALIGN_BAM" "$STATS_FILE"; do
        if [[ ! -s "$FILE" ]]; then
            FILE_STATUS="❌"
            break
        fi
    done

    # Get read counts if files exist
    if [[ -s "$DEDUP_BAM" && -s "$ALIGN_BAM" ]]; then
        PRE_DEDUP_READS=$(samtools flagstat "$ALIGN_BAM" | head -n 1 | awk '{print $1}')
        POST_DEDUP_READS=$(samtools flagstat "$DEDUP_BAM" | head -n 1 | awk '{print $1}')
    else
        PRE_DEDUP_READS="NA"
        POST_DEDUP_READS="NA"
    fi

    # Check dedup stats presence
    if [[ -s "$STATS_FILE" ]]; then
        STATS_STATUS="✔️"
    else
        STATS_STATUS="❌"
    fi

    # Print summary line
    printf "%-35s %-15s %-15s %-15s %-20s\n" "$SAMPLE" "$PRE_DEDUP_READS" "$POST_DEDUP_READS" "$STATS_STATUS" "$FILE_STATUS"
done

echo
echo "=============================="
echo "  Check & Summary complete!"
echo "=============================="
