#!/bin/bash

# Activate environment
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rna-tools

# Set directories
PROJECT_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data"
PREPROC_DIR="${PROJECT_DIR}/preproc"
INTERMED_DIR="${PROJECT_DIR}/intermediary_files"

echo "========================="
echo " Pipeline Summary Report "
echo "========================="
echo

# Header for report
echo -e "Sample\tPre-dedup Reads\tPost-dedup Reads"

# Loop through dedup BAM files
for DEDUP_BAM in ${PREPROC_DIR}/*.dedup.bam; do
    SAMPLE=$(basename "$DEDUP_BAM" .dedup.bam)

    # Corresponding pre-dedup BAM in intermediary_files
    ALIGN_BAM="${INTERMED_DIR}/${SAMPLE}.sorted.bam"

    # Check both files exist
    if [[ -f "$ALIGN_BAM" && -f "$DEDUP_BAM" ]]; then
        PRE_DEDUP_READS=$(samtools flagstat "$ALIGN_BAM" | head -n 1 | awk '{print $1}')
        POST_DEDUP_READS=$(samtools flagstat "$DEDUP_BAM" | head -n 1 | awk '{print $1}')
        echo -e "${SAMPLE}\t${PRE_DEDUP_READS}\t${POST_DEDUP_READS}"
    else
        echo -e "${SAMPLE}\tMISSING FILE(S)"
    fi
done

echo
echo "========================="
echo " Summary complete! "
echo "========================="
