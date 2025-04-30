#!/bin/bash

# Check pipeline outputs for each sample

BASE_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data"
PREPROC_DIR="$BASE_DIR/preproc"
LOG_DIR="$PREPROC_DIR/logs"

echo "========================="
echo " Checking pipeline output "
echo "========================="

# Loop through deduplicated BAM files
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
