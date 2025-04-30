#!/bin/bash

PREPROC_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/preproc"

echo "========================="
echo "Checking pipeline outputs"
echo "========================="

for fq in ${PREPROC_DIR}/*_clean.fastq.gz; do
    [ -e "$fq" ] || continue
    sample=$(basename "$fq")
    echo "Sample: $sample"

    # Check if gzipped fastq is valid
    gzip -t "$fq" && echo "  FASTQ OK" || echo "  FASTQ CORRUPTED"

    # Count reads if OK
    if gzip -t "$fq" >/dev/null 2>&1; then
        readcount=$(zcat "$fq" | echo $((`wc -l`/4)))
        echo "  Read count: $readcount"
    fi

    echo "-------------------------"
done

# Check JSON reports
for json in ${PREPROC_DIR}/*.json; do
    [ -e "$json" ] || continue
    sample=$(basename "$json")
    echo "Sample: $sample"
    total_reads=$(grep -m1 '"total_reads"' "$json" | head -1 | awk -F ':' '{print $2}' | tr -d ', ')
    echo "  total_reads: $total_reads"
    echo "-------------------------"
done

echo "========================="
echo " Check completed "
echo "========================="
