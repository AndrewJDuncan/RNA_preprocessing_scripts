#!/bin/bash

# Ensure Conda is initialised
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rna-tools

# Set variables
PROJECT_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data"
REF_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/references"
PREPROC_DIR="${PROJECT_DIR}/preproc"
INTERMEDIATE_DIR="${PROJECT_DIR}/intermediary_files"
LOG_DIR="${PREPROC_DIR}/logs"

# Create necessary directories if they don't exist
mkdir -p "$PREPROC_DIR"
mkdir -p "$INTERMEDIATE_DIR"
mkdir -p "$LOG_DIR"

# Process each sample
for R1 in "${PROJECT_DIR}"/rawdata/*_R1_001.fastq.gz; do
    SAMPLE=$(basename "$R1" | cut -d_ -f1-4)
    R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"

    echo "========================================="
    echo "Processing sample: $SAMPLE"
    echo "R1: $R1"
    echo "R2: $R2"
    echo "========================================="

    ## UMI Deduplication
    echo "Running UMI deduplication for $SAMPLE"
    umi_tools dedup --paired --stdin="$R1" --stdout="${INTERMEDIATE_DIR}/${SAMPLE}_dedup_R1.fastq.gz"
    umi_tools dedup --paired --stdin="$R2" --stdout="${INTERMEDIATE_DIR}/${SAMPLE}_dedup_R2.fastq.gz"

    ## rRNA Filtering (skipping if no reference)
    if [[ -f "${REF_DIR}/hg38_rRNA_mask.fa" ]]; then
        echo "Running rRNA filtering for $SAMPLE"
        bbmap.sh ref="${REF_DIR}/hg38_rRNA_mask.fa" \
            in1="${INTERMEDIATE_DIR}/${SAMPLE}_dedup_R1.fastq.gz" \
            in2="${INTERMEDIATE_DIR}/${SAMPLE}_dedup_R2.fastq.gz" \
            outu1="${INTERMEDIATE_DIR}/${SAMPLE}_rRNAclean_R1.fastq.gz" \
            outu2="${INTERMEDIATE_DIR}/${SAMPLE}_rRNAclean_R2.fastq.gz" \
            threads=4 stats="${INTERMEDIATE_DIR}/${SAMPLE}_rRNA_stats.txt"
    else
        echo "Skipping rRNA filtering for $SAMPLE (reference missing)"
        cp "${INTERMEDIATE_DIR}/${SAMPLE}_dedup_R1.fastq.gz" "${INTERMEDIATE_DIR}/${SAMPLE}_rRNAclean_R1.fastq.gz"
        cp "${INTERMEDIATE_DIR}/${SAMPLE}_dedup_R2.fastq.gz" "${INTERMEDIATE_DIR}/${SAMPLE}_rRNAclean_R2.fastq.gz"
    fi

    ## Paired-End Merging
    echo "Merging paired reads for $SAMPLE"
    pear -f "${INTERMEDIATE_DIR}/${SAMPLE}_rRNAclean_R1.fastq.gz" \
         -r "${INTERMEDIATE_DIR}/${SAMPLE}_rRNAclean_R2.fastq.gz" \
         -o "${INTERMEDIATE_DIR}/${SAMPLE}_merged" -j 4 > "${LOG_DIR}/${SAMPLE}_pear.out" 2> "${LOG_DIR}/${SAMPLE}_pear.err"

    ## Quality Trimming with FASTP (if assembled file exists)
    if [[ -f "${INTERMEDIATE_DIR}/${SAMPLE}_merged.assembled.fastq" ]]; then
        echo "Running FASTP trimming for $SAMPLE"
        fastp -i "${INTERMEDIATE_DIR}/${SAMPLE}_merged.assembled.fastq" \
              -o "${PREPROC_DIR}/${SAMPLE}_clean.fastq.gz" \
              --json "${PREPROC_DIR}/${SAMPLE}.json" \
              --html "${PREPROC_DIR}/${SAMPLE}.html" \
              --thread 4 > "${LOG_DIR}/${SAMPLE}_fastp.out" 2> "${LOG_DIR}/${SAMPLE}_fastp.err"
    else
        echo "Warning: no merged reads for $SAMPLE - skipping trimming"
    fi

    echo "Completed processing for $SAMPLE"
    echo "-----------------------------------------"

done

echo "========================="
echo "Pipeline run completed!"
echo "========================="
