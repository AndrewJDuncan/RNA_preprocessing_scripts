#!/bin/bash

# Ensure Conda is initialised
source ~/miniforge3/etc/profile.d/conda.sh

# Activate the desired environment
conda activate rna-tools

# Set variables
PROJECT_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data"
REF_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/references"
BBMAP_DIR="/raid/VIDRL-USERS/HOME/aduncan/bbmap"
preproc_dir="$PROJECT_DIR/preproc"
intermediary_dir="$PROJECT_DIR/intermediary_files"

mkdir -p "$preproc_dir"
mkdir -p "$intermediary_dir"

# Process each sample
for r1 in "$PROJECT_DIR/rawdata"/*_R1_001.fastq.gz; do
    sample=$(basename "$r1" | sed 's/_R1_001.fastq.gz//')
    r2="$PROJECT_DIR/rawdata/${sample}_R2_001.fastq.gz"

    echo "============================="
    echo "Processing sample: $sample"
    echo "R1: $r1"
    echo "R2: $r2"
    echo "============================="

    # PhiX removal skipped for now

    # UMI extraction
    echo "Extracting UMIs for $sample"
    umi_tools extract --bc-pattern=NNNNNNNN \
        --stdin="$r1" \
        --stdout="$intermediary_dir/${sample}_R1_extracted.fastq.gz" \
        --read2-in="$r2" \
        --read2-out="$intermediary_dir/${sample}_R2_extracted.fastq.gz"

    # Align with BBMap
    echo "Aligning reads for $sample"
    "$BBMAP_DIR/bbmap.sh" ref="$REF_DIR/hg38.fa" \
        in1="$intermediary_dir/${sample}_R1_extracted.fastq.gz" \
        in2="$intermediary_dir/${sample}_R2_extracted.fastq.gz" \
        out="$intermediary_dir/${sample}.sam"

    # Convert SAM to BAM
    samtools view -bS "$intermediary_dir/${sample}.sam" > "$intermediary_dir/${sample}.bam"

    # Sort BAM
    samtools sort "$intermediary_dir/${sample}.bam" -o "$intermediary_dir/${sample}.sorted.bam"

    # Index BAM
    samtools index "$intermediary_dir/${sample}.sorted.bam"

    # Deduplicate
    echo "Deduplicating BAM for $sample"
    umi_tools dedup -I "$intermediary_dir/${sample}.sorted.bam" -S "$preproc_dir/${sample}.dedup.bam"

    # BAM stats
    samtools flagstat "$preproc_dir/${sample}.dedup.bam" > "$preproc_dir/${sample}.dedup_stats.txt"

    echo "Finished processing $sample"
    echo "------------------------------------------------------"

done

echo "========================="
echo "Pipeline run completed!"
echo "========================="
