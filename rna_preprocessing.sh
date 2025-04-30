#!/bin/bash

# Activate Conda and environment
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rna-tools

# Define directories
BASE_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data"
RAW_DIR="$BASE_DIR/rawdata"
PREPROC_DIR="$BASE_DIR/preproc"
INTERMED_DIR="$BASE_DIR/intermediary_files"
LOG_DIR="$PREPROC_DIR/logs"
REF_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/references"
REF_GENOME="$REF_DIR/hg38.fa"

# Create directories if missing
mkdir -p "$PREPROC_DIR" "$INTERMED_DIR" "$LOG_DIR"

# Tools
BBMAP="/raid/VIDRL-USERS/HOME/aduncan/bbmap/bbmap.sh"
SAMTOOLS="samtools"
UMI_TOOLS="umi_tools"

# Process each sample
for fq1 in "$RAW_DIR"/*R1_001.fastq.gz; do
  sample=$(basename "$fq1" _R1_001.fastq.gz)
  fq2="$RAW_DIR/${sample}_R2_001.fastq.gz"

  echo "============================="
  echo "Processing sample: $sample"
  echo "============================="

  echo "Aligning reads to reference..."
  $BBMAP ref="$REF_GENOME" \
         in1="$fq1" in2="$fq2" \
         outm="$INTERMED_DIR/${sample}.sam" \
         1>"$LOG_DIR/${sample}_align.out" 2>"$LOG_DIR/${sample}_align.err"

  echo "Converting SAM to sorted BAM..."
  $SAMTOOLS view -bS "$INTERMED_DIR/${sample}.sam" | \
  $SAMTOOLS sort -o "$INTERMED_DIR/${sample}.sorted.bam"
  $SAMTOOLS index "$INTERMED_DIR/${sample}.sorted.bam"

  echo "Deduplicating with umi_tools..."
  $UMI_TOOLS dedup -I "$INTERMED_DIR/${sample}.sorted.bam" \
                   -S "$PREPROC_DIR/${sample}.dedup.bam" \
                   --log="$LOG_DIR/${sample}_dedup.log"

  echo "Generating BAM stats..."
  $SAMTOOLS flagstat "$PREPROC_DIR/${sample}.dedup.bam" > "$PREPROC_DIR/${sample}_stats.txt"

  echo "Finished processing $sample"
  echo "-----------------------------"
done

echo "========================="
echo "Pipeline run completed!"
echo "========================="
