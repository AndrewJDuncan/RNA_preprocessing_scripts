#!/bin/bash

# Set directories
RAW_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/rawdata"
INTERMED_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/intermediary_files"
PREPROC_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/preproc"
REF_GENOME="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/references/hg38.fa"

# Load environment if needed
source activate rna-tools

# Create output directories if not exist
mkdir -p "$INTERMED_DIR" "$PREPROC_DIR"

# List samples
samples=$(ls "$RAW_DIR"/*.bam | xargs -n 1 basename | sed 's/.bam//')

# Process each sample
for SAMPLE in $samples; do
  echo "=============================================================="
  echo "Processing sample: $SAMPLE"
  echo "=============================================================="

  # Sort BAM
  echo "Sorting BAM for $SAMPLE..."
  samtools sort -o "$INTERMED_DIR/${SAMPLE}.sorted.bam" "$RAW_DIR/${SAMPLE}.bam"

  if [[ ! -s "$INTERMED_DIR/${SAMPLE}.sorted.bam" ]]; then
    echo "ERROR: sorted BAM missing for $SAMPLE"
    exit 1
  fi

  # Index BAM
  echo "Indexing BAM for $SAMPLE..."
  samtools index "$INTERMED_DIR/${SAMPLE}.sorted.bam"

  if [[ ! -s "$INTERMED_DIR/${SAMPLE}.sorted.bam.bai" ]]; then
    echo "ERROR: BAM index missing for $SAMPLE"
    exit 1
  fi

  # Deduplicate with umi_tools
  echo "Deduplicating BAM for $SAMPLE..."
  umi_tools dedup \
    -I "$INTERMED_DIR/${SAMPLE}.sorted.bam" \
    -S "$PREPROC_DIR/${SAMPLE}.dedup.bam" \
    --log="$PREPROC_DIR/${SAMPLE}.dedup_log.txt"

  if [[ ! -s "$PREPROC_DIR/${SAMPLE}.dedup.bam" ]]; then
    echo "ERROR: deduplicated BAM missing for $SAMPLE"
    exit 1
  fi

  # BAM stats
  echo "Generating stats for $SAMPLE..."
  samtools flagstat "$PREPROC_DIR/${SAMPLE}.dedup.bam" > "$PREPROC_DIR/${SAMPLE}.dedup_stats.txt"

  if [[ ! -s "$PREPROC_DIR/${SAMPLE}.dedup_stats.txt" ]]; then
    echo "ERROR: stats.txt missing for $SAMPLE"
    exit 1
  fi

  echo "Finished processing $SAMPLE"
done

# Run pipeline output check script
echo "=============================================================="
echo "Pipeline run completed!"
echo "Running check_pipeline_output.sh..."
echo "=============================================================="

bash check_pipeline_output.sh

echo " ----------------------- "
echo "~~~ All done ~~~"
echo " ----------------------- "
