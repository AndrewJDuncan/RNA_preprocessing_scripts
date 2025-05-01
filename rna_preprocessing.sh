#!/bin/bash

# Stop on error
set -e

# Load conda environment system
source /raid/VIDRL-USERS/HOME/aduncan/miniforge3/etc/profile.d/conda.sh
conda activate rna-tools

# Dry-run mode flag
DRYRUN=false
if [[ "$1" == "--dry-run" ]]; then
    DRYRUN=true
    echo "🔍 Dry-run mode enabled — no commands will be executed."
fi

# Set directories
RAW_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/rawdata"
INTER_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/intermediary_files"
PREPROC_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/preproc"
REF_GENOME="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/references/hg38_index"

# Validate directories exist
echo "✅ Checking directories..."
for DIR in "$RAW_DIR" "$INTER_DIR" "$PREPROC_DIR"; do
    if [ ! -d "$DIR" ]; then
        echo "❌ Missing directory: $DIR"
        exit 1
    fi
done
echo "✅ All directories found."

# Confirm tools available
echo "✅ Checking tool paths..."
echo "umi_tools: $(which umi_tools)"
echo "hisat2:    $(which hisat2)"
echo "samtools:  $(which samtools)"
echo ""

# Make sure intermediary and output directories exist
mkdir -p "$INTER_DIR"
mkdir -p "$PREPROC_DIR"

# Loop through R1 files in raw data directory
for R1_FILE in "$RAW_DIR"/*_R1_001.fastq.gz; do
    BASENAME=$(basename "$R1_FILE")
    SAMPLE=$(echo "$BASENAME" | sed 's/_R1_001.fastq.gz//')
    R2_FILE="$RAW_DIR/${SAMPLE}_R2_001.fastq.gz"

    echo "============================"
    echo "Preparing sample: $SAMPLE"
    echo "============================"

    # Check input files exist
    if [[ ! -f "$R1_FILE" ]] || [[ ! -f "$R2_FILE" ]]; then
        echo "❌ Missing input files for sample $SAMPLE"
        continue
    fi

    if $DRYRUN; then
        echo "✅ Would extract UMIs from:"
        echo "   $R1_FILE and $R2_FILE"
        echo "✅ Would align with hisat2 using reference: $REF_GENOME"
        echo "✅ Would convert SAM to sorted BAM"
        echo "✅ Would index the sorted BAM"
        echo "✅ Would deduplicate BAM"
        echo "✅ Would generate flagstat stats"
    else
        # Step 1: Extract UMIs
        umi_tools extract \
            --extract-method=string \
            --bc-pattern=NNNNNNNN \
            --stdin="$R1_FILE" \
            --stdout="$INTER_DIR/${SAMPLE}_R1_extracted.fastq.gz" \
            --read2-in="$R2_FILE" \
            --read2-out="$INTER_DIR/${SAMPLE}_R2_extracted.fastq.gz" \
            --log="$INTER_DIR/${SAMPLE}_extract.log"

        # Step 2: Align with HISAT2
        hisat2 -x "$REF_GENOME" \
            -1 "$INTER_DIR/${SAMPLE}_R1_extracted.fastq.gz" \
            -2 "$INTER_DIR/${SAMPLE}_R2_extracted.fastq.gz" \
            -S "$INTER_DIR/${SAMPLE}.sam" \
            > "$PREPROC_DIR/${SAMPLE}_hisat2.out" \
            2> "$PREPROC_DIR/${SAMPLE}_hisat2.err"

        # Step 3: Convert SAM to sorted BAM
        samtools view -bS "$INTER_DIR/${SAMPLE}.sam" | samtools sort -o "$INTER_DIR/${SAMPLE}.sorted.bam"

        # Step 3.5: Index the sorted BAM so umi_tools can use it
        samtools index "$INTER_DIR/${SAMPLE}.sorted.bam"

        # Remove SAM file to save space
        rm "$INTER_DIR/${SAMPLE}.sam"

        # Step 4: Deduplicate using UMI-tools
        umi_tools dedup \
            --extract-umi-method=read_id \
            -I "$INTER_DIR/${SAMPLE}.sorted.bam" \
            -S "$PREPROC_DIR/${SAMPLE}.dedup.bam" \
            --log="$PREPROC_DIR/${SAMPLE}_dedup.log"

        # Step 5: Generate BAM statistics
        samtools flagstat "$PREPROC_DIR/${SAMPLE}.dedup.bam" > "$PREPROC_DIR/${SAMPLE}_dedup_stats.txt"

        echo "✅ Sample $SAMPLE complete."
    fi

    echo ""
done

echo "============================"
if $DRYRUN; then
    echo "✅ Dry-run completed successfully!"
else
    echo "✅ Pipeline run completed!"
fi
echo "============================"
