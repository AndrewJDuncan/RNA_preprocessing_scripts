#!/bin/bash

# Ensure Conda is initialised and environment is activated
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rna-tools

# Define base directory
BASE_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data"

# Define subdirectories
RAW_DIR="$BASE_DIR/rawdata"
PREPROC_DIR="$BASE_DIR/preproc"
LOG_DIR="$PREPROC_DIR/logs"
OUTPUT_DIR="$BASE_DIR/output"
REF_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/references"

# Create directories if they don't exist
mkdir -p "$PREPROC_DIR" "$LOG_DIR" "$OUTPUT_DIR"

# Tool paths
BBMAP="/raid/VIDRL-USERS/HOME/aduncan/bbmap/bbmap.sh"
FASTP="fastp"
UMI_TOOLS="umi_tools"
PEAR="pear"
SEQKIT="seqkit"

# Process each sample
for fq1 in "$RAW_DIR"/*R1_001.fastq.gz; do
    sample=$(basename "$fq1" _R1_001.fastq.gz)
    fq2="$RAW_DIR/${sample}_R2_001.fastq.gz"

    echo "Processing sample: $sample"
    echo "  R1: $fq1"
    echo "  R2: $fq2"

    # PhiX removal
    echo "Running PhiX removal for $sample"
    $BBMAP ref="$REF_DIR/hg38_masked.fa" \
        in1="$fq1" in2="$fq2" \
        out1="$PREPROC_DIR/${sample}_noPhiX_R1.fastq.gz" \
        out2="$PREPROC_DIR/${sample}_noPhiX_R2.fastq.gz" \
        1>"$LOG_DIR/${sample}_phix.out" 2>"$LOG_DIR/${sample}_phix.err"

    # UMI deduplication
    echo "Running UMI deduplication for $sample"
    $UMI_TOOLS dedup --paired \
        --in1="$PREPROC_DIR/${sample}_noPhiX_R1.fastq.gz" \
        --in2="$PREPROC_DIR/${sample}_noPhiX_R2.fastq.gz" \
        --out1="$PREPROC_DIR/${sample}_dedup_R1.fastq.gz" \
        --out2="$PREPROC_DIR/${sample}_dedup_R2.fastq.gz" \
        1>"$LOG_DIR/${sample}_umi.out" 2>"$LOG_DIR/${sample}_umi.err"

    # Merge paired reads
    echo "Merging paired reads for $sample"
    $PEAR -f "$PREPROC_DIR/${sample}_dedup_R1.fastq.gz" \
         -r "$PREPROC_DIR/${sample}_dedup_R2.fastq.gz" \
         -o "$PREPROC_DIR/${sample}_merged" \
         1>"$LOG_DIR/${sample}_pear.out" 2>"$LOG_DIR/${sample}_pear.err"

    if [[ -f "$PREPROC_DIR/${sample}_merged.assembled.fastq" ]]; then
        # Adapter trimming, polyA/T trimming, Q20 trimming
        echo "Running FASTP trimming for $sample"
        $FASTP -i "$PREPROC_DIR/${sample}_merged.assembled.fastq" \
               -o "$PREPROC_DIR/${sample}_clean.fastq.gz" \
               --qualified_quality_phred 20 --trim_poly_g --detect_adapter_for_pe \
               1>"$LOG_DIR/${sample}_fastp.out" 2>"$LOG_DIR/${sample}_fastp.err"

        # Generate stats if seqkit is available
        if command -v $SEQKIT &> /dev/null; then
            $SEQKIT stats "$PREPROC_DIR/${sample}_clean.fastq.gz" -j 4 | tee "$PREPROC_DIR/${sample}_stats.json"
        else
            echo "Warning: seqkit not found — skipping stats for $sample"
        fi
    else
        echo "No merged reads for $sample — skipping trimming"
    fi

    echo "Completed processing for $sample"
done

echo "========================="
echo "Pipeline run completed!"
echo "========================="
