#!/bin/bash
set -euo pipefail

# Define base directory
BASE_DIR="$HOME/projects/rna_pipeline/mgp_test_data"

# Define subdirectories
RAW_DIR="$BASE_DIR/rawdata"
PREPROC_DIR="$BASE_DIR/preproc"
LOG_DIR="$PREPROC_DIR/logs"
OUTPUT_DIR="$BASE_DIR/output"
REF_DIR="$BASE_DIR/references"

# Create directories if they don't exist
mkdir -p "$PREPROC_DIR" "$LOG_DIR" "$OUTPUT_DIR" "$REF_DIR"

# Tool paths
BBMAP="bbmap.sh"
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
echo "Running PhiX removal for sample: $sample"
$BBMAP filter ref="$REF_DIR/phix174_ill.ref.fa" \
    in1="$fq1" in2="$fq2" \
    out1="$PREPROC_DIR/${sample}_noPhiX_R1.fastq.gz" \
    out2="$PREPROC_DIR/${sample}_noPhiX_R2.fastq.gz" \
    1>"$LOG_DIR/${sample}_phix.out" 2>"$LOG_DIR/${sample}_phix.err"
echo "Finished PhiX removal for sample: $sample"

# UMI-based deduplication
echo "Running UMI deduplication for sample: $sample"
$UMI_TOOLS dedup --paired \
    --in1="$PREPROC_DIR/${sample}_noPhiX_R1.fastq.gz" \
    --in2="$PREPROC_DIR/${sample}_noPhiX_R2.fastq.gz" \
    --out1="$PREPROC_DIR/${sample}_dedup_R1.fastq.gz" \
    --out2="$PREPROC_DIR/${sample}_dedup_R2.fastq.gz" \
    1>"$LOG_DIR/${sample}_umi.out" 2>"$LOG_DIR/${sample}_umi.err"
echo "Finished UMI deduplication for sample: $sample"

    # rRNA count
    $BBMAP filter ref="$REF_DIR/rrna.fasta" \
        in1="$PREPROC_DIR/${sample}_dedup_R1.fastq.gz" \
        in2="$PREPROC_DIR/${sample}_dedup_R2.fastq.gz" \
        stats="$PREPROC_DIR/${sample}_rrna.txt" \
        1>"$LOG_DIR/${sample}_rrna.out" 2>"$LOG_DIR/${sample}_rrna.err"

    # Merge paired-end reads
    $PEAR -f "$PREPROC_DIR/${sample}_dedup_R1.fastq.gz" \
          -r "$PREPROC_DIR/${sample}_dedup_R2.fastq.gz" \
          -o "$PREPROC_DIR/${sample}_merged" \
          1>"$LOG_DIR/${sample}_pear.out" 2>"$LOG_DIR/${sample}_pear.err"

    # Adapter trimming, polyA/T trimming, Q20 trimming
    FASTP_INPUT="$PREPROC_DIR/${sample}_merged.assembled.fastq"
    FASTP_OUTPUT="$PREPROC_DIR/${sample}_clean.fastq.gz"

    $FASTP -i "$FASTP_INPUT" -o "$FASTP_OUTPUT" \
        --qualified_quality_phred 20 --trim_poly_g --detect_adapter_for_pe \
        1>"$LOG_DIR/${sample}_fastp.out" 2>"$LOG_DIR/${sample}_fastp.err"

    # Generate statistics
    $SEQKIT stats "$FASTP_OUTPUT" -j 4 | tee "$PREPROC_DIR/${sample}_stats.json"
done
