#!/bin/bash

# Ensure Conda is initialised
source ~/miniforge3/etc/profile.d/conda.sh

# Activate the desired environment
conda activate rna-tools

# Define base and sub-directories
BASE_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data"
RAW_DIR="$BASE_DIR/rawdata"
PREPROC_DIR="$BASE_DIR/preproc"
LOG_DIR="$PREPROC_DIR/logs"
OUTPUT_DIR="$BASE_DIR/output"
REF_DIR="$BASE_DIR/references"

# BBMap reference paths
PHIX_REF="/raid/VIDRL-USERS/HOME/aduncan/bbmap/resources/phix174_ill.ref.fa.gz"
RRNA_REF="$REF_DIR/rrna.fasta"

# Tool paths
BBMAP="/raid/VIDRL-USERS/HOME/aduncan/bbmap/bbmap.sh"
FASTP="fastp"
UMI_TOOLS="umi_tools"
PEAR="pear"
SEQKIT="seqkit"

# Create necessary directories
mkdir -p "$PREPROC_DIR" "$LOG_DIR" "$OUTPUT_DIR" "$REF_DIR"

echo "Pipeline starting..."

# Process each sample
for fq1 in "$RAW_DIR"/*R1_001.fastq.gz; do
    sample=$(basename "$fq1" _R1_001.fastq.gz)
    fq2="$RAW_DIR/${sample}_R2_001.fastq.gz"

    echo "Processing sample: $sample"

    # PhiX removal
    echo "Running PhiX removal for $sample"
    $BBMAP ref="$PHIX_REF" in="$fq1" in2="$fq2" \
        out="$PREPROC_DIR/${sample}_noPhiX_R1.fastq.gz" \
        out2="$PREPROC_DIR/${sample}_noPhiX_R2.fastq.gz" \
        1>"$LOG_DIR/${sample}_phix.out" 2>"$LOG_DIR/${sample}_phix.err"

    if [[ ! -s "$PREPROC_DIR/${sample}_noPhiX_R1.fastq.gz" ]]; then
        echo "  PhiX removal failed for $sample — skipping"
        continue
    fi

    # UMI-based deduplication
    echo "Running UMI deduplication for $sample"
    $UMI_TOOLS dedup --paired \
        --stdin="$PREPROC_DIR/${sample}_noPhiX_R1.fastq.gz" \
        --stdout="$PREPROC_DIR/${sample}_dedup_R1.fastq.gz" \
        1>"$LOG_DIR/${sample}_umi.out" 2>"$LOG_DIR/${sample}_umi.err"

    cp "$PREPROC_DIR/${sample}_noPhiX_R2.fastq.gz" "$PREPROC_DIR/${sample}_dedup_R2.fastq.gz"

    # rRNA filtering
    echo "Running rRNA filtering for $sample"
    $BBMAP ref="$RRNA_REF" in="$PREPROC_DIR/${sample}_dedup_R1.fastq.gz" \
        in2="$PREPROC_DIR/${sample}_dedup_R2.fastq.gz" \
        outu="$PREPROC_DIR/${sample}_rRNAclean_R1.fastq.gz" \
        outu2="$PREPROC_DIR/${sample}_rRNAclean_R2.fastq.gz" \
        stats="$PREPROC_DIR/${sample}_rrna_stats.txt" \
        1>"$LOG_DIR/${sample}_rrna.out" 2>"$LOG_DIR/${sample}_rrna.err"

    # Merge paired reads
    echo "Merging paired reads for $sample"
    $PEAR -f "$PREPROC_DIR/${sample}_rRNAclean_R1.fastq.gz" \
        -r "$PREPROC_DIR/${sample}_rRNAclean_R2.fastq.gz" \
        -o "$PREPROC_DIR/${sample}_merged" \
        1>"$LOG_DIR/${sample}_pear.out" 2>"$LOG_DIR/${sample}_pear.err"

    MERGED_FASTQ="$PREPROC_DIR/${sample}_merged.assembled.fastq"
    if [[ ! -s "$MERGED_FASTQ" ]]; then
        echo "  No merged reads for $sample — skipping trimming"
        continue
    fi

    # FASTP trimming
    echo "Running FASTP trimming for $sample"
    $FASTP -i "$MERGED_FASTQ" -o "$PREPROC_DIR/${sample}_clean.fastq.gz" \
        --qualified_quality_phred 20 --trim_poly_g --detect_adapter_for_pe \
        1>"$LOG_DIR/${sample}_fastp.out" 2>"$LOG_DIR/${sample}_fastp.err"

    # SeqKit stats
    echo "Generating stats for $sample"
    $SEQKIT stats "$PREPROC_DIR/${sample}_clean.fastq.gz" -j 4 | tee "$PREPROC_DIR/${sample}_stats.json"

    echo "Completed processing for $sample"
done

echo "========================="
echo "Pipeline run completed!"
echo "========================="
