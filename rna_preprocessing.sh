#!/bin/bash

# Ensure Conda is initialised
source ~/miniforge3/etc/profile.d/conda.sh

# Activate the desired environment
conda activate rna-tools

# Define base directory
BASE_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data"

# Define subdirectories
RAW_DIR="$BASE_DIR/rawdata"
PREPROC_DIR="$BASE_DIR/preproc"
INTERMEDIATE_DIR="$BASE_DIR/intermediary_files"
LOG_DIR="$PREPROC_DIR/logs"
REF_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/references"

# Create necessary directories
mkdir -p "$PREPROC_DIR" "$INTERMEDIATE_DIR" "$LOG_DIR"

# Tool paths
BBMAP="/raid/VIDRL-USERS/HOME/aduncan/bbmap/bbmap.sh"
FASTP="fastp"
UMI_TOOLS="umi_tools"
PEAR="pear"
SEQKIT="$CONDA_PREFIX/bin/seqkit"

# Process each sample
for fq1 in "$RAW_DIR"/*R1_001.fastq.gz; do
    sample=$(basename "$fq1" _R1_001.fastq.gz)
    fq2="$RAW_DIR/${sample}_R2_001.fastq.gz"

    echo "========================="
    echo "Processing sample: $sample"
    echo "  R1: $fq1"
    echo "  R2: $fq2"

    # PhiX removal
    echo "Running PhiX removal for $sample"
    $BBMAP ref="$REF_DIR/phix174_ill.ref.fa.gz" \
        in1="$fq1" in2="$fq2" \
        out1="$INTERMEDIATE_DIR/${sample}_noPhiX_R1.fastq.gz" \
        out2="$INTERMEDIATE_DIR/${sample}_noPhiX_R2.fastq.gz" \
        1>"$LOG_DIR/${sample}_phix.out" 2>"$LOG_DIR/${sample}_phix.err"

    # UMI-based deduplication
    echo "Running UMI deduplication for $sample"
    $UMI_TOOLS dedup --paired \
        --stdin="$INTERMEDIATE_DIR/${sample}_noPhiX_R1.fastq.gz" \
        --stdout="$INTERMEDIATE_DIR/${sample}_dedup_R1.fastq.gz" \
        --stdin2="$INTERMEDIATE_DIR/${sample}_noPhiX_R2.fastq.gz" \
        --stdout2="$INTERMEDIATE_DIR/${sample}_dedup_R2.fastq.gz" \
        1>"$LOG_DIR/${sample}_umi.out" 2>"$LOG_DIR/${sample}_umi.err"

    # rRNA filtering
    echo "Running rRNA filtering for $sample"
    $BBMAP ref="$REF_DIR/hg38_rRNA.fa.gz" \
        in1="$INTERMEDIATE_DIR/${sample}_dedup_R1.fastq.gz" \
        in2="$INTERMEDIATE_DIR/${sample}_dedup_R2.fastq.gz" \
        out1="$INTERMEDIATE_DIR/${sample}_rRNAclean_R1.fastq.gz" \
        out2="$INTERMEDIATE_DIR/${sample}_rRNAclean_R2.fastq.gz" \
        stats="$LOG_DIR/${sample}_rrna_stats.txt" \
        1>"$LOG_DIR/${sample}_rrna.out" 2>"$LOG_DIR/${sample}_rrna.err"

    # Merge paired-end reads
    echo "Merging paired reads for $sample"
    $PEAR -f "$INTERMEDIATE_DIR/${sample}_rRNAclean_R1.fastq.gz" \
        -r "$INTERMEDIATE_DIR/${sample}_rRNAclean_R2.fastq.gz" \
        -o "$INTERMEDIATE_DIR/${sample}_merged" \
        1>"$LOG_DIR/${sample}_pear.out" 2>"$LOG_DIR/${sample}_pear.err"

    # FASTP trimming
    MERGED_ASSEMBLED="$INTERMEDIATE_DIR/${sample}_merged.assembled.fastq"
    CLEAN_OUTPUT="$PREPROC_DIR/${sample}_clean.fastq.gz"
    if [[ -s "$MERGED_ASSEMBLED" ]]; then
        echo "Running FASTP trimming for $sample"
        $FASTP -i "$MERGED_ASSEMBLED" -o "$CLEAN_OUTPUT" \
            --qualified_quality_phred 20 --trim_poly_g --detect_adapter_for_pe \
            1>"$LOG_DIR/${sample}_fastp.out" 2>"$LOG_DIR/${sample}_fastp.err"

        # Generate stats
        echo "Generating stats for $sample"
        $SEQKIT stats "$CLEAN_OUTPUT" -j 4 | tee "$PREPROC_DIR/${sample}_stats.json"
    else
        echo "Warning: No merged reads for $sample — skipping FASTP trimming"
    fi

    echo "Completed processing for $sample"
    echo "-------------------------"
done

echo "========================="
echo "Pipeline run completed!"
echo "========================="
