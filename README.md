# Tools required
```
conda create -N RNA-seq

conda install -n RNA-seq -c bioconda fastqc
conda install -n RNA-seq -c bioconda fastp
conda install -n RNA-seq -c bioconda multiqc
conda install -n RNA-seq -c bioconda star
conda install -n RNA-seq -c bioconda samtools
conda install -n RNA-seq -c bioconda deeptools
conda install -n RNA-seq -c bioconda salmon

#For differential expression using DESeq2
conda create -N DEseq2 r-essentials r-base

conda install -N DEseq2 -c bioconda bioconductor-deseq2
conda install -N DEseq2 -c bioconda bioconductor-tximport 
conda install -N DEseq2 -c r r-ggplot2
```

Using hg38 from UCSC (indexed with hiseq2)

# Steps / tools used
1. Raw sequence QC - **FastQC**
2. UMI extraction - **UMI-tools**
3. Trimming - removal of low-quality reads and adapter sequences **FastP** ?trim galore
4. Repeat **FastQC** post-trim
5. ? use BBsplit for genome contaminant removal
6. rRNA count? - BBsplit?
7. Summary stats generation
8. STAR alignment --> then take these to the transcriptome and quantify them.  transcriptome on ensembl -- ftp://ftp.ensembl.org/pub/release-90/fasta/bos_taurus/
     Faster if straight to salmon - psuedoaligner

# Code used at each step

Initial QC
```
fastqc <sample>_1.fastq.gz -d . -o .
fastqc <sample>_2.fastq.gz -d . -o .
```


  --> combine these into one report with MultiQC


Then can get read count with:
```
totalreads=$(unzip -c <sample>_fastqc.zip <sample>_fastqc/fastqc_data.txt | grep 'Total Sequences' | cut -f 2)
echo $totalreads
```
  --> saves totalreads as a variable for later in pipeline


Then trimming...
```
fastp -i <sample>_R1.fastq.gz -I <sample>_R2.fastq.gz -o <sample>_R1.trimmed.fastq.gz -O <sample>_R2.trimmed.fastq.gz --detect_adapter_for_pe -XX 25 -j <sample>.fastp.json -h <sample>.fastp.html
# change XX to min read count. Clarify this w Eike...
```
  --> generates html report


Repeat fastqc for trimmed data
```
fastqc <sample>_R1.trimmed.fastq.gz -d . -o .
fastqc <sample>_R2.trimmed.fastq.gz -d . -o .
```


Alignment: 
```
GENOMEDIR=/.../references/hg38.fa
#need to write this bit
```



# Directory setup:
/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/rawdata/ 
→ raw FASTQs

/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/preproc/ 
→ final BAM/json/logs

/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/intermediary_files/ 
→ intermediate BAMs, fastqs

/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/references/hg38.fa 
→ reference genome


# References
- UCD bioinformatics core git and RNA-seq course
- nf-core/rnaseq https://github.com/nf-core/rnaseq
- Cebola Lab https://github.com/CebolaLab/RNA-seq
