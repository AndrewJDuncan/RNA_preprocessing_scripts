# Steps in preprocessing and tools used
 1. Raw sequence QC - **FastQC** // **MultiQC**

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

# Code used at each step
**FastQC** 

```fastqc <sample>_1.fastq.gz -d . -o .
fastqc <sample>_2.fastq.gz -d . -o .```

--> combine these into one report with MultiQC

Then can get read count with:
```totalreads=$(unzip -c <sample>_fastqc.zip <sample>_fastqc/fastqc_data.txt | grep 'Total Sequences' | cut -f 2)
echo $totalreads```
--> saves totalreads as a variable for later in pipeline





# Directory setup:
/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/rawdata/ 
→ raw FASTQs

/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/preproc/ 
→ final BAM/json/logs

/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/intermediary_files/ 
→ intermediate BAMs, fastqs

/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/references/hg38.fa 
→ reference genome
