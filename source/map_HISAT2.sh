#!/bin/bash

# Directories
FASTQ_DIR="./data"
RESULT_DIR="result"
FA_DIR="result/01-FA"
SAM_DIR="result/02-SAM"
SAM_SORT_DIR="result/03-SAM_SORT"
VCF_DIR="result/04-VCF"

mkdir -p ${RESULT_DIR} ${FA_DIR} ${SAM_DIR} ${SAM_SORT_DIR} ${VCF_DIR} 

# test for input
# BASE_NAME="demo"
echo "Processing File: ${BASE_NAME}"

# Define the FastQ files (R1 and R2)
fq1=${FASTQ_DIR}/${BASE_NAME}_R1_001.fastq.gz
fq2=${FASTQ_DIR}/${BASE_NAME}_R2_001.fastq.gz
SUID="S001"

echo "SUID: ${SUID}"
echo "Processing: ${BASE_NAME}"
echo "FastQ files R1: ${fq1}"
echo "FastQ files R2: ${fq2}"

echo " " 
echo "-------------------------------------------------"
echo "------- Step 1: Convert FQ files to FA ----------"
echo "-------------------------------------------------"
echo "Step 1 start at: $(date)"
start_time=$(date +%s)
seqkit fq2fa ${fq1} -o ${FA_DIR}/${SUID}_1.fa.gz
seqkit fq2fa ${fq2} -o ${FA_DIR}/${SUID}_2.fa.gz
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Step 1 finish at: $(date)"
echo "Step 1 used:  $hours hour(s), $minutes minute(s), and $seconds second(s)."

# HISAT2 index path



echo " " 
echo "-------------------------------------------------"
echo "------- Step 2: HISAT2 Mapping Mutil Mapping ----"
echo "-------------------------------------------------"
echo "Step 2 start at: $(date)"
start_time=$(date +%s)

HISAT2_INDEX="tool/grch38"
hisat2 -f -x ${HISAT2_INDEX} \
        -1 ${FA_DIR}/${SUID}_1.fa.gz \
        -2 ${FA_DIR}/${SUID}_2.fa.gz \
        -S ${SAM_DIR}/${SUID}.sam \
                                
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Step 2 finish at: $(date)"
echo "Step 2 used:  $hours hour(s), $minutes minute(s), and $seconds second(s)."

# delect the fa files
rm ${FA_DIR}/${SUID}_1.fa.gz
rm ${FA_DIR}/${SUID}_2.fa.gz

echo " " 
echo "-------------------------------------------------"
echo "------- Step 3: Sorting with SAMTOOLS ------------"
echo "-------------------------------------------------"
echo "Step 3 start at: $(date)"
start_time=$(date +%s)

samtools sort ${SAM_DIR}/${SUID}.sam \
         -o ${SAM_SORT_DIR}/${SUID}_sort.sam
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

echo "Step 3 finish at: $(date)"
echo "Step 3 used:  $hours hour(s), $minutes minute(s), and $seconds second(s)."

# delect the sam files
rm ${SAM_DIR}/${SUID}.sam

echo " " 
echo "-------------------------------------------------"
echo "------- Step 4: VCF file ------------------------"
echo "-------------------------------------------------"
echo "Step 4 start at: $(date)"
start_time=$(date +%s)
REFERENCE_GENOME="/projects/jmschil/RNAseq_all/00-Index/GRCh38.fna"


bcftools mpileup \
    -f ${REFERENCE_GENOME} \
    ${SAM_SORT_DIR}/${SUID}_sort.sam \
    --max-depth 8000 | bcftools call \
    -vm > ${VCF_DIR}/${SUID}.vcf

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Step 4 finish at: $(date)"
echo "Step 4 used:  $hours hour(s), $minutes minute(s), and $seconds second(s)."

# delect the sam_sort files
rm ${SAM_SORT_DIR}/${SUID}_sort.sam



