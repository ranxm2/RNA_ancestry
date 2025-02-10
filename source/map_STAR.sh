#!/bin/bash

cd /projects/compbio/users/xran2/Courtney/STAR_all/

mkdir -p logs logs/01-Pipeline_mutil/

# Directories
FASTQ_DIR="../Raw_data"
RESULT_DIR="result_mutil"
FQ_DIR="result_mutil/01-FQ"
SAM_DIR="result_mutil/02-SAM"
SAM_SORT_DIR="result_mutil/02-SAM_SORT"
VCF_DIR="result_mutil/04-VCF"

STAR_REF_INDEX="/projects/compbio/users/xran2/tool/STAR_tool/index/"
REFERENCE_GENOME="/projects/compbio/users/xran2/tool/bcftool/simplified_GRCh38.fna"

export PATH=$PATH:/projects/compbio/users/xran2/tool/STAR_tool/STAR-2.7.11b/source
export PATH=$PATH:/projects/compbio/users/xran2/tool/bcftool/bcftools-1.21
export PATH=$PATH:/projects/compbio/users/xran2/tool/RNA_map/samtools-1.20

mkdir -p ${RESULT_DIR} ${FQ_DIR} ${SAM_DIR} ${SAM_SORT_DIR} ${VCF_DIR} 

echo "Processing File index: ${SLURM_ARRAY_TASK_ID}"
input_file="00-Ref/RNASeq_link_all.csv"
line=$(tail -n +2 "$input_file" | sed -n "${SLURM_ARRAY_TASK_ID}p")
IFS=, read -r SUID _ base_name_raw <<< "$line"
BASE_NAME=$(echo "$base_name_raw" | tr -d '\r' | xargs)

# Define the FastQ files (R1 and R2)
fq1=${FASTQ_DIR}/${BASE_NAME}_R1_001.fastq.gz
fq2=${FASTQ_DIR}/${BASE_NAME}_R2_001.fastq.gz

echo "SUID: ${SUID}"
echo "Processing: ${BASE_NAME}"
echo "FastQ files R1: ${fq1}"
echo "FastQ files R2: ${fq2}"

echo " " 
echo "-------------------------------------------------"
echo "------- Step 1: Unzip FastQ files ---------------"
echo "-------------------------------------------------"
echo "Step 2 start at: $(date)"
start_time=$(date +%s)

cp ${fq1} ${FQ_DIR}
cp ${fq2} ${FQ_DIR}
gunzip ${FQ_DIR}/${BASE_NAME}_R1_001.fastq.gz
gunzip ${FQ_DIR}/${BASE_NAME}_R2_001.fastq.gz
fq1_unzip=${FQ_DIR}/${BASE_NAME}_R1_001.fastq
fq2_unzip=${FQ_DIR}/${BASE_NAME}_R2_001.fastq
                                
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Step 2 finish at: $(date)"
echo "Step 2 used:  $hours hour(s), $minutes minute(s), and $seconds second(s)."

echo " " 
echo "-------------------------------------------------"
echo "------- Step 2: STAR Mapping Mutil Mapping ----"
echo "-------------------------------------------------"
echo "Step 2 start at: $(date)"
start_time=$(date +%s)

STAR --runThreadN 20 \
     --genomeDir ${STAR_REF_INDEX} \
     --readFilesIn ${fq1_unzip} ${fq2_unzip} \
     --outSAMtype SAM\
     --outFileNamePrefix ./${SAM_DIR}/${SUID}.

                            
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Step 2 finish at: $(date)"
echo "Step 2 used:  $hours hour(s), $minutes minute(s), and $seconds second(s)."

rm ${fq1_unzip} ${fq2_unzip}pan

echo " " 
echo "-------------------------------------------------"
echo "------- Step 3: Sorting with SAMTOOLS ------------"
echo "-------------------------------------------------"
echo "Step 3 start at: $(date)"
start_time=$(date +%s)

samtools sort ${SAM_DIR}/${SUID}.Aligned.out.sam \
        -o ${SAM_SORT_DIR}/${SUID}_sort.sam

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

echo "Step 3 finish at: $(date)"
echo "Step 3 used:  $hours hour(s), $minutes minute(s), and $seconds second(s)."

rm ${SAM_DIR}/${SUID}.Aligned.out.sam


echo " " 
echo "-------------------------------------------------"
echo "------- Step 4: VCF file ------------------------"
echo "-------------------------------------------------"
echo "Step 4 start at: $(date)"
start_time=$(date +%s)

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

rm ${SAM_SORT_DIR}/${SUID}_sort.sam