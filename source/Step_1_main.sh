#!/bin/bash

# create a new directory for RNA_map
mkdir -p tool
cd tool

# #-------------------------------------------#
# #                                           #
# #      Method 1: use conda to install       #
# #                                           #
# #-------------------------------------------#

# # Check current conda env
# conda env list

# # create a new env called rna_map
# conda create --name rna_map
# conda activate rna_map
# conda install -c bioconda hisat2 subread 
# conda install -c bioconda samtools

# # Check for installation 
# hisat2 --version
# featureCounts -v
# samtools -v

#-------------------------------------------#
#                                           #
#      Method 2: install local              #
#                                           #
#-------------------------------------------#

# Download HISAT2 
wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download
unzip download
#add /projects/compbio/users/xran2/tool/RNA_map/hisat2-2.2.1 to PATH
export PATH="$(pwd)/hisat2-2.2.1:$PATH"
hisat2 --version


# Download featureCounts
wget https://github.com/ShiLab-Bioinformatics/subread/releases/download/2.0.2/subread-2.0.2-Linux-x86_64.tar.gz
tar -xzvf subread-2.0.2-Linux-x86_64.tar.gz
export PATH="$(pwd)/subread-2.0.2-Linux-x86_64/bin:$PATH"
featureCounts -v

# Download samtools
wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2
tar -xjvf samtools-1.20.tar.bz2
export PATH="$(pwd)/samtools-1.20:$PATH"
samtools 

wget https://github.com/shenwei356/seqkit/releases/latest/download/seqkit_linux_amd64.tar.gz
tar -xvzf seqkit_linux_amd64.tar.gz.1
chmod +x seqkit
export PATH="$(pwd):$PATH"



#-------------------------------------------#
#                                           #
#      Part 2: Make reference               #
#                                           #
#-------------------------------------------#

# Download Human reference data
# Human reference (GRCh38) 
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xzvf grch38_genome.tar.gz

cd grch38
./make_grch38.sh
hisat2-build genome.fa genome
cd ..


# Download Human reference gtf file
# Version 106
wget https://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz  # Download GTF file version 106
gunzip Homo_sapiens.GRCh38.106.gtf.gz

# # Version 113
# wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz # Download GTF file version 113
# gunzip Homo_sapiens.GRCh38.113.gtf.gz

# Create Gene ID and Gene Name Mapping file
awk '$3 == "gene" {print $10 "\t" $14}' Homo_sapiens.GRCh38.106.gtf | sed 's/"//g' | sed 's/;//g' | sort > gene_id_name_map_sorted.txt
# awk '$3 == "gene" {print $10 "\t" $14}' Homo_sapiens.GRCh38.113.gtf | sed 's/"//g' | sed 's/;//g' | sort > gene_id_name_map_sorted.txt


cd ..




tool_path="./tool/bcftool"

mkdir -p $tool_path
cd $tool_path


#-------------------------------------------#
#                                           #
#      Part1 : install tool                 #
#                                           #
#-------------------------------------------#

# Download bcftools
wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2
tar -xjvf bcftools-1.21.tar.bz2

# Install bcftools
cd bcftools-1.21
./configure
make

# export current path to PATH
export PATH="$(pwd):$PATH"

# Check for installation
bcftools --version


#-------------------------------------------#
#                                           #
#      Part 2: Make reference               #
#                                           #
#-------------------------------------------#


cd ../

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz
gunzip GCF_000001405.26_GRCh38_genomic.fna.gz

sed -E 's/^(>NC_[0-9]+\.[0-9]+) .*chromosome ([0-9XYM]+).*/>\2/' \
    GCF_000001405.26_GRCh38_genomic.fna \
    > simplified_GRCh38.fna


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
SUID="S001" # this is the sample ID, user can change it to any name
echo "Processing File: ${BASE_NAME}"

# Define the FastQ files (R1 and R2)
fq1=${FASTQ_DIR}/${BASE_NAME}_R1_001.fastq.gz
fq2=${FASTQ_DIR}/${BASE_NAME}_R2_001.fastq.gz


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





echo " " 
echo "-------------------------------------------------"
echo "------- Step 2: HISAT2 Mapping Mutil Mapping ----"
echo "-------------------------------------------------"
echo "Step 2 start at: $(date)"
start_time=$(date +%s)

# HISAT2 index path
HISAT2_INDEX="./tool/grch38/genome"
hisat2 -f -x ${HISAT2_INDEX} \
        -1 ${FA_DIR}/${SUID}_1.fa.gz \
        -2 ${FA_DIR}/${SUID}_2.fa.gz \
        -S ${SAM_DIR}/${SUID}.sam 
                                
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
REFERENCE_GENOME="./tool/bcftool/GCF_000001405.26_GRCh38_genomic.fna"

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




