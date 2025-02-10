#!/bin/bash


tool_path="/projects/compbio/users/xran2/tool/bcftool"

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
export PATH=$PATH:${tool_path}/bcftools-1.21

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
