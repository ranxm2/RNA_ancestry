#!/bin/bash

tool_path="/projects/compbio/users/xran2/tool/STAR_tool"
cd $tool_path

# -----------------------------------------#
# ----- STEP 1: Download STAR package -----#
# -----------------------------------------#
wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz
tar -xzf 2.7.11b.tar.gz
cd STAR-2.7.11b

cd source
make STAR

# add current path to PATH
export PATH=$PATH:/projects/compbio/users/xran2/tool/STAR_tool/STAR-2.7.11b/source


# -------------------------------------------#
# ----- STEP 2: Download reference data -----#
# -------------------------------------------#

# download annotation file
# Link: https://www.ensembl.org/Homo_sapiens/Info/Index?redirect=no
cd $tool_path
wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
gunzip Homo_sapiens.GRCh38.113.gtf.gz

# download fa file
# Link: https://www.ensembl.org/Homo_sapiens/Info/Index?redirect=no
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz


mkdir -p index

STAR --runThreadN 6 \
     --runMode genomeGenerate \
     --genomeDir index \
     --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --sjdbGTFfile Homo_sapiens.GRCh38.113.gtf \
     --sjdbOverhang 100 


