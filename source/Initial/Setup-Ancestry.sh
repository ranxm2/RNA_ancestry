#!/bin/bash

#-------------------------------------------#
#                                           #
#           Step 1: Install Plink           #
#                                           #
#-------------------------------------------#





mkdir -p ref/1000Gene_ref logs/06-SNP
cd /projects/compbio/users/xran2/Courtney/RNAseq_all/ref

# Run PLINK command for each chromosome
plink --bfile PLINK/ALL.chr${CHR}.GRCh38 --keep PLINK/allref1010_ids.txt --make-bed --out 1000Gene_ref/chr${CHR}

export PATH=$PATH:/projects/jmschil/

# get individuals that we want for our reference panel--loop through 1 to 22 chromosomes with a line similar to this:
# Loop through chromosomes 1 to 22
for i in {1..22}; do
    echo "Processing chromosome ${i}"
    # Get individuals for the reference panel for each chromosome
    plink --bfile PLINK/ALL.chr${i}.GRCh38 --keep allref1010_ids.txt --make-bed --out 1000Gene_ref/chr${i}

    # Name the SNPs in the reference set for each chromosome
    plink --bfile 1000Gene_ref/chr${i} --set-missing-var-ids @:#[b38]\$1\$2 --make-bed --out 1000Gene_ref/chr${i}_named
    echo "Done"
done

# After looping through chromosomes, merge all chromosomes into a whole genome reference
ls 1000Gene_ref/chr{2..22}_named.bed > merge_list.txt
echo "Merging all chromosomes into a whole genome reference"
plink --bfile 1000Gene_ref/chr1_named --merge-list merge_list.txt --make-bed --out reference_wgs
echo "Done"
plink version



