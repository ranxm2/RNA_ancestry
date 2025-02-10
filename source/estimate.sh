#!/bin/bash

echo "Start Processing at: $(date)"
over_all_start_time=$(date +%s)

cd /projects/compbio/users/xran2/Courtney/HISAT2_all/
mkdir -p logs/02-Estimete-ADMIXTURE_mutil-loop/

# get the suid list from ../ref/vcf_id.csv
echo "Processing File index: ${SLURM_ARRAY_TASK_ID}"
input_file="00-Ref/RNASeq_link_all.csv"
line=$(tail -n +2 "$input_file" | sed -n "${SLURM_ARRAY_TASK_ID}p")
IFS=, read -r SUID _ base_name_raw <<< "$line"
BASE_NAME=$(echo "$base_name_raw" | tr -d '\r' | xargs)

echo "SUID: ${SUID}"

mkdir -p result_mutil/06-Summary_result
mkdir -p result_mutil/05-Estimate-ADMIXTURE-loop
mkdir -p result_mutil/05-Estimate-ADMIXTURE-loop/${SUID}

# save the slurm id and the suid to the csv file,
echo "${SLURM_ARRAY_JOB_ID},${SLURM_ARRAY_TASK_ID},${SUID}" >>  result_mutil/06-Summary_result/suid_slurm_id-loop.csv

cd result_mutil/05-Estimate-ADMIXTURE-loop/${SUID}
cp ../../04-VCF/${SUID}.vcf .

export PATH=$PATH:/projects/jmschil/

echo "Processing: SUID: ${SUID}"
echo " " 
echo "-------------------------------------------------"
echo "------- Step 1: Convert VCF to BED/BIM/FAM ------"
echo "-------------------------------------------------"
echo "Step 1 start at: $(date)"
start_time=$(date +%s)


# Convert VCF to BED/BIM/FAM files, save to demo
plink --vcf ${SUID}.vcf  \
      --set-missing-var-ids @:#[b38]\$1\$2  \
      --allow-extra-chr   \
      --make-bed \
      --double-id \
      --out 01-${SUID}

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Step 1 finish at: $(date)"
echo "Step 1 used:  $hours hour(s), $minutes minute(s), and $seconds second(s)."

# change the sex of the samples
sed -i 's/ 0 -9$/ 2 -9/' 01-${SUID}.fam

# check the number of SNPs
echo "Number of SNPs in ${SUID} :$(cat 01-${SUID}.bim | wc -l)"

# delete the chr in the bim file
awk '{ sub(/^chr/, "", $2); print }' 01-${SUID}.bim > 01-${SUID}_fixed.bim
mv 01-${SUID}_fixed.bim 01-${SUID}.bim


echo " " 
echo "-----------------------------------------------------------------"
echo "--- Step 2: Exctract samples SNP overlap with 1000Genome ref ----"
echo "-----------------------------------------------------------------"
echo "Step 2 start at: $(date)"
start_time=$(date +%s)

# exctract the samples from the reference dataset
plink  --bfile 01-${SUID} \
      --extract /projects/jmschil/RNAseq_all/ref/reference_wgs.bim \
      --make-bed \
      --allow-extra-chr \
      --chr 1-22 \
      --out 02-extracted
echo "Number of sample SNPs in extracted dataset: $(cat 02-extracted.bim | wc -l)"

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Step 2 finish at: $(date)"
echo "Step 2 used:  $hours hour(s), $minutes minute(s), and $seconds second(s)."


# add "./" to the family ID
awk '{ $1="./"$1; print }' 02-extracted.fam > 02-extracted_updated.fam
mv 02-extracted_updated.fam 02-extracted.fam


echo " " 
echo "------------------------------------------------------"
echo "--- Step 3: Exctract ref SNP overlap with samples ----"
echo "------------------------------------------------------"
echo "Step 3 start at: $(date)"
start_time=$(date +%s)

# extra the reference samples
plink  --bfile /projects/jmschil/RNAseq_all/ref/reference_wgs \
      --extract 02-extracted.bim \
      --make-bed \
      --allow-extra-chr \
      --chr 1-22 \
      --out 03-REF-extracted

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

echo "Step 3 finish at: $(date)"
echo "Step 3 used:  $hours hour(s), $minutes minute(s), and $seconds second(s)."

# echo the number of SNPs in the extracted dataset
echo "Number of SNPs in reference extracted dataset: $(cat 03-REF-extracted.bim | wc -l)"


echo " " 
echo "-------------------------------------------------"
echo "------- Step 4: Merge samples and reference -----"
echo "-------------------------------------------------"
echo "Step 4 start at: $(date)"
start_time=$(date +%s)


# use the merged reference dataset to extract the samples from the demo dataset
plink --bfile 02-extracted \
      --extract 03-REF-extracted.bim \
        --make-bed \
        --allow-extra-chr \
        --out 04-final_extracted

plink --bfile 03-REF-extracted \
      -bmerge 04-final_extracted.bed 04-final_extracted.bim 04-final_extracted.fam \
        --make-bed \
        --out 05-final_data

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Step 4 finish at: $(date)"
echo "Step 4 used:  $hours hour(s), $minutes minute(s), and $seconds second(s)."

cut -d' ' -f2 05-final_data.fam > fam_ids.txt

sort fam_ids.txt > sorted_fam_ids.txt
sort -k1,1 /projects/jmschil/RNAseq_all/ref/1000G_key.txt > sorted_1000G_key.txt
join -1 1 -2 1 sorted_fam_ids.txt sorted_1000G_key.txt > merged.txt
cut -d' ' -f2 merged.txt > 05-final_super.pop
cut -d' ' -f3 merged.txt > 05-final_sub.pop
sed -i 's/"//g' 05-final_super.pop
sed -i 's/"//g' 05-final_sub.pop
sed -i '1s/^/\n/' 05-final_super.pop
sed -i '1s/^/\n/' 05-final_sub.pop

# Run admixture
export PATH=$PATH:/projects/compbio/users/xran2/Courtney/RNAseq_all/tool

echo " " 
echo "-------------------------------------------------"
echo "------- Step 5: QC with Different MAF and HWE ---"
echo "-------------------------------------------------"

# Loop through different QC parameter combinations
for maf in 0.05 0.01 0.001; do
  for hwe in 0.05 0.01 0.001; do
      echo " " 
      echo "-------------------------------------------------"
      echo "-----     QC with MAF: $maf and HWE: $hwe    ----"
      echo "-------------------------------------------------"
      echo " " 
      echo "Start QC at: $(date) with MAF: $maf, HWE: $hwe"
      start_time=$(date +%s)

      plink --bfile 05-final_data \
            --make-bed \
            --maf $maf \
            --hwe $hwe \
            --allow-extra-chr \
            --chr 1-22 \
            --out final_data_MAF${maf}_HWE${hwe}
      
      echo "MAF: $maf, HWE: $hwe - Number of SNPs: $(cat final_data_MAF${maf}_HWE${hwe}.bim | wc -l)" >> QC_summary.txt

      end_time=$(date +%s)
      elapsed_time=$((end_time - start_time))
      hours=$((elapsed_time / 3600))
      minutes=$(((elapsed_time % 3600) / 60))
      seconds=$((elapsed_time % 60))

      echo "Finish QC at: $(date) with MAF: $maf, HWE: $hwe"
      echo "QC used:  $hours hour(s), $minutes minute(s), and $seconds second(s)."
      echo " " 
      echo " "       
      echo " "       
      echo "Start ADMIXTURE 4 at: $(date) with MAF: $maf, HWE: $hwe"
      start_time=$(date +%s)

      cp 05-final_sub.pop final_data_MAF${maf}_HWE${hwe}.pop
      admixture --cv final_data_MAF${maf}_HWE${hwe}.bed --supervised 4 | tee log_MAF${maf}_HWE${hwe}_4_sup.out
      
      end_time=$(date +%s)
      elapsed_time=$((end_time - start_time))
      hours=$((elapsed_time / 3600))
      minutes=$(((elapsed_time % 3600) / 60))
      seconds=$((elapsed_time % 60))
      echo "Finish ADMIXTURE 4 at: $(date) with MAF: $maf, HWE: $hwe"
      echo "ADMIXTURE 4 used:  $hours hour(s), $minutes minute(s), and $seconds second(s)."
      echo " "       
      echo " "       
      echo " " 
      echo "Start ADMIXTURE 20 at: $(date) with MAF: $maf, HWE: $hwe"
      start_time=$(date +%s)

      cp 05-final_super.pop final_data_MAF${maf}_HWE${hwe}.pop
      admixture --cv final_data_MAF${maf}_HWE${hwe}.bed --supervised 20 | tee log_MAF${maf}_HWE${hwe}_20_sup.out

      end_time=$(date +%s)
      elapsed_time=$((end_time - start_time))
      hours=$((elapsed_time / 3600))
      minutes=$(((elapsed_time % 3600) / 60))
      seconds=$((elapsed_time % 60))
      echo "Finish ADMIXTURE 20 at: $(date) with MAF: $maf, HWE: $hwe"
      echo "ADMIXTURE 20 used:  $hours hour(s), $minutes minute(s), and $seconds second(s)."
      echo " "       
      echo " "       
      echo " "       
      echo " " 
  done
done

over_all_end_time=$(date +%s)
elapsed_time=$((over_all_end_time - over_all_start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
      
echo " "       
echo " "       
echo " " 
echo "End Processing at: $(date)"
echo "Total used:  $hours hour(s), $minutes minute(s), and $seconds second(s)."
