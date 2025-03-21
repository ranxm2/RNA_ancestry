# RNA Ancestry Project
This project focuses on using RNA sequence data to infer ancestry.

## Pipeline Overview
![Pipeline Overview](images/overview.jpg)

## Tools Required

## Pipeline Installation
```bash
# Download the pipeline
git clone https://github.com/ranxm2/RNA_ancestry.git
cd RNA_ancestry
```

After cloning the repository, please run the following code to initialize the tools:

```bash
chmod +x ./source/Initial/Setup-HISAT2.sh
bash ./source/Initial/Setup-HISAT2.sh
```

##  Pipeline Input
The pipeline requires two fastq files as input. Please specify the path to the fastq files in the following code:

```bash
fastq_1="data/demo_R1_001.fastq.gz"
fastq_2="data/demo_R2_001.fastq.gz"
ls ${fastq_1} ${fastq_2}
```

## Step 1 : Variant Calling
If you want to use HISAT2 as the aligner, the paired fastq files should be have similar names. The code to run the HISAT2 is in `source/map_HISAT2.sh`. User can run the following code to map the fastq files:

```bash
./source/map_HISAT2.sh "data/demo"
```

If you want to use STAR as the aligner, please use the following code:

```bash
./source/map_STAR.sh  "data/demo"
```

This will give a result of VCF file recording the variants. We use S001 for sample ID. User can modify the sample ID in the code at `source/map_HISAT2.sh`. The result should be in the following format:

```bash
>ls "result/04-VCF"
S001.vcf
```

## Step 2  Ancestry Estimation

After the variant calling steps, user need to specify the reference genome for the ancestry estimation. The code to build reference is in `source/Initial/Setup-Ancestry.sh`. User can select the reference genome from the [1000 Genomes Project](https://www.internationalgenome.org/). User can select the target population and download the reference genome. 

Affer the reference genome is built, user can use the ADMIXTURE to estimate the ancestry proportion. Details of the ADMIXTURE can be found in the [ADMIXTURE](https://dalexander.github.io/admixture/). 

The code to run the ADMIXTURE is in `source/estimate.sh`. User can run the following code to estimate the ancestry proportion:


## Result Output

Result of the estimated result should be in this format:

```bash
Ancestry       Proportion   N_SNPs    
African        75.5%        240001    
European       13.2%        240001    
Asian          11.3%        240001    
```