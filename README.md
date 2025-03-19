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
```

## Step 1 : Variant Calling
If you want to use HISAT2 as the aligner, please use the following code:

```bash
./source/map_HISAT2.sh ${fastq_1} ${fastq_2}
```

If you want to use STAR as the aligner, please use the following code:

```bash
./source/map_STAR.sh ${fastq_1} ${fastq_2}
```

## Step 2  Ancestry Estimation

After the variant calling steps, user need to specify the reference genome for the ancestry estimation. The code to build reference is in `source/Initial/Setup-HISAT2.sh`.

Affer the reference genome is built, user can use the ADMIXTURE to estimate the ancestry proportion. Details of the ADMIXTURE can be found in the [ADMIXTURE](https://dalexander.github.io/admixture/). 

The code to run the ADMIXTURE is in `source/Ancestry/Ancestry.sh`. User can run the following code to estimate the ancestry proportion:



## Result Output

Result should be similar like thgis:

```bash
Ancestry       Proportion   N_SNPs    
African        75.5%        240001    
European       13.2%        240001    
Asian          11.3%        240001    


```