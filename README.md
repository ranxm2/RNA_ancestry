# RNA Ancestry Project
This project focuses on using RNA sequence data to infer ancestry.

## Pipeline Overview
![Pipeline Overview](images/overview.jpg)

## Tools Required

## Pipeline Installation
```bash
# Download the pipeline
git clone https://github.com/ranxm2/RNA_ancestry.git
```

After cloning the repository, please run the following code to initialize the tools:

```bash
./source/Initial/setup_all.sh
```

##  Pipeline Input
```bash
fastq_1="demo_R1_001.fastq.gz"
fastq_2="demo_R2_001.fastq.gz"
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

After the variant calling steps, use the flowing code to estimation ancestry:


```bash
./source/estimate.sh
```


Result should be similar like thgis:

```bash
Ancestry       Proportion   N_SNPs    Estimation_Corr (95% CI)
African        75.5%        240001    0.78 (0.68-0.88)
European       13.2%        240001    0.62 (0.23-0.97)
Asian          11.3%        240001    0.32 (0.18-0.39)


```