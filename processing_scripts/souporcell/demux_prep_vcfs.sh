#!/bin/sh
#SBATCH --signal=USR2
#SBATCH --ntasks=1
#SBATCH --mem=80GB
#SBATCH --time=2-0:00:00
#SBATCH --output=/home/%u/projects/CAR-T/00_new_2025/logs/prep_vcfs_demx.job.%j
#SBATCH --job-name=prep_vcfs_demx                  # Job name
#SBATCH --mail-type=NONE     # Mail events (NONE, BEGIN, END, FAIL, ALL)


source activate demuxlet
module load BCFtools
module load VCFtools


VCF_PATH=/scratch/aoill/projects/CAR-T/00_new_2025/new_data/vcfs
# 514
PEDCAR_0002_514=PEDCAR_0002_1_PB_Whole_C1_TETE2.bwa.deepvariant.pass.db.vcf.gz
# 574
PEDCAR_0003_574=PEDCAR_0003_1_PB_Whole_C1_TETE2.bwa.deepvariant.pass.db.vcf.gz
# 689
PEDCAR_0004_689=PEDCAR_0004_1_PB_Whole_C1_TETE2.bwa.deepvariant.pass.db.vcf.gz 
# 692
PEDCAR_0005_692=PEDCAR_0005_1_PB_Whole_C1_TETE2.bwa.deepvariant.pass.db.vcf.gz
# 705
PEDCAR_0006_705=PEDCAR_0006_1_PB_Whole_C1_TETE2.bwa.deepvariant.pass.db.vcf.gz
# 716 
PEDCAR_0007_716=PEDCAR_0007_1_PB_Whole_C1_TETE2.bwa.deepvariant.pass.db.vcf.gz 


################################################################################
# Batch 1
# UPNs: 514, 689, 705
################################################################################
BATCH=1
echo "Batch $BATCH" 

# Merge samples into one vcf
bcftools merge $VCF_PATH/$PEDCAR_0002_514 $VCF_PATH/$PEDCAR_0004_689 $VCF_PATH/$PEDCAR_0006_705 -Oz -o $VCF_PATH/batch$BATCH.vcf.gz

# Get number of sites
echo "number snps before filtering" 
zcat $VCF_PATH/batch$BATCH.vcf.gz | grep -v "#" | wc -l
# 163750

# This isnt working (output name missing batch number)...add brackets
# Keep sites with MAF >0.05, bialleleic sites, and no missing data
vcftools --gzvcf $VCF_PATH/batch${BATCH}.vcf.gz \
--maf 0.05 \
--min-alleles 2 --max-alleles 2 \
--max-missing 1.0 \
--recode \
--stdout | bgzip -c > $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz

#zcat $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz | grep -v "#" | wc -l
# 29376

# calculate allele frequency and add to vcf
bcftools +fill-tags -Oz -o $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz \
  $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz \
  -- -t AF
tabix -p vcf $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz

echo "number snps after filtering" 
zcat $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz  | grep -v "#" | wc -l
# 29376

################################################################################
# Batch 2
# UPNs: 574, 689, 716
################################################################################
BATCH=2
echo "Batch $BATCH" 

# Merge samples into one vcf
bcftools merge $VCF_PATH/$PEDCAR_0003_574 $VCF_PATH/$PEDCAR_0004_689 $VCF_PATH/$PEDCAR_0007_716 -Oz -o $VCF_PATH/batch$BATCH.vcf.gz

# Get number of sites
zcat $VCF_PATH/batch$BATCH.vcf.gz | grep -v "#" | wc -l
# 166952

# Keep sites with MAF >0.05, bialleleic sites, and no missing data
vcftools --gzvcf $VCF_PATH/batch${BATCH}.vcf.gz \
--maf 0.05 \
--min-alleles 2 --max-alleles 2 \
--max-missing 1.0 \
--recode \
--stdout | bgzip -c > $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz

#zcat $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz | grep -v "#" | wc -l


# calculate allele frequency and add to vcf
bcftools +fill-tags -Oz -o $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz \
  $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz \
  -- -t AF
tabix -p vcf $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz

zcat $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz  | grep -v "#" | wc -l
# 29189


################################################################################
# Batch 3
# UPNs: 689, 692, 705
################################################################################
BATCH=3

# Merge samples into one vcf
bcftools merge $VCF_PATH/$PEDCAR_0004_689 $VCF_PATH/$PEDCAR_0005_692 $VCF_PATH/$PEDCAR_0006_705 -Oz -o $VCF_PATH/batch$BATCH.vcf.gz

# Get number of sites
zcat $VCF_PATH/batch$BATCH.vcf.gz | grep -v "#" | wc -l
# 162856

# Keep sites with MAF >0.05, bialleleic sites, and no missing data
vcftools --gzvcf $VCF_PATH/batch${BATCH}.vcf.gz \
--maf 0.05 \
--min-alleles 2 --max-alleles 2 \
--max-missing 1.0 \
--recode \
--stdout | bgzip -c > $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz

#zcat $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz | grep -v "#" | wc -l


# calculate allele frequency and add to vcf
bcftools +fill-tags -Oz -o $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz \
  $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz \
  -- -t AF
tabix -p vcf $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz

zcat $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz  | grep -v "#" | wc -l
# 30554


################################################################################
# Batch 4
# UPNs: 514, 689, 716
################################################################################
BATCH=4

# Merge samples into one vcf
bcftools merge $VCF_PATH/$PEDCAR_0002_514 $VCF_PATH/$PEDCAR_0004_689 $VCF_PATH/$PEDCAR_0007_716 -Oz -o $VCF_PATH/batch$BATCH.vcf.gz

# Get number of sites
zcat $VCF_PATH/batch$BATCH.vcf.gz | grep -v "#" | wc -l
# 164599

# Keep sites with MAF >0.05, bialleleic sites, and no missing data
vcftools --gzvcf $VCF_PATH/batch${BATCH}.vcf.gz \
--maf 0.05 \
--min-alleles 2 --max-alleles 2 \
--max-missing 1.0 \
--recode \
--stdout | bgzip -c > $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz

#zcat $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz | grep -v "#" | wc -l


# calculate allele frequency and add to vcf
bcftools +fill-tags -Oz -o $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz \
  $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz \
  -- -t AF
tabix -p vcf $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz

zcat $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz  | grep -v "#" | wc -l
#28694

################################################################################
# Batch 5
# UPNs: 574, 689, 692
################################################################################
BATCH=5

# Merge samples into one vcf
bcftools merge $VCF_PATH/$PEDCAR_0003_574 $VCF_PATH/$PEDCAR_0004_689 $VCF_PATH/$PEDCAR_0005_692 -Oz -o $VCF_PATH/batch$BATCH.vcf.gz

# Get number of sites
zcat $VCF_PATH/batch$BATCH.vcf.gz | grep -v "#" | wc -l
# 166275

# Keep sites with MAF >0.05, bialleleic sites, and no missing data
vcftools --gzvcf $VCF_PATH/batch${BATCH}.vcf.gz \
--maf 0.05 \
--min-alleles 2 --max-alleles 2 \
--max-missing 1.0 \
--recode \
--stdout | bgzip -c > $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz

#zcat $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz | grep -v "#" | wc -l


# calculate allele frequency and add to vcf
bcftools +fill-tags -Oz -o $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz \
  $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz \
  -- -t AF
tabix -p vcf $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz

zcat $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz  | grep -v "#" | wc -l
# 28965


################################################################################
# Batch 6
# UPNs: 705, 716
################################################################################
BATCH=6

# Merge samples into one vcf
bcftools merge $VCF_PATH/$PEDCAR_0006_705 $VCF_PATH/$PEDCAR_0007_716 -Oz -o $VCF_PATH/batch$BATCH.vcf.gz

# Get number of sites
zcat $VCF_PATH/batch$BATCH.vcf.gz | grep -v "#" | wc -l
# 140622

# Keep sites with MAF >0.05, bialleleic sites, and no missing data
vcftools --gzvcf $VCF_PATH/batch${BATCH}.vcf.gz \
--maf 0.05 \
--min-alleles 2 --max-alleles 2 \
--max-missing 1.0 \
--recode \
--stdout | bgzip -c > $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz

#zcat $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz | grep -v "#" | wc -l


# calculate allele frequency and add to vcf
bcftools +fill-tags -Oz -o $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz \
  $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz \
  -- -t AF
tabix -p vcf $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz

zcat $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz  | grep -v "#" | wc -l
# 37823


################################################################################
# Batch 7
# UPNs: 514, 689, 705, 716
################################################################################
BATCH=7

# Merge samples into one vcf
bcftools merge $VCF_PATH/$PEDCAR_0002_514 $VCF_PATH/$PEDCAR_0004_689 $VCF_PATH/$PEDCAR_0006_705 $VCF_PATH/$PEDCAR_0007_716 -Oz -o $VCF_PATH/batch$BATCH.vcf.gz

# Get number of sites
zcat $VCF_PATH/batch$BATCH.vcf.gz | grep -v "#" | wc -l
# 183721

# Keep sites with MAF >0.05, bialleleic sites, and no missing data
vcftools --gzvcf $VCF_PATH/batch${BATCH}.vcf.gz \
--maf 0.05 \
--min-alleles 2 --max-alleles 2 \
--max-missing 1.0 \
--recode \
--stdout | bgzip -c > $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz

#zcat $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz | grep -v "#" | wc -l


# calculate allele frequency and add to vcf
bcftools +fill-tags -Oz -o $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz \
  $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz \
  -- -t AF
tabix -p vcf $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz

zcat $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz  | grep -v "#" | wc -l
# 24950

################################################################################
# Batch 8
# UPNs: 574, 689, 692, 716
################################################################################
BATCH=8

# Merge samples into one vcf
bcftools merge $VCF_PATH/$PEDCAR_0003_574 $VCF_PATH/$PEDCAR_0004_689 $VCF_PATH/$PEDCAR_0005_692 $VCF_PATH/$PEDCAR_0007_716 -Oz -o $VCF_PATH/batch$BATCH.vcf.gz

# Get number of sites
zcat $VCF_PATH/batch$BATCH.vcf.gz | grep -v "#" | wc -l
# 185558

# Keep sites with MAF >0.05, bialleleic sites, and no missing data
vcftools --gzvcf $VCF_PATH/batch${BATCH}.vcf.gz \
--maf 0.05 \
--min-alleles 2 --max-alleles 2 \
--max-missing 1.0 \
--recode \
--stdout | bgzip -c > $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz

#zcat $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz | grep -v "#" | wc -l


# calculate allele frequency and add to vcf
bcftools +fill-tags -Oz -o $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz \
  $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz \
  -- -t AF
tabix -p vcf $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz

zcat $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz  | grep -v "#" | wc -l
# 24726

################################################################################
# Batch 9
# UPNs: 689, 705, 716
################################################################################
BATCH=9

# Merge samples into one vcf
bcftools merge $VCF_PATH/$PEDCAR_0004_689 $VCF_PATH/$PEDCAR_0006_705 $VCF_PATH/$PEDCAR_0007_716 -Oz -o $VCF_PATH/batch$BATCH.vcf.gz

# Get number of sites
zcat $VCF_PATH/batch$BATCH.vcf.gz | grep -v "#" | wc -l
# 165143

# Keep sites with MAF >0.05, bialleleic sites, and no missing data
vcftools --gzvcf $VCF_PATH/batch${BATCH}.vcf.gz \
--maf 0.05 \
--min-alleles 2 --max-alleles 2 \
--max-missing 1.0 \
--recode \
--stdout | bgzip -c > $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz

#zcat $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz | grep -v "#" | wc -l


# calculate allele frequency and add to vcf
bcftools +fill-tags -Oz -o $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz \
  $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz \
  -- -t AF
tabix -p vcf $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz

zcat $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz  | grep -v "#" | wc -l
# 30432


################################################################################
# Batch 10
# UPNs: 514, 689, 692
################################################################################
BATCH=10

# Merge samples into one vcf
bcftools merge $VCF_PATH/$PEDCAR_0002_514 $VCF_PATH/$PEDCAR_0004_689 $VCF_PATH/$PEDCAR_0005_692 -Oz -o $VCF_PATH/batch$BATCH.vcf.gz

# Get number of sites
zcat $VCF_PATH/batch$BATCH.vcf.gz | grep -v "#" | wc -l
# 164192

# Keep sites with MAF >0.05, bialleleic sites, and no missing data
vcftools --gzvcf $VCF_PATH/batch${BATCH}.vcf.gz \
--maf 0.05 \
--min-alleles 2 --max-alleles 2 \
--max-missing 1.0 \
--recode \
--stdout | bgzip -c > $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz

#zcat $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz | grep -v "#" | wc -l


# calculate allele frequency and add to vcf
bcftools +fill-tags -Oz -o $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz \
  $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz \
  -- -t AF
tabix -p vcf $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz

zcat $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz  | grep -v "#" | wc -l
# 28631

################################################################################
# Batch 11
# UPNs: 574, 689, 705
################################################################################
BATCH=11

# Merge samples into one vcf
bcftools merge $VCF_PATH/$PEDCAR_0003_574 $VCF_PATH/$PEDCAR_0004_689 $VCF_PATH/$PEDCAR_0006_705 -Oz -o $VCF_PATH/batch$BATCH.vcf.gz

# Get number of sites
zcat $VCF_PATH/batch$BATCH.vcf.gz | grep -v "#" | wc -l
# 166098

# Keep sites with MAF >0.05, bialleleic sites, and no missing data
vcftools --gzvcf $VCF_PATH/batch${BATCH}.vcf.gz \
--maf 0.05 \
--min-alleles 2 --max-alleles 2 \
--max-missing 1.0 \
--recode \
--stdout | bgzip -c > $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz

#zcat $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz | grep -v "#" | wc -l


# calculate allele frequency and add to vcf
bcftools +fill-tags -Oz -o $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz \
  $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing.vcf.gz \
  -- -t AF
tabix -p vcf $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz

zcat $VCF_PATH/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz  | grep -v "#" | wc -l
# 30246
