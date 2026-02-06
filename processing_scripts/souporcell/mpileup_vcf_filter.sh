#!/bin/sh
#SBATCH --signal=USR2
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --time=1:00:00
#SBATCH --output=/home/%u/rstudio/vcf_filter.job.%j
#SBATCH --job-name=VCFfilter                  # Job name
#SBATCH --mail-type=ALL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=rmudunuri@tgen.org

source /home/rmudunuri/.bashrc
conda activate demux
module load BCFtools


BATCH=1 # this is the VCF batch name
echo "Starting batch${BATCH}"
bcftools filter \
	-R /scratch/rmudunuri/souporcell/mpileup/filtered_bed_outputs/b1_mpileup3k.bed \ 
	-Oz -o /scratch/rmudunuri/souporcell/mpileup/filtered_vcfs/batch${BATCH}_maf0.05_biallelic_no_missing_AF_3k.vcf.gz \
	/scratch/aoill/projects/CAR-T/00_new_2025/new_data/vcfs/batch${BATCH}_maf0.05_biallelic_no_missing_AF.vcf.gz
