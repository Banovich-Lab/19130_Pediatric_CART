#!/bin/sh
#SBATCH --signal=USR2
#SBATCH --ntasks=1
#SBATCH --mem=100GB
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=/home/%u/rstudio/mpileup.job.%j
#SBATCH --job-name=mpileup                  # Job name
#SBATCH --mail-type=ALL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=rmudunuri@tgen.org

source /home/rmudunuri/.bashrc
conda activate demux
module load BCFtools
module load SAMtools

FASTA=/scratch/rmudunuri/souporcell/pedcart/V5_IL13OP/fasta/genome.fa


#############################
#Batch 1
#############################
echo "Starting Batch1"
BAM1=/scratch/rmudunuri/souporcell/mpileup/b1_cram.txt
VCF1=/scratch/aoill/projects/CAR-T/00_new_2025/new_data/vcfs/batch1_maf0.05_biallelic_no_missing_AF.vcf.gz
OUT1=/scratch/rmudunuri/souporcell/mpileup/b1_mpileup

echo "BED generation started"
bcftools query -f '%CHROM\t%POS0\t%END\n' $VCF8 > /scratch/rmudunuri/souporcell/mpileup/batch1.bed
echo "BED generation finished"

echo "mpileup started"
samtools mpileup nthreads=8 \
	--bam-list $BAM1 \
	--fasta-ref $FASTA \
	--positions /scratch/rmudunuri/souporcell/mpileup/batch1.bed \
	--output $OUT1
echo "mpileup finished"