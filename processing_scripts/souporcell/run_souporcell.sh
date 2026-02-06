#!/bin/sh
#SBATCH --partition=compute
#SBATCH --signal=USR2
#SBATCH --ntasks=1
#SBATCH --mem=25GB
#SBATCH --time=1-0:00:00
#SBATCH --cpus-per-task=16
#SBATCH --output=/home/%u/rstudio/souporcell.job.%j
#SBATCH --job-name=run_souporcell         # Job name
#SBATCH --mail-type=ALL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=rmudunuri@tgen.org

source /home/rmudunuri/.bashrc
module load singularity
module load SAMtools

echo "Starting Batch1"

cd /scratch/rmudunuri/souporcell/pedcart/1k_reads_filt/batch1/

singularity exec /scratch/rmudunuri/souporcell/Demuxafy.sif souporcell_pipeline.py \
-i possorted_genome_bam.bam \
-b barcodes.tsv -f /scratch/rmudunuri/souporcell/pedcart/V5_IL13OP/fasta/genome.fa -t 16 -o souporcell_data_test -k 3 \
--known_genotypes /scratch/rmudunuri/souporcell/mpileup/filtered_vcfs/batch11_maf0.05_biallelic_no_missing_AF_1k.vcf \
--known_genotypes_sample_names  PEDCAR_0003_1_PB_Whole_C1 PEDCAR_0004_1_PB_Whole_C1  PEDCAR_0006_1_PB_Whole_C1 

singularity exec /scratch/rmudunuri/souporcell/Demuxafy.sif bash souporcell_summary.sh souporcell_data_test/clusters.tsv > souporcell_data_test/souporcell_summary.tsv

echo "Finished Batch1"
