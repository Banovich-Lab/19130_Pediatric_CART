!/bin/sh
#SBATCH --partition=compute
#SBATCH --signal=USR2
#SBATCH --ntasks=1
#SBATCH --mem=400GB
#SBATCH --time=7-0:00:00
#SBATCH --cpus-per-task=16
#SBATCH --output=/home/%u/cramtobam.job.%j
#SBATCH --job-name=cramtobam                  # Job name
#SBATCH --mail-type=ALL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=rmudunuri@tgen.org

source /home/rmudunuri/.bashrc
module load singularity
module load SAMtools

echo "Starting Batch1"

cd /scratch/rmudunuri/souporcell/pedcart/1k_reads_filt/batch1/

samtools view -b -T /scratch/rmudunuri/souporcell/pedcart/V5_IL13OP/fasta/genome.fa \
-o possorted_genome_bam.bam possorted_genome_bam.cram

echo "Finished Batch1"