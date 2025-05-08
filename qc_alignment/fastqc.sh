#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --array=1-1
#SBATCH --output=logs/fastqc_%A_%a.out
#SBATCH --error=logs/fastqc_%A_%a.err

module purge
module load fastqc/0.11.9

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

R1="${SAMPLE}_R1.fastq"
R2="${SAMPLE}_R2.fastq"

mkdir -p fastqc_out

fastqc \
  --threads ${SLURM_CPUS_PER_TASK} \
  -o fastqc_out \
  ${R1} ${R2}

