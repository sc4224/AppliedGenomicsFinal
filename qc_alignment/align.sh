#!/bin/bash
#SBATCH --job-name=star_align
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-1
#SBATCH --output=logs/star_%A_%a.out
#SBATCH --error=logs/star_%A_%a.err

module purge
module load star/intel/2.7.11a

GENOME_DIR="./reference_index"

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

R1="./trimmed_fastq/${SAMPLE}/${SAMPLE}_R1_paired.fastq"
R2="./trimmed_fastq/${SAMPLE}/${SAMPLE}_R2_paired.fastq"

OUTDIR="star_out/${SAMPLE}"
mkdir -p ${OUTDIR}

STAR \
  --runThreadN ${SLURM_CPUS_PER_TASK} \
  --genomeDir ${GENOME_DIR} \
  --readFilesIn ${R1} ${R2} \
  --outFileNamePrefix ${OUTDIR}/ \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts

