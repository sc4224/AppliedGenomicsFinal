#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=1-1
#SBATCH --output=logs/trimmomatic_%A_%a.out
#SBATCH --error=logs/trimmomatic_%A_%a.err

module purge
module load trimmomatic/0.39

TRIMMOMATIC_JAR=/share/apps/trimmomatic/0.39/trimmomatic-0.39.jar
ADAPTERS=/share/apps/trimmomatic/0.39/adapters/TruSeq3-PE.fa

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

R1="${SAMPLE}_R1.fastq"
R2="${SAMPLE}_R2.fastq"

OUTDIR="trimmed_fastq/${SAMPLE}"
mkdir -p "${OUTDIR}"

java -jar "${TRIMMOMATIC_JAR}" PE \
  -threads ${SLURM_CPUS_PER_TASK} \
  -phred33 \
  "${R1}" "${R2}" \
  "${OUTDIR}/${SAMPLE}_R1_paired.fastq" \
  "${OUTDIR}/${SAMPLE}_R1_unpaired.fastq" \
  "${OUTDIR}/${SAMPLE}_R2_paired.fastq" \
  "${OUTDIR}/${SAMPLE}_R2_unpaired.fastq" \
  ILLUMINACLIP:"${ADAPTERS}":2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

