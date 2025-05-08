#!/bin/bash
#SBATCH --job-name=star_build
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-1
#SBATCH --output=logs/star_build_%A_%a.out
#SBATCH --error=logs/star_build_%A_%a.err

module purge
module load star/intel/2.7.11a

STAR \
  --runMode genomeGenerate \
  --runThreadN 8 \
  --genomeDir ./reference_index \
  --genomeFastaFiles GCF_000001635.27_GRCm39_genomic.fna \
  --sjdbGTFfile      genomic.gtf \
  --sjdbOverhang     99

