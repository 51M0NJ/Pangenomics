#!/bin/bash
#SBATCH --job-name=geneproc_%a
#SBATCH --output=logs/geneproc_%a.out
#SBATCH --error=logs/geneproc_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=31G
#SBATCH --array=0-486

mamba activate Coinfinder
mkdir -p /net/nfs-ip34/home/def-jacquesp/jeanneau/9-Pangenomics/2_Produce_cooccurrence_matrix/logs
python /net/nfs-ip34/home/def-jacquesp/jeanneau/9-Pangenomics/2_Produce_cooccurrence_matrix/3_process_chunk.py $SLURM_ARRAY_TASK_ID
