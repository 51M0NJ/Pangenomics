#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=31G

mamba activate Coinfinder
python /nfs3_ib/nfs-ip34/home/def-jacquesp/jeanneau/9-Pangenomics/2_Produce_cooccurrence_matrix/1_chunk_files.py
