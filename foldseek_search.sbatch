#!/bin/bash
PDB=$1
name=$(basename ${PDB})
sbatch <<EOF
#!/bin/sh
#SBATCH --job-name=$name
#SBATCH --cpus-per-task=48
#SBATCH --mem=270g
#SBATCH --partition=norm
#SBATCH --time=12:00:00
#SBATCH --constraint=e7543
#SBATCH --output=${name}_slurm_2.out

ml foldseek


foldseek-10-941cd33 easy-search \
-e 9999 \
--format-output query,target,evalue,bits \
--threads 48 \
--max-seqs 10000 \
$PDB /fdb/foldseek/9-427df8a/Alphafold-UniProt50 foldseek_scop40_vs_afdb50_max_seq_10000_2.tsv tmp


EOF

