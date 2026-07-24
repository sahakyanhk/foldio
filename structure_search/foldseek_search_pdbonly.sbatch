#!/bin/bash

#foldseek_search.sbatch <pdbdir/pdbfile> <output_file> 

PDB=$(realpath "$1")
output=$(realpath "$2") 
name=$(basename "$PDB")

sbatch <<EOF
#!/bin/sh
#SBATCH --job-name=$name
#SBATCH --cpus-per-task=24
#SBATCH --mem=30g
#SBATCH --partition=norm
#SBATCH --time=00:15:00
#SBATCH --gres=lscratch:30
#SBATCH --output=${name}_slurm.out


$(dirname $(realpath $0))/../bin/foldseek easy-search \
-e 9999 \
--format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits \
--threads 96 \
--max-seqs 1000 \
$PDB /fdb/foldseek/9-427df8a/PDB ${output%.*}.tsv /lscratch/\$SLURM_JOB_ID

EOF

