#!/bin/bash

#foldseek_search.sbatch <pdbdir/pdbfile> <output_file> 

PDB=$(realpath "$1")
output=$(realpath "$2") 
name=$(basename "$PDB")

sbatch <<EOF
#!/bin/sh
#SBATCH --job-name=$name
#SBATCH --cpus-per-task=96
#SBATCH --mem=350g
#SBATCH --partition=largemem
#SBATCH --time=8:00:00
#SBATCH --gres=lscratch:50
#SBATCH --output=${name}_slurm.out


$(dirname $(realpath $0))/../bin/foldseek easy-search \
-e 9999 \
--format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits \
--threads 96 \
--max-seqs 1000 \
$PDB /fdb/foldseek/9-427df8a/PDB ${output%.*}_pdbonly.tsv /lscratch/\$SLURM_JOB_ID



$(dirname $(realpath $0))/../bin/foldseek easy-search \
-e 9999 \
--format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits \
--threads 96 \
--max-seqs 1000 \
$PDB /fdb/foldseek/9-427df8a/Alphafold-UniProt50 ${output%.*}_afdb50.tsv /lscratch/\$SLURM_JOB_ID

cat ${output%.*}_pdbonly.tsv ${output%.*}_afdb50.tsv > $output

rm ${output%.*}_pdbonly.tsv ${output%.*}_afdb50.tsv


EOF

