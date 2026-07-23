#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")"

mkdir -p ../reseek_search

run_reseek() {

query=$1
db=$2

qname=$(basename ${query%.*})
dbname=$(basename ${db%.*})

outname="../reseek_search/${qname}_${dbname}"


sbatch <<EOF
#!/bin/sh
#SBATCH --job-name=$dbname
#SBATCH --cpus-per-task=96
#SBATCH --mem=40g
#SBATCH --partition=norm
#SBATCH --time=48:00:00
#SBATCH --constraint=e9454
#SBATCH --output=${outname}.slurm

/bin/time -v -o ${outname}.time reseek-v2.5  \
	  -search $query \
	  -db $db \
	  -output ${outname}.tsv \
	  -fast \
	  -dbsize 10000 \
	  -evalue 9999 \
	  -threads 96



EOF

}


# run_reseek ../xdata/scop40.bca ../xdata/scop40.bca
# run_reseek ../xdata/scop40.bca ../xdata/X10_scop40.bca
# run_reseek ../xdata/scop40.bca ../xdata/X100_scop40.bca
# run_reseek ../xdata/scop40.bca ../xdata/X1000_scop40.bca
# run_reseek ../xdata/scop40.bca ../xdata/X5000_scop40.bca

for n in 1 10 100 1000 5000
do
	run_reseek ../xdata/scop40a.bca ../xdata/X${n}_scop40a.bca	# TPs
	run_reseek ../xdata/scop40b.bca ../xdata/X${n}_scop40a.bca	# FPs
done
