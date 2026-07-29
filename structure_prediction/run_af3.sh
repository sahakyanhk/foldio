#!/bin/bash
set -eou pipefail

# #predict a single protein in a single fasta. This takes a fasta with a single sequence, and the predicted results will be in the <output_dir>
# ./run_af3.sh <input.fasta> <output_dir>

# #multiple fasta files with a protein in each. If you put two proteins or protein/RNA it will predict complexes
# mkdir -p AF3runs # all outputs will go here
# for fasta in <your_dir_with_fastas/*>; do
# 	fasta_name=$(basename ${fasta%.*})
# 	./run_af3.sh $fasta AF3runs/$fasta_name
# done

# #once all jobs are finished, collect all *cif files and extract plddts
# mkdir -p cif_files; cp AF3runs/*/*/*_model.cif cif_files # make a dir "cif_files" and collect all predicted models there
# #run this loop to extract the mean plddt from .cif files and save in "plddt_summary.txt" 
# for cif in cif_files/*cif; do
# 	echo $(basename $cif .cif) \
#  		 $(awk '$1=="ATOM" { sum += $15; count++ } END { if (count > 0) print sum / count }' $cif);
# done > plddt_summary.txt


if [ $# -ne 2 ]; then
    echo "usage: $(basename "$0") <input.fasta> <output_dir>" >&2
    exit 1
fi

# resolve everything to absolute paths so the script (and the sbatch jobs it
# submits) work no matter what directory it is called from
inputFASTA=$(readlink -f "$1")
outdir=$2

if [ ! -f "$inputFASTA" ]; then
    echo "error: input fasta not found: $1" >&2
    exit 1
fi

mkdir -p "$outdir"
outdir=$(readlink -f "$outdir")

name=$(basename "${inputFASTA%.*}")
json1="${outdir}/${name}.json"
json2="${outdir}/${name}/${name}_data.json"


modeldir="/data/$USER/.af3_model"
modelpath="${modeldir}/af3.bin.zst"

if [ ! -f "$modelpath" ]; then
    mkdir -p "$modeldir"
    curl -o "$modelpath" https://storage.googleapis.com/alphafold3/af3.bin.zst
fi

# Step 1: Make MSA
jobid1=$(sbatch --parsable <<EOF
#!/bin/bash
#SBATCH --job-name=${name}_prep
#SBATCH --ntasks=16
#SBATCH --mem=32g
#SBATCH --partition=norm
#SBATCH --time=8:00:00
#SBATCH --output=${outdir}/af3_1.out

ml alphafold3

af3 convert --guess --output-dir "$outdir" --seeds 1,2,3,4,5 "$inputFASTA"
af3 run --model_dir "$modeldir" --output_dir="$outdir" --json_path "$json1" --norun_inference
EOF
)

# Step 2: Predict structure
sbatch --dependency=afterok:$jobid1 <<EOF
#!/bin/bash
#SBATCH --job-name=${name}_infer
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g
#SBATCH --partition=gpu
#SBATCH --gres=lscratch:150,gpu:a100:1
#SBATCH --time=2:00:00
#SBATCH --output=${outdir}/af3_2.out

ml alphafold3

af3 run --norun_data_pipeline --json_path "$json2" --model_dir "$modeldir" --output_dir="$outdir"

EOF

