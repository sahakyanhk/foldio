#!/bin/bash
# Submit an ESMFold2 prediction job.
#   ./esmfold2.sbatch <input> <output_dir> [walltime]
#
# <input>:
#   * a FASTA file (.fasta/.fa/.fas/.faa) ... each record folded as a monomer
#   * a JSON file (.json) ................... queries (monomers/complexes)
#   * anything else is treated as a LIST file: one MSA path per line; the first
#     sequence of each MSA is folded using that MSA (aligned FASTA / a3m / sto).
# Everything runs in a SINGLE job: esmfold2.py loads the model once, streams the
# queries, and skips already-done outputs, so a timed-out job just needs resubmit.
set -euo pipefail

input="$(realpath "$1")"
output="$(realpath -m "$2")"
walltime="${3:-0-01:00:00}"

scriptdir="$(cd "$(dirname "$0")" && pwd)"
ext="${input##*.}"

case "${ext,,}" in
    fasta|fa|fas|faa)
        fold_input="-ifasta $input" ;;
    json)
        fold_input="-json $input" ;;
    a3m|sto|stockholm|sth|stk|aln|afa|a2m)
        echo "ERROR: '$input' looks like a single MSA. For MSA prediction pass a" >&2
        echo "       LIST file (one MSA path per line), e.g.  ls *.a3m > list.txt" >&2
        exit 1 ;;
    *)
        fold_input="-inputmsa $input" ;;   # list of MSA paths, one per line
esac

mkdir -p "$output"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=$(basename "$input")
#SBATCH --cpus-per-task=2
#SBATCH --mem=16g
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --time=$walltime
#SBATCH --output=$output/esmfold2.out

set -euo pipefail
source ~/.bashrc
micromamba activate esmfold2

python "$scriptdir/esmfold2.py" $fold_input -output "$output"
EOF
