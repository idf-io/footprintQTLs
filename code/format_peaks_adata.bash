#!/usr/bin/env bash

# Cluster execute: format_anndata.ipynb

set -euo pipefail


PROJECT_DIR=".."
PROJECT_DIR="$(realpath "$PROJECT_DIR")"
cd $PROJECT_DIR

source "code/glob_vars.bash" # MAIN_ENV

JOB_ID="format_anndata_$(date +"%Y-%m-%d")"

bsub << EOF
#!/usr/bin/env bash
#BSUB -R "rusage[mem=100G]"
#BSUB -q medium
#BSUB -J $JOB_ID
#BSUB -cwd $PROJECT_DIR

#BSUB -o ${PROJECT_DIR}/code/bsub/logs/${JOB_ID}_%I.out
#BSUB -e ${PROJECT_DIR}/code/bsub/logs/${JOB_ID}_%I.err

set -euo pipefail

cd $PROJECT_DIR

head -n50 "${HOME}/.bash_profile" | tail -n10

source "${HOME}/.bash_profile"
echo $(pwd -P)
echo $PROJECT_DIR
load-micromamba
micromamba activate $MAIN_ENV

jupyter nbconvert --to notebook --execute --allow-errors --inplace code/format_peaks_adata.ipynb
EOF
