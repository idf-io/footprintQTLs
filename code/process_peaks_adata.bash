#!/usr/bin/env bash

# Cluster execute: process_anndata.ipynb

set -euo pipefail


PROJECT_DIR=".."
PROJECT_DIR="$(realpath "$PROJECT_DIR")"
cd $PROJECT_DIR

source "code/glob_vars.bash" # MAIN_ENV

JOB_ID="process_peaks-adata_$(date +"%Y-%m-%d")"

bsub << EOF
#!/usr/bin/env bash
#BSUB -R "rusage[mem=400G]"
#BSUB -q highmem
#BSUB -J $JOB_ID
#BSUB -cwd $PROJECT_DIR

#BSUB -o ${PROJECT_DIR}/code/bsub/logs/${JOB_ID}_%I.out
#BSUB -e ${PROJECT_DIR}/code/bsub/logs/${JOB_ID}_%I.err

set -euo pipefail

cd $PROJECT_DIR

source "${HOME}/.bash_profile"
load-micromamba
micromamba activate $MAIN_ENV

jupyter nbconvert --to notebook --execute --allow-errors --inplace code/process_peaks_adata.ipynb
EOF
