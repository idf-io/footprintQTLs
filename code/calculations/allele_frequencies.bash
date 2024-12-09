#!/usr/bin/env bash
#
# Calculate genotype allele frequencies.
#
# Input:
#     - GENOTYPES.TSV 
#
# Output:
#     - GENOTYPE_DIR/allele_frequencies.tsv
#     - GENOTYPE_DIR/allele_frequencies_stats.tsv
#
#
#BSUB -R "rusage[mem=200G]"
#BSUB -q long
#BSUB -cwd "/omics/groups/OE0540/internal_temp/users/fichtner/projects/footprintQTL/code/calculations"
#BSUB -J "allele_frequency"
#BSUB -o "/omics/groups/OE0540/internal_temp/users/fichtner/projects/footprintQTL/code/bsub/logs/%I.out"
#BSUB -e "/omics/groups/OE0540/internal_temp/users/fichtner/projects/footprintQTL/code/bsub/logs/%I.err"

echo "This script is very computation intensive. Recommended to submit as cluster job. Continuing with the script..."


### Setup ###

set -euo pipefail

pwd
source "${HOME}/.bash_profile"
source "../glob_vars.bash" # MAIN_ENV
load-micromamba
micromamba activate "$MAIN_ENV"


python <<EOF
import os
import sys
import pandas as pd

os.chdir('../..')

sys.path.append('code')
from glob_vars import GENOTYPES_PROCESSED_TSV, GENOTYPE_DIR
from helpers.python.utils import create_dir


gt = pd.read_csv(GENOTYPES_PROCESSED_TSV, sep='\t', header=0, index_col=0)


# Allele frequencies across donors

allele_freqs = gt.T.apply(pd.value_counts, normalize=True).T
allele_freqs = allele_freqs.round(3)

allele_freqs_out = f'{GENOTYPE_DIR}/metadata/allele_frequencies.tsv'
create_dir(allele_freqs_out)
allele_freqs.to_csv(allele_freqs_out, sep='\t', na_rep='NaN')
print(f'Saved allele frequencies at: {allele_freqs_out}')


# Summary of allele frequencies

allele_freqs_stats = allele_freqs.describe()
allele_freqs_stats = allele_freqs_stats.round(3)

allele_freqs_stats_out = f'{GENOTYPE_DIR}/metadata/allele_frequencies_stats.tsv'
create_dir(allele_freqs_stats_out)
allele_freqs_stats.to_csv(allele_freqs_stats_out, sep='\t')
print(f'Saved allele frequency stats at: {allele_freqs_stats_out}')
EOF
