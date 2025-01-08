#!/usr/bin/env bash

# Execute process_anndata.ipynb

### Setup ##

# Env

set -euo pipefail

PROJECT_DIR=".."
PROJECT_DIR="$(realpath "$PROJECT_DIR")"
cd $PROJECT_DIR

# Vars

source "code/glob_vars.bash" # MAIN_ENV, ATAC_CHROM_ACCESS_DIR

USE_CLUSTER="false"

while getopts "c" opt; do

	case "$opt" in

		c )
			USE_CLUSTER="true"
			;;

		* )
			echo "Wrong flag"
			;;

	esac

done


### Script ###

ipynb_file_template="code/process_peaks_adata.ipynb"
ipynb_file_exec="${ATAC_CHROM_ACCESS_DIR}/adata/process_peaks_adata.ipynb"

mkdir -p "$(dirname "$ipynb_file_exec")"

cat "$ipynb_file_template" | \
	sed '/^    "    PROJECT_DIR = '\''manual'\''/c\    "    PROJECT_DIR = '\'"${PROJECT_DIR}"\''\\n",' > "$ipynb_file_exec"


case "$USE_CLUSTER" in

	false )

		source "${HOME}/.bash_profile"
		load-micromamba
		micromamba activate $MAIN_ENV

		jupyter nbconvert --to notebook --execute --allow-errors --inplace "$ipynb_file_exec"
		;;

	true )

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

jupyter nbconvert --to notebook --execute --allow-errors --inplace "$ipynb_file_exec"
EOF
		;;

esac
