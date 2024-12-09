#!/usr/bin/env bash

# Run pipeline: make_matrix-eqtl_input_footprint-qtls.ipynb on all cell types.
#
# Requires:
# 	- ct-specific processed and pre-annotated footprint anndatas
# 	- genotype vcf, tsv and pcs

set -eou pipefail

# Setup variables and env

USE_CLUSTER=false

while getopts "c" opt; do

	case "$opt" in

		c ) USE_CLUSTER=true ;;
		? ) echo "Usage: $0 [-c]"; exit 1 ;;
	
	esac

done

PROJECT_DIR=".."
PROJECT_DIR="$(realpath ${PROJECT_DIR})"
cd $PROJECT_DIR

source code/glob_vars.bash # FOOTPRINTS_DIR, MATRIX_EQTL_INPUT_DIR


if [[ "$USE_CLUSTER" == "false" ]]; then

	source $HOME/.bash_profile
	load-micromamba
	micromamba activate main04

fi

DATE="$(date +"%Y-%m-%d")"




# Run function over ct-specific notebooks

while IFS= read -r -d '' adatas_long; do

	# Make relevant variables
	adata="$(basename $adatas_long)"
	adata_short="${adata%_processed.h5ad}"
	cell_type="$(echo "$adata_short" | cut -d '_' -f 2)"
	ipynb_out="${MATRIX_EQTL_INPUT_DIR}/${cell_type}/make_matrix-eqtl_input_footprint-qtls.ipynb"


	# Remove old files and create dir

	if [[ -d "${MATRIX_EQTL_INPUT_DIR}/${cell_type}" ]]; then

		rm -rf "${MATRIX_EQTL_INPUT_DIR}/${cell_type}"

	fi

	mkdir -p "$(dirname $ipynb_out)"


	# Create ct-specific notebook

	echo "Processing file: ${ipynb_out}"

	cat "code/make_matrix-eqtl_input_footprint-qtls.ipynb" |
	sed '/^    "    PROJECT_DIR = '\''manual'\''/c\    "    PROJECT_DIR = '\'"${PROJECT_DIR}"\''\\n",' |
	sed '/^    "cell_type = str(/c\    "cell_type = '\'"$cell_type"\''\\n",' |
	sed '/^    "RUN_ID = /c\    "RUN_ID = '\'"meqtl_io_${DATE}_${CT_MAP_ID}_${DATASET}"\''\\n",' > "$ipynb_out"


	case "$USE_CLUSTER" in

		false )

			jupyter nbconvert --to notebook --execute --allow-error --inplace "$ipynb_out"
			;;

		true )

			job_id="meqtl_io_${DATE}_${cell_type}_BORGS"
			bsub <<EOF
#!/usr/bin/env bash
#BSUB -R "rusage[mem=300G]"
#BSUB -q highmem
#BSUB -cwd ${PROJECT_DIR}
#BSUB -J ${job_id}
#BSUB -o ${PROJECT_DIR}/code/bsub/logs/${job_id}.out
#BSUB -e ${PROJECT_DIR}/code/bsub/logs/${job_id}.err

set -euo pipefail

cd $PROJECT_DIR

source "$HOME/.bash_profile"
load-micromamba
micromamba activate main04

jupyter nbconvert --to notebook --execute --allow-errors --inplace "$ipynb_out"
EOF

			;;

	esac

done < <(find "${FOOTPRINTS_DIR}" -mindepth 1 -maxdepth 1 \( -type f -o -type l \) -iname "*_processed.h5ad" -print0)
