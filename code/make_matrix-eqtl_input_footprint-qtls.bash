#!/usr/bin/env bash

# Run pipeline: make_matrix-eqtl_input_footprint-qtls.ipynb on all cell types.
#
# Requires:
# 	- ct-specific processed and pre-annotated footprint anndatas
# 	- genotype tsv and pcs

### Setup ###

set -eou pipefail

## Variables

MODE="unset"
USE_CLUSTER=false

while getopts "m:c" opt; do

	case "$opt" in

		m ) MODE=$OPTARG ;;
		c ) USE_CLUSTER=true ;;
		? ) echo "Usage: $0 [-c] cluster [-m mode {single-test, peak-tests, bulk-tests}"; exit 1;;

	esac

done

if [[ "$MODE" == "unset" ]]; then

	echo "Usage: $0 [-c] cluster [-m mode {single-tests, peak-tests, bulk-tests}"
	exit 1

fi

DATE="$(date +"%Y-%m-%d")"


## Environment

PROJECT_DIR=".."
PROJECT_DIR="$(realpath ${PROJECT_DIR})"
cd $PROJECT_DIR

source code/glob_vars.bash # FOOTPRINTS_DIR, MATRIX_EQTL_INPUT_DIR, MAIN_ENV, CT_MAP_ID

if [[ "$USE_CLUSTER" == "false" ]]; then

	source $HOME/.bash_profile
	load-micromamba
	micromamba activate main06

fi


### SCRIPT ###

while IFS= read -r -d '' algorithm_path; do

	algorithm="$(basename "$algorithm_path")"

	echo "$algorithm"


	while IFS= read -r -d '' peak_set_path; do

		peak_set="$(basename "$peak_set_path")"

		if [[ "$peak_set" == *old* ]]; then

			continue

		fi

		echo -e "\t$peak_set"


		cmd_list=()

		while IFS= read -r -d '' cell_type_path; do

			cell_type="$(basename "$cell_type_path")"

			echo -e "\t\t$cell_type"

			
			# Create ipynb notebook & customize
			
			ipynb_out="${MATRIX_EQTL_INPUT_DIR}/footprints/${algorithm}/${peak_set}/${CT_MAP_ID}/${cell_type}/${MODE}/make_matrix-eqtl_input_footprint-qtls.ipynb"
			category="${algorithm}/${peak_set}/${CT_MAP_ID}/${cell_type}"

			mkdir -p "$(dirname "$ipynb_out")"

			cat "code/make_matrix-eqtl_input_footprint-qtls.ipynb" |
				sed '/^    "    PROJECT_DIR = '\''manual'\''/c\    "    PROJECT_DIR = '\'"${PROJECT_DIR}"\''\\n",' |
				sed '/^    "category = str(/c\    "category = '\'"$category"\''\\n",' |
				sed '/^    "cell_type = str(/c\    "cell_type = '\'"$cell_type"\''\\n",' |
				sed '/^    "mode = /c\    "mode = '\'"$MODE"\''\\n",' > "$ipynb_out"


			## Run

			case "$USE_CLUSTER" in

				false )
					
					# Remove old files
					find "$(dirname "$ipynb_out")" -mindepth 1 -depth ! -name "$(basename "$ipynb_out")" -exec rm -rf {} \;

					# Execute notebook
					jupyter nbconvert --to notebook --execute --allow-error --inplace "$ipynb_out"

					;;

				true )

					job_id="meqtl_i_fp_${DATE}_${algorithm:0:4}_${peak_set}_${cell_type}"
					bsub <<EOF
#!/usr/bin/env bash
#BSUB -R "rusage[mem=100G]"
#BSUB -q long
#BSUB -cwd ${PROJECT_DIR}
#BSUB -J ${job_id}
#BSUB -o ${PROJECT_DIR}/code/bsub/logs/${job_id}.out
#BSUB -e ${PROJECT_DIR}/code/bsub/logs/${job_id}.err

set -euo pipefail

cd $PROJECT_DIR

source "$HOME/.bash_profile"
load-micromamba
micromamba activate $MAIN_ENV
			
# Remove old files
find "$(dirname "$ipynb_out")" -mindepth 1 -depth ! -name "$(basename "$ipynb_out")" -exec rm -rf {} \\;

# Execute notebook
jupyter nbconvert --to notebook --execute --allow-errors --inplace "$ipynb_out"
EOF

					;;

			esac

		done < <(find "${peak_set_path}/${CT_MAP_ID}" -mindepth 1 -maxdepth 1 -type d -print0)

	done < <(find "${algorithm_path}" -mindepth 1 -maxdepth 1 -type d -print0)

done < <(find "${FOOTPRINTS_DIR}" -mindepth 1 -maxdepth 1 -type d -print0)
