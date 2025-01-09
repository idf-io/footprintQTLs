#!/usr/bin/env bash

# Make Matrix eQTL input for caQTLs
#
# Requires:

### Setup ###

set -eou pipefail

## Variables

MODE="unset"
MODE="bulk-tests"
USE_CLUSTER=false

while getopts "m:c" opt; do

	case "$opt" in

		#m ) MODE=$OPTARG ;;
		c ) USE_CLUSTER=true ;;
		? ) echo "Usage: $0 [-c] cluster"; exit 1;;

	esac

done

if [[ "$MODE" == "unset" ]]; then

	echo "Usage: $0 [-c] cluster [-m mode {single-tests, peak-tests, bulk-tests}"
	exit 1

fi


## Environment

PROJECT_DIR=".."
PROJECT_DIR="$(realpath ${PROJECT_DIR})"
cd $PROJECT_DIR

source code/glob_vars.bash # ATAC_CHROM_ACCESS_DIR #FOOTPRINTS_DIR, MATRIX_EQTL_INPUT_DIR, MAIN_ENV, CT_MAP_ID

if [[ "$USE_CLUSTER" == "false" ]]; then

	source $HOME/.bash_profile
	load-micromamba
	micromamba activate $MAIN_ENV

fi


### SCRIPT ###

while IFS= read -r -d '' cell_type_path; do

	cell_type="$(basename "$cell_type_path")"

	echo -e "$cell_type"

	
	# Create ipynb notebook & customize
	
	ipynb_out="${MATRIX_EQTL_INPUT_DIR}/chromatin-accessibility//${CT_MAP_ID}/${cell_type}/${MODE}/make_matrix-eqtl_input_chromatin-accessibility-qtls.ipynb"
	ipynb_in="code/make_matrix-eqtl_input_chromatin-accessibility-qtls.ipynb"
	mkdir -p "$(dirname "$ipynb_out")"

	cat "$ipynb_in" |
		sed '/^    "    PROJECT_DIR = '\''manual'\''/c\    "    PROJECT_DIR = '\'"${PROJECT_DIR}"\''\\n",' |
		sed '/^    "cell_type = str(/c\    "cell_type = '\'"$cell_type"\''\\n",' |
		sed '/^    "mode = /c\    "mode = '\'"$MODE"\''\\n",' > "$ipynb_out" > "$ipynb_out"


	## Run

	case "$USE_CLUSTER" in

		false )
			
			# Remove old files
			find "$(dirname "$ipynb_out")" -mindepth 1 -depth ! -name "$(basename "$ipynb_out")" ! -name "phenotype.tsv" -exec rm -rf {} \;

			# Execute notebook
			jupyter nbconvert --to notebook --execute --allow-error --inplace "$ipynb_out"

			;;

		true )

			job_id="meqtl_i_ca_$(date +"%Y-%m-%d")_${cell_type}"
			bsub <<EOF
#!/usr/bin/env bash
#BSUB -R "rusage[mem=100G]"
#BSUB -q medium
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

done < <(find "${ATAC_CHROM_ACCESS_DIR}/adata/subset/${CT_MAP_ID}" -mindepth 1 -maxdepth 1 -type d -print0)
