#!/usr/bin/env bash

# Select peaks for caQTL testing
#
# Requires:
#

### Setup ###

set -eou pipefail

## Variables

USE_CLUSTER=false

while getopts "c" opt; do

	case "$opt" in

		c ) USE_CLUSTER=true ;;
		? ) echo "Usage: $0 [-c] cluster"; exit 1;;

	esac

done

## Environment

PROJECT_DIR=".."
PROJECT_DIR="$(realpath ${PROJECT_DIR})"
cd $PROJECT_DIR

source code/glob_vars.bash # MAIN_ENV, CT_MAP_ID, SELECT_PEAKS_TSV_DIR, ATAC_CHROM_ACCESS_DIR

if [[ "$USE_CLUSTER" == "false" ]]; then

	source $HOME/.bash_profile
	load-micromamba
	micromamba activate ${MAIN_ENV}

fi


### SCRIPT ###

while IFS= read -r -d '' cell_type_path; do

	cell_type="$(basename "$cell_type_path")"

	echo -e "Processing: $cell_type"


	# Cell-type specific thresholds, manuaally set a-posteriori
	

	case "$cell_type" in


		DL-EN )
			min_mean_acells="0.0125"
			max_mean_acells="3"
			min_donors="0.3"
			min_mean_adonors="0.02"
			min_score="9"
			;;

		Glia )
			min_mean_acells="0.0125"
			max_mean_acells="3"
			min_donors="0.4"
			min_mean_adonors="0.02"
			min_score="9"
			;;

		Midbrain-EN )
			min_mean_acells="0.0125"
			max_mean_acells="3"
			min_donors="0.25"
			min_mean_adonors="0.02"
			min_score="9"
			;;

		Neural-progenitors )
			min_mean_acells="0.0125"
			max_mean_acells="4.5"
			min_donors="0.35"
			min_mean_adonors="0.02"
			min_score="9"
			;;

		UL-EN )
			min_mean_acells="0.0125"
			max_mean_acells="3"
			min_donors="0.35"
			min_mean_adonors="0.02"
			min_score="9"
			;;

		* )
			echo "Cell-type not found! Aborting."
			exit 1
			;;

	esac

	
	# Create ipynb notebook & customize
	
	ipynb_in="code/select_peaks_and_create_ca-qtl_adata.ipynb"
	ipynb_out="${SELECT_PEAKS_TSV_DIR}/${cell_type}/select_peaks_and_create_ca-qtl_adata.ipynb"
	mkdir -p "$(dirname "$ipynb_out")"

	cat "$ipynb_in" |
		sed '/^    "cell_type = str(/c\    "cell_type = '\'"$cell_type"\''"' |
		sed '/^    "min_mean_acells = /c\    "min_mean_acells = '"${min_mean_acells}"'\\n",' |
		sed '/^    "max_mean_acells = /c\    "max_mean_acells = '"${max_mean_acells}"'\\n",' |
		sed '/^    "min_donors = /c\    "min_donors = '"${min_donors}"'\\n",' |
		sed '/^    "min_mean_adonors = /c\    "min_mean_adonors = '"${min_mean_adonors}"'\\n",' |
		sed '/^    "min_score = /c\    "min_score = '"${min_score}"'\\n"' |
	    sed '/^    "    PROJECT_DIR = '\''manual'\''/c\    "    PROJECT_DIR = '\'"${PROJECT_DIR}"\''\\n",' > "$ipynb_out"


	## Run

	case "$USE_CLUSTER" in

		false )
			
			# Remove old files
			find "$(dirname "$ipynb_out")" -mindepth 1 -depth ! -name "$(basename "$ipynb_out")" -exec rm -rf {} \;

			# Execute notebook
			jupyter nbconvert --to notebook --execute --allow-error --inplace "$ipynb_out"

			;;

		true )

			job_id="select_peaks_$(date +"%Y-%m-%d")_${cell_type}"
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

done < <(find "${ATAC_CHROM_ACCESS_DIR}/adata/subset/${CT_MAP_ID}" -mindepth 1 -maxdepth 1 -type d -print0)
