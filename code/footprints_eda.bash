#!/usr/bin/env bash

# Run cell-level exploratory data analysis notebooks (footprints_eda.ipynb) on all cell types
#
# Requires: ct-specific processed or pre-annotated anndatas

### SETUP ###

set -eou pipefail


# Args

USE_CLUSTER=false

while getopts "c" opt; do

	case "$opt" in

		c ) USE_CLUSTER=true ;;
		? ) echo "Usage: $0 [-c]"; exit 1 ;;
	
	esac

done

# Environment

PROJECT_DIR=".."
PROJECT_DIR="$(realpath ${PROJECT_DIR})"
cd $PROJECT_DIR

source code/glob_vars.bash # FOOTPRINTS_DIR, FOOTPRINTS_EDA, MAIN_ENV

mkdir -p "$FOOTPRINTS_EDA"


if [[ "$USE_CLUSTER" == "false" ]]; then

	source $HOME/.bash_profile
	load-micromamba
	micromamba activate $MAIN_ENV

fi


### SCRIPT ###

# Run function over ct-specific notebooks
while IFS= read -r -d '' adata_full; do

	# Make relevant variables
	
	adata="$(basename \"$adata_full\")"
	adata_short="${adata%_processed.h5ad}"

	cell_type="$(cut -d '_' -f 2 <<< \"${adata_short}\")"

	ipynb_out="${FOOTPRINTS_EDA}/eda_${cell_type}.ipynb"

	
	# Clean old run files and create dir
	if [[ -f "${ipynb_out}" ]]; then

		rm "${ipynb_out}" 

	fi

	if [[ -d "${FOOTPRINTS_EDA}/${cell_type}" ]]; then

		rm -rf "${FOOTPRINTS_EDA}/${cell_type}"

	fi

	mkdir -p "${FOOTPRINTS_EDA}/${cell_type}"


	# Create ct-specific notebook

	echo "Processing file: ${adata}"

	cat "code/footprints_eda.notebook.ipynb" |
	sed '/^    "    PROJECT_DIR = '\''manual'\''/c\    "    PROJECT_DIR = '\'"${PROJECT_DIR}"\''\\n",' |
	sed '/^    "adata_path = /c\    "adata_path = '\'"${adata_full}"\''\\n",' |
	sed '/^    "cell_type = /c\    "cell_type = '\'"$cell_type"\''\\n",' > "$ipynb_out"


	case "$USE_CLUSTER" in

		false )

			jupyter nbconvert --to notebook --execute --allow-error --inplace "$ipynb_out"
			;;

		true )

			job_id="footprints_eda_$(date '+%Y-%m-%d')_${adata_short}_BORGS"
			bsub <<EOF
#!/usr/bin/env bash
#BSUB -R "rusage[mem=10G]"
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

jupyter nbconvert --to notebook --execute --allow-errors --inplace "$ipynb_out"
EOF

			;;

	esac

done < <(find "${FOOTPRINTS_DIR}" -mindepth 1 -maxdepth 1 \( -type f -o -type l \) -iname "*_processed.h5ad" -print0)
