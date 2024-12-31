#!/usr/bin/env bash

# Run exploratory data analysis notebooks on all cell types
#
# Input: ct-specific processed or pre-annotated anndatas
# Output:

### SETUP ###

set -eou pipefail


## Args

USE_CLUSTER=false

while getopts "c" opt; do

	case "$opt" in

		c ) USE_CLUSTER=true ;;
		? ) echo "Usage: $0 [-c]"; exit 1 ;;
	
	esac

done

## Environment

PROJECT_DIR=".."
PROJECT_DIR="$(realpath ${PROJECT_DIR})"
cd $PROJECT_DIR

source code/glob_vars.bash # FOOTPRINTS_DIR, FOOTPRINTS_EDA, MAIN_ENV, CT_MAP_ID

if [[ "$USE_CLUSTER" == "false" ]]; then

	source $HOME/.bash_profile
	load-micromamba
	micromamba activate $MAIN_ENV

fi


### SCRIPT ###

while IFS= read -r -d '' algorithm_path; do

	algorithm="$(basename "$algorithm_path")"
	if [[ "$algorithm" != "counts" ]]; then
		continue
	fi

	echo "$algorithm"


	while IFS= read -r -d '' peak_set_path; do

		peak_set="$(basename "$peak_set_path")"

		if [[ "$peak_set" == *old* ]]; then

			continue

		fi

		if [[ "$peak_set" != "ca-qtls_variant-centred_15bp" ]]; then
			continue
		fi
		echo -en "\t$peak_set\n"


		cell_type_dir="${peak_set_path}/${CT_MAP_ID}"
		ipynb_out="${cell_type_dir}/eda_all-cell-types.ipynb"


		## Remove old run file and create dir

		if [[ -f "$ipynb_out" ]]; then

			rm "$ipynb_out"

		fi

		mkdir -p "$(dirname "$ipynb_out")"


		## Make new ipynb file and modify variables
		
		cat "code/footprints_eda_all-cts.notebook.ipynb" |
			sed '/^    "    PROJECT_DIR = '\''manual'\''/c\    "    PROJECT_DIR = '\'"${PROJECT_DIR}"\''\\n",' |
			sed '/^    "cell_type_dir = str(/c\    "cell_type_dir = '\'"$cell_type_dir"\''"' > "$ipynb_out"


		case "$USE_CLUSTER" in

			false )

				jupyter nbconvert --to notebook --execute --allow-error --inplace "$ipynb_out"
				;;

			true )

				job_id="fp_eda_all-cts_$(date '+%Y-%m-%d')_${algorithm:0:4}_${peak_set}"
				bsub <<EOF
#!/usr/bin/env bash
#BSUB -R "rusage[mem=10G]"
#BSUB -q short
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

	done < <(find "${algorithm_path}" -mindepth 1 -maxdepth 1 -type d -print0)

done < <(find "${FOOTPRINTS_DIR}" -mindepth 1 -maxdepth 1 -type d -print0)
