#!/usr/bin/env bash

# Process ct-specific footprint anndatas by calculating the descriptive stats on var-level, 
# highly variable peaks, PCA, UMAP and clustering.
#
# Input:
# 	- footprint anndatas 
# 		- adata.obs['n_cells'] annotation
# Output:
# 	- footprint anndatas processed (see above)

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

source code/glob_vars.bash # FOOTPRINTS_DIR, CT_MAP_ID

if [[ "$USE_CLUSTER" == "false" ]]; then

	source $HOME/.bash_profile
	load-micromamba
	micromamba activate main05

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

		echo -en "\t$peak_set"

		case "$peak_set" in

			ca-qtls )
				min_peak_counts=0
				;;

			ca-qtls_variant-centred_15bp )
				min_peak_counts=10
				;;

			ca-qtls_variant-centred_25bp )
				min_peak_counts=13
				;;

			ca-qtls_variant-centred_51bp )
				min_peak_counts=20
				;;

			ca-qtls_variant-centred_101bp )
				min_peak_counts=40
				;;

			ca-qtls_pm1k_variant-centred_15bp )
				min_peak_counts=10
				;;

			ca-qtls_pm1k_variant-centred_25bp )
				min_peak_counts=13
				;;

			ca-qtls_pm1k_variant-centred_51bp )
				min_peak_counts=20
				;;

			ca-qtls_pm1k_variant-centred_101bp )
				min_peak_counts=40
				;;

			* )
				echo -e "\tNo min_peak_count hard coded -> set to 0"
				min_peak_counts=0
				;;

		esac

		echo -ne " -> min counts per peak = ${min_peak_counts}\n"


		cmd_list=()

		while IFS= read -r -d '' cell_type_path; do

			cell_type="$(basename "$cell_type_path")"

			echo -e "\t\t$cell_type"


			cmd_main="python code/helpers/python/process_footprint_adata.py \
						-a '${FOOTPRINTS_DIR}/${algorithm}/${peak_set}/${CT_MAP_ID}/${cell_type}/footprints_pre-annotated.h5ad' \
						-o '${FOOTPRINTS_DIR}/${algorithm}/${peak_set}/${CT_MAP_ID}/${cell_type}/footprints_processed.h5ad' \
						-m '200' \
						-p '${min_peak_counts}'"

			cmd_main="$(echo $cmd_main | tr '\t' ' ')"


			case "$USE_CLUSTER" in

				false )

					eval "$cmd_main"
					;;

				true )

					cmd_list+=("$cmd_main")
					;;

			esac

		done < <(find "${peak_set_path}/${CT_MAP_ID}" -mindepth 1 -maxdepth 1 -type d -print0)


		## Jobs

		if [[ "$USE_CLUSTER" == "true" ]]; then

			job_id="p_fp_a_$(date '+%Y-%m-%d')_${algorithm:0:4}_${peak_set}"

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
micromamba activate main05

$(for cmd in "${cmd_list[@]}"; do

	echo "$cmd"

done)
exit
EOF

		fi

	done < <(find "${algorithm_path}" -mindepth 1 -maxdepth 1 -type d -print0)

done < <(find "${FOOTPRINTS_DIR}" -mindepth 1 -maxdepth 1 -type d -print0)
