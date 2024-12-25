#!/usr/bin/env bash

# Annotate (ct-specific) anndatas with:
#    - donor, donor_id
# 	- nr. of cells per donor-ct
# 	- nr. of insertions per donor-ct
# 	- nr. of fragments per donor-ct
# 	- locations of peaks
# 	- metadata from reference anndata
#
# Requires:
# 	- anndata
# 	- path with donor-ct grouped fragment files, to extract:
# 		- n_insertions per donor
# 		- n_fragmetns per donor
# 		- n_cells per donor: 1st option to extract this data
# 		- path format: "%d_myDS.tsv.gz" where %d represents the donors that will be matched t the anndata obs
# 	- reference anndata with other .obs and .vars annotations, to extract:
# 		- n_cells per donor, 2nd option to extract this data
# 		- var and obs metadata
#
# Assumptions:
# 	- anndata object contains an .obs column with the ct that each observation(donor) belongs to.
# 	- in requirement 2: donor-ct grouped fragment files, all files which match *ct* (1 lvl deep) are used to match cell-type from footprint anndata and grouped fragment-file anndata

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

source code/glob_vars.bash # FOOTPRINTS_DIR, ATAC_PEAKS_PROCESSED_H5AD, GROUPED_FRAG_FILES_DIR, MAIN_ENV, CT_MAP_ID

if [[ "$USE_CLUSTER" == "false" ]]; then

        source $HOME/.bash_profile
        load-micromamba
        micromamba activate $MAIN_ENV

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


			cmd_main="python code/helpers/python/pre-annotate_footprint_anndata.py \
						-a '${FOOTPRINTS_DIR}/${algorithm}/${peak_set}/${CT_MAP_ID}/${cell_type}/footprints_raw.h5ad' \
						-o '${FOOTPRINTS_DIR}/${algorithm}/${peak_set}/${CT_MAP_ID}/${cell_type}/footprints_pre-annotated.h5ad' \
						-f '$GROUPED_FRAG_FILES_DIR/${cell_type}/%d.tsv.gz' \
						-r '$ATAC_PEAKS_PROCESSED_H5AD' \
						-u \
						-k 'grouping_col:donor,filter_col:cell_type, filter_key:${cell_type}, obs_map_col:index'"
						
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
		
			job_id="pre-annotate_footprints_$(date '+%Y-%m-%d')_${algorithm:0:4}_${peak_set}"

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

$(for cmd in "${cmd_list[@]}"; do

	echo "$cmd"

done)
EOF

		fi

	done < <(find "${algorithm_path}" -mindepth 1 -maxdepth 1 -type d -print0)

done < <(find "${FOOTPRINTS_DIR}" -mindepth 1 -maxdepth 1 -type d -print0)


