#!/usr/bin/env bash

# Annotate (ct-specific) anndatas with:
#       - donor, donor_id
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

source code/glob_vars.bash # FOOTPRINTS_DIR, ATAC_PEAKS_H5AD_NEW, GROUPED_FRAG_FILES_DIR, MAIN_ENV

if [[ "$USE_CLUSTER" == "false" ]]; then

        source $HOME/.bash_profile
        load-micromamba
        micromamba activate $MAIN_ENV

fi



### SCRIPT ###

while IFS= read -r -d '' filename; do

	adata="$(basename "$filename")"
	adata_short="${adata%_raw.h5ad}"

	echo "Processing file: ${adata}"

	# Extract cell-type from anndata
	cell_type="$(h5dump -d "/obs/cell_type/categories" "$filename" |
        	awk '/^   DATA {/,/}/' |
        	grep -e '   ([0-9]*):' |  
		awk -F '"' '{ for(i=2; i<=NF; i+=2) print $i}')"

	n_cell_types=$(printf "%s\n" "$cell_type" | wc -l)

	if [[ "$n_cell_types" -ne 1 ]]; then

		echo "More than one cell type! <${cell_type}> $n_cell_types"
		exit 1

	fi


	case "$USE_CLUSTER" in

		false )

			python code/helpers/python/pre-annotate_footprint_anndata.py \
				-a "${FOOTPRINTS_DIR}/${adata}" \
				-o "${FOOTPRINTS_DIR}/${adata_short}_pre-annotated.h5ad" \
				-f "$GROUPED_FRAG_FILES_DIR/%d_${cell_type}.tsv.gz" \
				-r "$ATAC_PEAKS_H5AD_NEW" \
				-u \
				-k "grouping_col:donor,filter_col:cell_type, filter_key:${cell_type}, obs_map_col:index"
				#-F "cell_type:${cell_type}" \
			;;

		true )

			job_id="pre-annotate_footprints_$(date '+%Y-%m-%d')_${adata_short}_BORGS"
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

python code/helpers/python/pre-annotate_footprint_anndata.py \
	-a "${FOOTPRINTS_DIR}/${adata}" \
	-o "${FOOTPRINTS_DIR}/${adata_short}_pre-annotated.h5ad" \
	-f "$GROUPED_FRAG_FILES_DIR/%d_${cell_type}.tsv.gz" \
	-r "$ATAC_PEAKS_H5AD_NEW" \
	-u \
	-k "grouping_col:donor,filter_col:cell_type_custom, filter_key:${cell_type}, obs_map_col:index"
EOF
			#-F "cell_type:${cell_type}" \
			;;

		* )
			:
			;;
	
	esac

done < <(find "${FOOTPRINTS_DIR}" -mindepth 1 -maxdepth 1 \( -type f -o -type l \) -iname "*_raw.h5ad" -print0)
