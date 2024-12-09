#!/usr/bin/env bash

# Join files within each donor_ct directory (contains sample%donor_ct.tsv.gz files)
#
# Input:
# 	- CT_MAP_KEY: json of cell-type grouping scheme
# 	- SAMPLS_DIR: directory of donor_ct directories with sample%donor_ct.tsv.gz to aggregate to donor_ct.tsv.gz
# 	- OUT_DIR: directory to create the donor_ct.tsv.gz files in
#
# Output:
# 	- donor_ct.tsv.gz files

set -euo pipefail

USE_CLUSTER="${1:-False}"

case $USE_CLUSTER in

	False ) echo "Local execution." ;;
	True ) echo "Cluster execution via jobs." ;;
	* ) echo "Argument 1 [True | False] not <${USE_CLUSTER}>."; exit 1 ;;

esac


# Set paths
PROJECT_PATH=".."
PROJECT_PATH="$(realpath ${PROJECT_PATH})"
cd $PROJECT_PATH

source "code/helpers/bash/utils.bash"
source "code/helpers/bash/join_frag_files_by_groups.bash"

MODE="${2:-borgs}" # hca brain organoids

case $MODE in

	borgs )
		DATASET="hca_brain-organoids"
		CT_MAP_JSON="config/cell-type_groupings/${DATASET}/approach_2024-09-12.json"
		CT_MAP_KEY="$(basename "${CT_MAP_JSON%.json}")"
		SAMPLES_DIR="data/datasets/${DATASET}/atac-seq/fragment-files"
		OUT_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/fragment-files/grouped/${CT_MAP_KEY}"
		# CELLIDS_GROUP_FILES_DIR="data/intermediate-data/datasets/${DATASET}/annotations_general/group_cell-ids/${CT_MAP_KEY}"
		;;

	toy )
		DATASET="toy_atac-seq"
		CT_MAP_JSON="/home/fichtner/data/toy_datasets/atac-seq_fragments_simple/ct_map.json"
		CT_MAP_KEY="$(basename "${CT_MAP_JSON%.json}")"
		SAMPLES_DIR="/home/fichtner/data/toy_datasets/atac-seq_fragments_simple/"
		OUT_DIR="/home/fichtner/data/toy_datasets/atac-seq_fragments_simple/grouped"
		# CELLIDS_GROUP_FILES_DIR="/home/fichtner/data/toy_datasets/atac-seq_fragments_simple/gci"
		;;

	borgs_small_10k )
		DATASET="hca_brain-organoids_small_10k"
		CT_MAP_JSON="config/cell-type_groupings/hca_brain-organoids/approach_2024-09-12.json"
		CT_MAP_KEY="$(basename "${CT_MAP_JSON%.json}")"
		SAMPLES_DIR="data/datasets/hca_brain-organoids_small/10k"
		OUT_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/fragment-files/grouped/${CT_MAP_KEY}"
		# CELLIDS_GROUP_FILES_DIR="data/intermediate-data/datasets/hca_brain-organoids/annotations_general/group_cell-ids/${CT_MAP_KEY}"
		;;

	borgs_small_1m )
		DATASET="hca_brain-organoids_small_1m"
		CT_MAP_JSON="config/cell-type_groupings/hca_brain-organoids/approach_2024-09-12.json"
		CT_MAP_KEY="$(basename "${CT_MAP_JSON%.json}")"
		SAMPLES_DIR="data/datasets/hca_brain-organoids_small/1m"
		OUT_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/fragment-files/grouped/${CT_MAP_KEY}"
		# CELLIDS_GROUP_FILES_DIR="data/intermediate-data/datasets/hca_brain-organoids/annotations_general/group_cell-ids/${CT_MAP_KEY}"
		;;

	* )
		echo "Bug in code. Check MODE variable."
		exit 1
		;;

esac

cluster_jobs=()

while IFS= read -r -d '' group_dir; do

	group_name="$(basename "$group_dir")"
	final_out="${OUT_DIR}/${group_name}.tsv"

	echo "Gathering all fragments of group: ${group_name}"

	case $USE_CLUSTER in

		False )
			join_frag_files_by_groups \
				"$group_dir" \
				"$final_out"
			;;

		True )
			job_script_path="${PROJECT_PATH}/code/bsub/logs/join_frags_bgroups_$(date '+%Y-%m-%d')_${group_name}.bash"

			join_frag_files_by_groups \
				"$group_dir" \
				"$final_out" \
				"True" \
				"$PROJECT_PATH" > "$job_script_path"

			job_id="$(bsub < "$job_script_path" | cut -d' ' -f2 | sed 's/[<>]//g')"

			cluster_jobs+=("$job_id")

			;;

		* )
			echo "Argument 1 <use cluster?> must be [True | False] not <${USE_CLUSTER}"
			exit 1
			;;
	
	esac

done < <(find "$OUT_DIR/tmp/" -mindepth 1 -maxdepth 1 -type d -print0)

# If using cluster, wait until jobs have finished

if [[ $USE_CLUSTER == "True" ]]; then

	echo "n_jobs: ${#cluster_jobs[@]}"
	echo "join_frag_files_by_groups"
	echo "Waiting for jobs to finish..."

	if [ ${#cluster_jobs[@]} -eq 0 ]; then

        	echo "No jobs were submitted. Exiting."
		exit 1

	fi

        while true; do

                if check_jobs "${cluster_jobs[@]}"; then

                        echo "All jobs finished."
                        break

                else

                        sleep $(expr 60 \* 1)

                fi

        done

fi
