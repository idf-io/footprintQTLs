#!/usr/bin/env bash
set -euo pipefail

source "helpers/bash/split_frag_file_by_groups.bash"
source "helpers/bash/utils.bash"

source ~/.bash_profile
load-conda
conda activate ian

USE_CLUSTER="${1:-"False"}"

case $USE_CLUSTER in

	False ) echo "Local execution." ;;
	True ) echo "Cluster execution via jobs." ;;
	* ) echo "Argument 1 [True | False] not <${USE_CLUSTER}>." exit 1 ;;

esac


# Set paths
PROJECT_PATH=".."
PROJECT_PATH="$(realpath ${PROJECT_PATH})"
cd $PROJECT_PATH

MODE="${2:-"borgs"}" # hca brain organoids

# ... for different modes
case $MODE in

	borgs )
		DATASET="hca_brain-organoids"
		CT_MAP_JSON="config/cell-type_groupings/${DATASET}/approach_2024-09-12.json"
		CT_MAP_KEY="$(basename "${CT_MAP_JSON%.json}")"
		SAMPLES_DIR="data/datasets/${DATASET}/atac-seq/fragment-files"
		OUT_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/fragment-files/grouped/${CT_MAP_KEY}"
		CELLIDS_GROUP_FILES_DIR="data/intermediate-data/datasets/${DATASET}/annotations_general/group_cell-ids/${CT_MAP_KEY}"
		;;

	toy )
		DATASET="toy_atac-seq"
		CT_MAP_JSON="/home/fichtner/data/toy_datasets/atac-seq_fragments_simple/ct_map.json"
		CT_MAP_KEY="$(basename "${CT_MAP_JSON%.json}")"
		SAMPLES_DIR="/home/fichtner/data/toy_datasets/atac-seq_fragments_simple/"
		OUT_DIR="/home/fichtner/data/toy_datasets/atac-seq_fragments_simple/grouped"
		CELLIDS_GROUP_FILES_DIR="/home/fichtner/data/toy_datasets/atac-seq_fragments_simple/gci"
		;;

	borgs_small_10k )
		DATASET="hca_brain-organoids_small_10k"
		CT_MAP_JSON="config/cell-type_groupings/hca_brain-organoids/approach_2024-09-12.json"
		CT_MAP_KEY="$(basename "${CT_MAP_JSON%.json}")"
		SAMPLES_DIR="data/datasets/hca_brain-organoids_small/10k"
		OUT_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/fragment-files/grouped/${CT_MAP_KEY}"
		CELLIDS_GROUP_FILES_DIR="data/intermediate-data/datasets/hca_brain-organoids/annotations_general/group_cell-ids/${CT_MAP_KEY}"
		;;

	borgs_small_1m )
		DATASET="hca_brain-organoids_small_1m"
		CT_MAP_JSON="config/cell-type_groupings/hca_brain-organoids/approach_2024-09-12.json"
		CT_MAP_KEY="$(basename "${CT_MAP_JSON%.json}")"
		SAMPLES_DIR="data/datasets/hca_brain-organoids_small/1m"
		OUT_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/fragment-files/grouped/${CT_MAP_KEY}"
		CELLIDS_GROUP_FILES_DIR="data/intermediate-data/datasets/hca_brain-organoids/annotations_general/group_cell-ids/${CT_MAP_KEY}"
		;;

	* )
		echo "Bug in code. Check MODE variable."
		exit 1
		;;

esac

# Handle if output folder not empty
if [[ -n "$(find $OUT_DIR -maxdepth 1 -mindepth 1 -print -quit 2>/dev/null)" ]]; then

	while true; do

		read -r -p "Output folder not empty:
${OUT_DIR}
Delete contents? [Y,n]: " del

		case $del in

			[Yy]* ) rm -rf "${OUT_DIR:?}"; break ;;
			[Nn]* ) echo "I won't work with confounding folders and files."; exit 1 ;;
			* ) echo "Please answer [Y | n].";;

		esac

	done

fi

mkdir -p $OUT_DIR


# Gather group cell-ids
#python code/gather_group_cell-ids.py

# Split files per sample


# If using cluster
cluster_jobs=()

while IFS= read -r -d '' sample; do
	
	echo "Processing sample: $sample"

	if [[ $USE_CLUSTER == "False" ]]; then

		split_frag_file_by_groups \
			"$sample" \
			$OUT_DIR \
			$CELLIDS_GROUP_FILES_DIR \
			"True"

	elif [[ $USE_CLUSTER == "True" ]]; then

		# For debugging, log eval code
		split_frag_file_by_groups \
			"$sample" \
			"$OUT_DIR" \
			"$CELLIDS_GROUP_FILES_DIR" \
			"True" \
			"True" \
			"$PROJECT_PATH" > "${PROJECT_PATH}/code/bsub/logs/split_frags_bgroups_$(date '+%Y-%m-%d')_$(basename "${sample%.tsv.gz}").bash"

		bsub < <(split_frag_file_by_groups \
				"$sample" \
				"$OUT_DIR" \
				"$CELLIDS_GROUP_FILES_DIR" \
				"True" \
				"True" \
				"$PROJECT_PATH")


					
		#cluster_jobs+=("$(bsub < <(split_frag_file_by_groups \
					#"$sample" \
					#"$OUT_DIR" \
					#"$CELLIDS_GROUP_FILES_DIR" \
					#"True" \
					#"True" \
					#"$PROJECT_PATH") \
				#| cut -d ' ' -f2 | sed 's/[<>]//g')")



	else

		echo "Argument 1 [True | False], not <${USE_CLUSTER}>"
		exit 1
	fi

done < <(find "${SAMPLES_DIR}" -mindepth 1 -maxdepth 1 \( -type l -o -type f \) -name "*.tsv.gz" -print0)
