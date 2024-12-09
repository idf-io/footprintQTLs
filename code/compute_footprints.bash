#!/usr/bin/env bash
set -euo pipefail

source "$HOME/.bash_profile"
load-conda
conda activate ian

PROJECT_PATH=".."
PROJECT_PATH="$(realpath ${PROJECT_PATH})"
cd $PROJECT_PATH

USE_CLUSTER="${1:-false}"

case $USE_CLUSTER in

	false ) echo "Local execution." ;;
	true ) echo "Cluster execution via jobs." ;;
	* ) echo "Wrong argument 1: use_cluster [true | false] not <${USE_CLUSTER}>."; exit 1 ;;

esac

MODE="${2:-borgs}" # hca brain organoids

case $MODE in

	borgs )
		DATASET="hca_brain-organoids"
		CT_MAP_JSON="config/cell-type_groupings/${DATASET}/approach_2024-09-12.json"
		CT_MAP_KEY="$(basename "${CT_MAP_JSON%.json}")"

		FOOTPRINT_APPROACH="js_divergence"
		COV_FILES_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/coverages/chrombpnet/${CT_MAP_KEY}/low_mem_50G_28-10" # TODO: correct if necessary
		OUT_DIR="results/datasets/${DATASET}/atac-seq/footprints/${FOOTPRINT_APPROACH}/${CT_MAP_KEY}"
		PEAK_FILES_DIR="data/intermediate-data/datasets/${DATASET}/matrix-eQTL_io_5caPCs/chromatin_accessibility/peak_ca/${CT_MAP_KEY}"
		;;

	borgs_small_10k )
		echo "Current implementation runs equally as fast as with all initial reads AKA borgs. 10k is therefore unnecessary. For a better ~downsampled~ version of the borgs dataset, define a new approach."
		DATASET="hca_brain-organoids_small_10k"
		CT_MAP_JSON="config/cell-type_groupings/hca_brain-organoids/approach_2024-09-12.json"
		CT_MAP_KEY="$(basename "${CT_MAP_JSON%.json}")"

		FOOTPRINT_APPROACH="js_divergence"
		COV_FILES_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/coverages/chrombpnet/${CT_MAP_KEY}"
		OUT_DIR="results/datasets/${DATASET}/atac-seq/footprints/${FOOTPRINT_APPROACH}/${CT_MAP_KEY}"
		PEAK_FILES_DIR="data/intermediate-data/datasets/hca_brain-organoids/matrix-eQTL_io_5caPCs/chromatin_accessibility/peak_ca/${CT_MAP_KEY}"
		;;

	toy )
		FOOTPRINT_APPROACH="js_divergence"
		COV_FILES_DIR="${HOME}/data/toy_datasets/atac-seq_coverages_for_footprinting/coverages"
		OUT_DIR="${HOME}/data/toy_datasets/atac-seq_coverages_for_footprinting/footprints"
		PEAK_FILES_DIR="${HOME}/data/toy_datasets/atac-seq_coverages_for_footprinting/peak_refs"
		;;

	* )
		echo "Wrong argument 2: mode [borgs | borgs_small_10k | toy] not <${MODE}>."
		exit 1
		;;

esac

case $FOOTPRINT_APPROACH in

	jensen-shannon_divergece | js_divergence | js_div | jsd )
		# Jensen-Shannon divergence
		fp_approach="js_div"
		;;

	* )
		echo "Footprint computation approach not implemented: [js_divergence] not <$FOOTPRINT_APPROACH>"
		exit 1
		;;

esac


## Script ##

# ASSUMPTION:
#   - COV_FILES_DIR:depth1 is only populated with directories (representing experimental groups) containing fragment_files.
#   - Names of dirs in COV_FILES_DIR are the same as the title of the respective peak files in .tsv.

while IFS= read -r -d '' condition_dir; do 

	condition_short="$(basename "$condition_dir")"

	if [[ "$condition_short" == "Discard" ]]; then

		continue

	fi

	echo "Processing condition: $condition_dir"

	case $USE_CLUSTER in

		false )
			python code/helpers/python/compute_footprints.py \
				-c "$condition_dir" \
				-o "$OUT_DIR/footprints_${condition_short}_raw.h5ad" \
				-p "$PEAK_FILES_DIR/$condition_short.tsv" \
				-H \
				-a "$FOOTPRINT_APPROACH"
				-A "cell_type,${condition_short}"
			;;

		true )

			job_id="compute_footprints_$(date '+%Y-%m-%d')_${condition_short}_BORGS"
			bsub <<EOF
#!/usr/bin/env bash
#BSUB -R "rusage[mem=30G]"
#BSUB -q medium
#BSUB -cwd ${PROJECT_PATH}
#BSUB -J ${job_id}
#BSUB -o ${PROJECT_PATH}/code/bsub/logs/${job_id}.out
#BSUB -e ${PROJECT_PATH}/code/bsub/logs/${job_id}.err

set -euo pipefail

cd $PROJECT_PATH

source "$HOME/.bash_profile"
load-conda
conda activate ian

python code/helpers/python/compute_footprints.py \
	-c "$condition_dir" \
	-o "$OUT_DIR/footprints_${condition_short}_raw.h5ad" \
	-p "$PEAK_FILES_DIR/$condition_short.tsv" \
	-H \
	-a "$FOOTPRINT_APPROACH" \
	-A "cell_type,${condition_short}"
EOF
			;;

		* )
			echo "Wrong logic."
			exit 1
			;;
	
	esac

done < <(find "${COV_FILES_DIR}/" -mindepth 1 -maxdepth 1 -type d -print0)
