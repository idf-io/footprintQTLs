#!/usr/bin/env bash
# Description: For each fragment file in a folder create the corresponding coverages of the pooled fragments.
# Arguments:
# 	cluster - str - opt - Use cluster to send jobs? [True | False] (default: False)
# 	mode - str - opt - Which data to run on.
# 			   Manually code this.
# 			   Currently [borgs | borgs_small_10k | toy]
# 	tool - str - opt - Tool used to create coverages.
# 			   [ CBN/cbn/chrompbnet | SAFT/saft/sc-atac-frag-tools wiggle_tools/wt]
# 			   (default=cbn)
set -euo pipefail

PROJECT_PATH=".."
PROJECT_PATH="$(realpath ${PROJECT_PATH})"
cd $PROJECT_PATH

source "code/helpers/bash/compute_coverage_from_frag_file.bash"

source $HOME/.bash_profile
load-conda
conda activate ian

USE_CLUSTER="${1:-false}"

case $USE_CLUSTER in

	false ) echo "Local execution." ;;
	true ) echo "Cluster execution via jobs." ;;
	* ) echo "Argument 1 [True | False] not <${USE_CLUSTER}>."; exit 1 ;;

esac


TOOL="${3:-"CBN"}" # ChromBPNet

case "$TOOL" in

	CBN|cbn|chrombpnet )
		TOOL_DESCR="chrombpnet"
		;;

	SAFT|saft|sc-atac-frag-tools )
		TOOL_DESCR="sc-atac-fragment-tools"
		;;

	wiggle_tools|wt )
		echo "Wiggle tools not implemented yet"
		exit 1
		TOOL_DESCR="wiggle_tools"
		;;

	* ) echo "Wrong argument 3: Footprint tool [CBN | SAFT] not <$TOOL>."; exit 1 ;;

esac

CBN_BIN="${HOME}/.conda/envs/ian/lib/python3.9/site-packages/chrombpnet/helpers/preprocessing/reads_to_bigwig.py"
MODE="${2:-"borgs"}" # hca brain organoids


case $MODE in

	borgs )
		DATASET="hca_brain-organoids"
		CT_MAP_JSON="config/cell-type_groupings/${DATASET}/approach_2024-09-12.json"
		CT_MAP_ID="$(basename "${CT_MAP_JSON%.json}")"
		FRAG_FILES_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/fragment-files/grouped/${CT_MAP_ID}"
		eval "OUT_DIR=\"data/intermediate-data/datasets/${DATASET}/atac-seq/coverages/${TOOL_DESCR}/${CT_MAP_ID}\""

		ps=""
		ms=""

		;;

	borgs_small_10k )

		DATASET="hca_brain-organoids_small_10k"
		CT_MAP_JSON="config/cell-type_groupings/hca_brain-organoids/approach_2024-09-12.json"
		CT_MAP_ID="$(basename "${CT_MAP_JSON%.json}")"
		FRAG_FILES_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/fragment-files/grouped/${CT_MAP_ID}"
		eval "OUT_DIR=\"data/intermediate-data/datasets/${DATASET}/atac-seq/coverages/${TOOL_DESCR}/${CT_MAP_ID}\""

		ps="+4"
		ms="-5"

		;;

	toy )
		DATASET="toy_atac-seq_fragments_for_coverage"
		CT_MAP_ID="approach_placeholder" # Only used to name files
		FRAG_FILES_DIR="$HOME/data/toy_datasets/atac-seq_fragments_for_coverage"
		eval "OUT_DIR=\"$HOME/data/toy_datasets/atac-seq_fragments_for_coverage/coverages/${TOOL_DESCR}\""

		ps="+4"
		ms="-5"

		;;

	* )
		echo "Wrong argument 2: MODE [borgs | borgs_small_10k | toy] not <$MODE>."
		exit 1
		;;

esac



## IO setup

# Output already contains files or dirs
if [[ -n "$(find "$OUT_DIR" -mindepth 1 -maxdepth 1 \( -type f -o -type l \) -print -quit 2>/dev/null)" ]]; then

	while true; do

		echo "Output directory <${OUT_DIR}> already contains files or folders."
		read -r -p "Delete contents? [Y,n]: " del

		case "$del" in

			[Yy]* ) rm -rf "${OUT_DIR:?}"; break ;;
			[Nn]* ) echo "I won't work with confounding folders or files."; exit 1 ;;
			* ) echo "Please answer [Y | n]." ;;

		esac

	done

fi

mkdir -p "$OUT_DIR"


## Main



while IFS= read -r -d '' frag_file; do

	echo "Processing: ${frag_file}"

	frag_name="$(basename "$frag_file")"
	frag_name_short="${frag_name%.tsv.gz}"

	main_cmd="frag_file_to_bw_chrombpnet \
		\"${CBN_BIN}\" \
		\"${FRAG_FILES_DIR}/\${frag_name}\" \
		\"${OUT_DIR}/\${frag_name_short}\" \
		\"data/GRCh38-p14/hg38.fa\" \
		\"data/GRCh38-p14/hg38.chrom.sizes\" \
		\"ATAC\" \
		\"${ps}\" \
		\"${ms}\""

	case $USE_CLUSTER in 

		false )

			case $TOOL in

				CBN|cbn|chrombpnet )
					eval "$main_cmd"
					;;

				SAFT|saft|sc-atac-frag-tools )
					#TODO: implement
					echo "Not implemented yet"
					#scatac_fragment_tools bigwig \
						#-i 
						#-o
						#-c
						#-x \
						#-normalize False \
						#-scaling 1.0
						#–cut-sites False # Use 1 bp Tn5 cut sites (start and end of each fragment) instead of 
						# 		   whole fragment length for coverage calculation
					;;

				* )
					echo "Error in logic"
					exit 1
					;;

			esac

			;;
			
		true )

			case $TOOL in

				CBN|cbn|chrombpnet )

					main_cmd="$main_cmd \
						\"true\" \
						\"$PROJECT_PATH\""

					job_script_path="${PROJECT_PATH}/code/bsub/logs/compute_coverage_$(date '+%Y-%m-%d')_${frag_name_short}.bash"

					eval "$main_cmd" > "$job_script_path"

					job_id="$(bsub < "$job_script_path" | cut -d' ' -f2 | sed 's/[<>]//g')"

					cluster_jobs+=("$job_id")
					;;

				SAFT|saft|sc-atac-frag-tools )
					echo "Not yet implemented"
					exit 1
					;;

				* )
					"Logic mistake"
					exit 1
					;;

			esac

			;;
			
		* )
			echo "Wrong logic"
			exit 1
			;;

	esac

done < <(find "$FRAG_FILES_DIR/" -mindepth 1 -maxdepth 1 \( -type f -o -type l \) -iregex ".*.tsv\(.gz\)?$" -print0)
