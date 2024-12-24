#!/usr/bin/env bash
#
# Compute footprints of grouped fragment files
#
# Arguments:
#   - cluster (str): Optional. Use cluster to send jobs? {True, False} (default: False)
#   - mode (str): Optional. Which data to run on.
#		Manually coded if not using the glob_env variable DATASET
#		{default, borgs_small_10k, toy}
#		(default = default)
#
# Notes:
# 	- Peak files in peak files directory will only be considered if starting with "peaks_"
#

### Setup ##


set -euo pipefail


## Environment

PROJECT_DIR=".."
PROJECT_DIR="$(realpath ${PROJECT_DIR})"
cd $PROJECT_DIR

source "code/glob_vars.bash" # MAIN_ENV, DATASET, CT_MAP_ID,
							 # GROUPED_BIGWIG_FILES_DIR, FOOTPRINTS_DIR,
							 # SELECT_PEAKS_TSV_DIR, ALGORITHMS, 


## Args

ALGORITHM="js_divergence"                                                   
MODE="default"
USE_CLUSTER="false"                                                                                              
                   
while getopts "a:m:c" opt; do
             
    case "$opt" in       
                    
        t )                                            
                                                       
            case "$OPTARG" in
                                            
                js_divergence )
               
                    ALGORITHM="chrombpnet"
                    ;;      
                                                                                                                  
                counts )
                           
                    ALGORITHM="sc-atac-fragment-tools"
                    ;;
                                                                                                         
                wasserstein_dist )     

                    echo "Wasserstein distance not implemented yet"
                    exit 1
                    ALGORITHM="wiggle_tools"                                          
                    ;;

                ? )

                    echo "Wrong argument -a: Algorithm [js_divergence | counts | wasserstein_dist] not <$OPTARG>."
                    exit 1
                    ;;

            esac
            ;;

        m )

            MODE="$OPTARG"
            ;;

        c )

            USE_CLUSTER="true"
            ;;

        ? )

            echo "Usage: $0 [-c] cluster [-a] algorithm {js_divergence, counts, wasserstein_dist}. Not a valid argument."
            exit 1
            ;;

    esac

done




## ... Environment

if [[ "$USE_CLUSTER" == "false" ]]; then

    source $HOME/.bash_profile
    load-micromamba
    micromamba activate $MAIN_ENV

fi


case $MODE in

	default )
        # DATASET
        # CT_MAP_ID
        # GROUPED_FRAG_FILES_DIR
		COV_FILES_DIR="$GROUPED_BIGWIG_FILES_DIR"
		OUT_DIR="${FOOTPRINTS_DIR}"
		PEAK_FILES_DIR="${SELECT_PEAKS_TSV_DIR}"
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



### Script ###

while IFS= read -r -d '' cell_type_dir; do 

	cell_type="$(basename "$cell_type_dir")"

	echo "${cell_type}"

	if [[ "$cell_type" == "Discard" ]]; then

		continue

	fi


	while IFS= read -r -d '' peak_set_file_path; do

		peak_set_file="$(basename "$peak_set_file_path")"
		peak_set="${peak_set_file%.bed}"
		peak_set="${peak_set/peaks_/}"

		echo -e "\t${peak_set}"



		for algorithm in "${ALGORITHMS[@]}"; do

			echo -e "\t\t$algorithm"


			## Vars

			out_file="${OUT_DIR}/${algorithm}/${peak_set}/${CT_MAP_ID}/${cell_type}/footprints_raw.h5ad"

			mkdir -p "$(dirname "$out_file")"


			cmd_main="python code/helpers/python/compute_footprints.py \
						-c '${cell_type_dir}' \
						-o '${out_file}' \
						-p '${peak_set_file_path}' \
						-a '${algorithm}' \
						-A 'cell_type,${cell_type}'"

			cmd_main="$(echo "$cmd_main")"


			## Run

			case $USE_CLUSTER in

				false )

					eval "$cmd_main"
					;;

				true )

					job_id="compute_fp_$(date '+%Y-%m-%d')_${peak_set}_${algorithm:0:4}_${cell_type}"
					bsub <<EOF
#!/usr/bin/env bash
#BSUB -R "rusage[mem=30G]"
#BSUB -q medium
#BSUB -cwd ${PROJECT_DIR}
#BSUB -J ${job_id}
#BSUB -o ${PROJECT_DIR}/code/bsub/logs/${job_id}.out
#BSUB -e ${PROJECT_DIR}/code/bsub/logs/${job_id}.err

set -euo pipefail

cd ${PROJECT_DIR}

source "${HOME}/.bash_profile"
load-micromamba
micromamba activate ${MAIN_ENV}

${cmd_main}
EOF
					;;
			
			esac

		done

	done < <(find "${PEAK_FILES_DIR}/${cell_type}" -mindepth 1 -maxdepth 1 \( -type f -o -type l \) -iname "peaks_*" -print0)

done < <(find "${COV_FILES_DIR}" -mindepth 1 -maxdepth 1 -type d -print0)
