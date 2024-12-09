# Process ct-specific footprint anndatas by calculating the descriptive stats on var-level, 
# highly variable peaks, PCA, UMAP and clustering.
#
# Requires:
# 	- footprint anndatas 
# 		- adata.obs['n_cells'] annotation

### SETUP ###
set -euo pipefail

# Get args

USE_CLUSTER=false

while getopts "c" opt; do

	case "$opt" in

		c ) USE_CLUSTER=true ;;
		? ) echo "Usage: $0 [-c]"; exit 1 ;;
	
	esac

done

# Environment

PROJECT_PATH=".."
PROJECT_PATH="$(realpath "$PROJECT_PATH")"
cd $PROJECT_PATH

source code/glob_vars.bash # FOOTPRINTS_DIR, MAIN_ENV


case "$USE_CLUSTER" in

	false )
        	source $HOME/.bash_profile
        	load-micromamba
        	micromamba activate $MAIN_ENV
		;;

	true )
		njobs=0

		jobs_file="code/bsub/logs/process_footprint_adata_$(date +%F)_BORGS.commands"

		if [[ -f "$jobs_file" ]]; then

			echo "jobs_file <$jobs_file> already exists. Delete and proceed? [y/n]"
			read -p "> " delete

			case "$delete" in

				[Yy]* )
					rm "$jobs_file"
					;;

				[Nn]* )
					echo "Cannot proceed then."
					exit 1
					;;

				*)
					echo "Wrong input"
					exit 1
					;;
			esac

		fi

		touch "$jobs_file"

		;;

esac


### SCRIPT ###

# Loop over anndata files
while IFS= read -r -d '' adata_full; do

	adata="$(basename ${adata_full})"
	adata_short="${adata%_pre-annotated.h5ad}"
	adata_out="${FOOTPRINTS_DIR}/${adata_short}_processed.h5ad"

	echo "Processing file: ${adata_full}"

	case "$USE_CLUSTER" in

		false )
			python code/helpers/python/process_footprint_adata.py \
				-a "$adata_full" \
				-o "$adata_out" \
				-m "200"
			;;

		true )
			njobs=$(( $njobs + 1 ))

			cmd="
			python code/helpers/python/process_footprint_adata.py \
				-a \"$adata_full\" \
				-o \"${adata_out}\" \
				-m \"200\""

			echo $cmd >> "$jobs_file"
			;;

		* )
			:
			;;

	esac

done < <(find ${FOOTPRINTS_DIR} -mindepth 1 -maxdepth 1 \( -type f -o -type l \) -iname "*_pre-annotated.h5ad" -print0)


# Job array
if [[ "$USE_CLUSTER" == "true" ]]; then

	job_id="process_footprint_adata_$(date +%F)_"
	bsub << EOF
#!/usr/bin/env bash
#BSUB -R "rusage[mem=10G]"
#BSUB -q medium
#BSUB -cwd ${PROJECT_PATH}
#BSUB -J "${job_id}[1-$njobs]"
#BSUB -o ${PROJECT_PATH}/code/bsub/logs/${job_id}.%I.out
#BSUB -e ${PROJECT_PATH}/code/bsub/logs/${job_id}.%I.err

set -euo pipefail

cd $PROJECT_PATH

source "$HOME/.bash_profile"
load-micromamba
micromamba activate $MAIN_ENV

sed -n "\${LSB_JOBINDEX}p" ${jobs_file} | bash
echo "Finished script"
EOF

fi
