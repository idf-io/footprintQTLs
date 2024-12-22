#!/usr/bin/env bash
#
# Compute atac-seq coverages of group files
#
# Arguments:
# 	cluster - str - opt - Use cluster to send jobs? [True | False] (default: False)
# 	mode - str - opt - Which data to run on.
# 			   Manually code this.
# 			   Currently [borgs | borgs_small_10k | toy]
# 	tool - str - opt - Tool used to create coverages.
# 			   [ CBN/cbn/chrompbnet | SAFT/saft/sc-atac-frag-tools wiggle_tools/wt]
# 			   (default=cbn)


### Setup

set -euo pipefail


## Environment

PROJECT_DIR=".."
PROJECT_DIR="$(realpath ${PROJECT_DIR})"
cd $PROJECT_DIR

source "code/glob_vars.bash" # DATASET, MAIN_ENV, GROUPED_GROUPED_FRAG_FILES_DIR, GROUPED_BIGWIG_FILES_DIR, TOOLS, REF_GENOME_FASTA, CHROM_SIZES


## Args

TOOL="chrombpnet"
MODE="default"
USE_CLUSTER="false"

while getopts "t:m:c" opt; do

	case "$opt" in

		t )

			case "$OPTARG" in

				CBN|cbn|chrombpnet )

					TOOL="chrombpnet"
					;;

				SAFT|saft|sc-atac-frag-tools )

					TOOL="sc-atac-fragment-tools"
					;;

				wiggle_tools|wt )

					echo "Wiggle tools not implemented yet"
					exit 1
					TOOL="wiggle_tools"
					;;

				? )

					echo "Wrong argument -t: Footprint tool [CBN | SAFT] not <$OPTARG>."
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

			echo "Usage: $0 [-c] cluster [-t] tool {cbn,saft,wt}. <${USE_CLUSTER}> not a valid argument."
			exit 1
			;;

	esac

done


case $MODE in

	default )

		if [[ "$TOOL" != "$TOOL" ]]; then

			echo "Tool in command and glob_var should be the same. <${TOOL}> != <${TOOLS}>"
			exit 1

		fi

		# DATASET
		# CT_MAP_ID
		# GROUPED_FRAG_FILES_DIR
		eval "OUT_DIR=\"$GROUPED_BIGWIG_FILES_DIR\""

		ps=""
		ms=""

		;;

	borgs_small_10k )

		DATASET="hca_brain-organoids_small_10k"
		# CT_MAP_ID
		GROUPED_FRAG_FILES_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/fragment-files/grouped/${CT_MAP_ID}"
		eval "OUT_DIR=\"data/intermediate-data/datasets/${DATASET}/atac-seq/coverages/${TOOL}/${CT_MAP_ID}\""

		ps="+4"
		ms="-5"

		;;

	toy )
		DATASET="toy_atac-seq_fragments_for_coverage"
		CT_MAP_ID="approach_placeholder" # Only used to name files
		GROUPED_FRAG_FILES_DIR="$HOME/data/toy_datasets/atac-seq_fragments_for_coverage"
		eval "OUT_DIR=\"$HOME/data/toy_datasets/atac-seq_fragments_for_coverage/coverages/${TOOL_DESCR}\""

		ps="+4"
		ms="-5"

		;;

	* )
		echo "Invalid mode argument: [default | borgs_small_10k | toy] not <$MODE>."
		exit 1
		;;

esac


## ... Environment

if [[ "$USE_CLUSTER" == "false" ]]; then

	source $HOME/.bash_profile
	load-micromamba
	micromamba activate $MAIN_ENV

	source "code/helpers/bash/compute_coverage_from_frag_file.bash"

fi


## Variables

CBN_BIN="${HOME}/.conda/envs/ian/lib/python3.9/site-packages/chrombpnet/helpers/preprocessing/reads_to_bigwig.py"



### SCRIPT ###


while IFS= read -r -d '' cell_type_dir; do

	cell_type="$(basename "$cell_type_dir")"


	## Clean old files

	if [[ -n "$(find "${OUT_DIR}/${cell_type}" -mindepth 1 -maxdepth 1 -type f -print -quit 2>/dev/null)" ]]; then

		rm -rf "${OUT_DIR}"

	fi

	mkdir -p "${OUT_DIR}/${cell_type}"


	## Cluster settings

	if [[ "$USE_CLUSTER" == "true" ]]; then

		njobs=0
		jobs_file="code/bsub/logs/compute_coverage_footprints_$(date -I)_${cell_type}.commands"

		if [[ -f "$jobs_file" ]]; then

			rm "$jobs_file"

		fi

	fi


	files_array=() # *

	while IFS= read -r -d '' frag_file; do

		echo "Processing: ${frag_file}"


		## Variables
		
		frag_name="$(basename "$frag_file")"

		if [[ "$frag_name" == *.tsv.gz ]]; then

			donor="${frag_name%.tsv.gz}"

		elif [[ "$frag_name" == *.tsv ]]; then

			donor="${frag_name%.tsv}"

		else

			echo "Wrong input format: must be .tsv or .tsv.gz"
			exit 1

		fi

		out_file="${OUT_DIR}/${cell_type}/${donor}"


		## * Avoid double processing .gz files
		
		if [[ "${#files_array[@]}" -gt 0 ]] && printf '%s\n' "${files_array[@]}" | grep -q "^${frag_file}$"; then

			echo "File ${frag_file} is already in the array. Skipping."
			continue


		else

			files_array+=("$frag_file")

			if [[ "$frag_file" == *.gz ]]; then

				files_array+=("${frag_file%.gz}")

			else

				files_array+=("${frag_file}.gz")

			fi

		fi


		cmd_cbn="frag_file_to_bw_chrombpnet \
			\"${CBN_BIN}\" \
			\"${frag_file}\" \
			\"${out_file}\" \
			\"${REF_GENOME_FASTA}\" \
			\"${CHROM_SIZES}\" \
			\"ATAC\" \
			\"${ps}\" \
			\"${ms}\""

		#scatac_fragment_tools bigwig \
			#-i 
			#-o
			#-c
			#-x \
			#-normalize False \
			#-scaling 1.0
			#–cut-sites False # Use 1 bp Tn5 cut sites (start and end of each fragment) instead of 
			# 		   whole fragment length for coverage calculation


		## Run
		
		case $USE_CLUSTER in 

			false )

				case "$TOOL" in

					chrombpnet )

						echo "$cmd_cbn"
						eval "$cmd_cbn"
						;;

					sc-atac-frag-tools )

						#TODO: implement
						echo "Not implemented yet"
						;;

				esac

				;;
				
			true )

				njobs=$(( njobs + 1 ))


				case "$TOOL" in

					chrombpnet )

						echo "$cmd_cbn" >> "$jobs_file"
						;;

					sc-atac-frag-tools )
						
						echo "Not yet implemented"
						exit 1
						;;

				esac

				;;

		esac

	done < <(find "$cell_type_dir/" -mindepth 1 -maxdepth 1 -type f -iregex ".*.tsv\(.gz\)?$" -print0)


	## Job array

	if [[ "$USE_CLUSTER" == "true" ]]; then

		job_id="compute_coverage_footprints_$(date -I)_${cell_type}"

		bsub << EOF
#!/usr/bin/env bash
#BSUB -R "rusage[mem=10G]"
#BSUB -q long
#BSUB -cwd ${PROJECT_DIR}
#BSUB -J "${job_id}"
#BSUB -o ${PROJECT_DIR}/code/bsub/logs/${job_id}.%I.out
#BSUB -e ${PROJECT_DIR}/code/bsub/logs/${job_id}.%I.err

set -euo pipefail

cd $PROJECT_DIR

source "$HOME/.bash_profile"
load-micromamba
micromamba activate $MAIN_ENV


## Run
bash <(printf "source code/helpers/bash/compute_coverage_from_frag_file.bash\n%s" "\$(cat ${jobs_file})")
EOF

	fi

done < <(find "$GROUPED_FRAG_FILES_DIR/" -mindepth 1 -maxdepth 1 -type d -print0)
