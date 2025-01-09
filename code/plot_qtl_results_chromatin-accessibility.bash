#!/usr/bin/env bash

# Plot footprintQTL results
#
# Requires: 
# 	- qtls file or folder
#
# Arguments:
# 	- mode: -m {single-tests, bulk-tests, peak-tests}
#	- cluster: -c


### SETUP ###

set -eou pipefail


## Args

MODE="unset"
MODE="bulk-tests"
USE_CLUSTER=false

while getopts "m:c" opt; do

	case "$opt" in

		#m ) MODE=$OPTARG ;;
		c ) USE_CLUSTER=true ;;
		? ) echo "Usage: $0 [-c] [-m mode {bulk-tests,single-tests}"; exit 1 ;;
	
	esac

done

case "$MODE" in

	single-tests | bulk-tests | peak-tests )

		:
		;;

	unset )

		echo "Usage: $0 [-c] [-m mode {bulk-tests,single-tests}"
		exit 1
		;;

	? )

		echo "Usage: $0 [-c] [-m mode {bulk-tests,single-tests}"
		exit 1
		;;

esac


## Environment

PROJECT_DIR=".."
PROJECT_DIR="$(realpath ${PROJECT_DIR})"
cd $PROJECT_DIR

source code/glob_vars.bash # MATRIX_EQTL_OUTPUT_DIR, CT_MAP_ID

if [[ "$USE_CLUSTER" == "false" ]]; then

	source "${HOME}/.bash_profile"
	load-r-440

fi


## User variables

ALPHA=0.01 # Aygun-2023 caQTLS: 0.05



### SCRIPT ###

while IFS= read -r -d '' cell_type_path; do

	cell_type="$(basename "$cell_type_path")"

	echo -e "\t\t$cell_type"


	in_dir="${cell_type_path}/${MODE}"
	out_dir="${cell_type_path}/${MODE}/plots"


	## Clean old run files and create dir
	
	if [[ -f "$out_dir" ]]; then
	
		rm -r "${in_dir}/qtls_all_fdr.tsv" 
	
	fi

	if [[ -d "$out_dir" ]]; then
	
		rm -rf "$out_dir" 
	
	fi
	
	mkdir -p "$out_dir"


	if [[ "$USE_CLUSTER" == "true" ]]; then

		njobs=0
		jobs_file="code/bsub/logs/plot_ca-qtls_$(date -I)_${cell_type}.commands"

	fi

	## Run function

	case "$MODE" in

		bulk-tests )

			in_file="${in_dir}/tests/qtls_all.tsv"

			cmd_main="Rscript --verbose 'code/plot_qtl_results.R' \
					'$in_file' \
					'$out_dir' \
					'$MODE' \
					'$ALPHA'"

			case "$USE_CLUSTER" in

				false )

					eval "$cmd_main"

					;;

				true )

					njobs=$(( njobs + 1 ))
					echo $cmd_main >> "$jobs_file"
					;;

			esac

			;;


		single-tests | peak-tests )
		
			## Create collected qtls.tsv

			in_file="${in_dir}/qtls_all.tsv"
			
			# Clean old file
			if [[ -f "$in_file" ]]; then
			
				rm "$in_file" 
			
			fi

			tests_bed="${cell_type_path}/${MODE}/tests_snp_peak_pairs.bed"

			cat "${in_dir}/tests/"* > "$in_file"

			n_significant="$(cat "$in_file" | wc -l)"
			n_tests="$(cat "$tests_bed" | wc -l)"
			n_missing_tests=$((n_tests - n_significant))

			# snp, gene, fdr_meqtl, pvalue, beta,
			# snp_chr, snp_pos, snp_ref, snp_alt, peak_chr, peak_start, peak_end, peak_len
			non_significant_placeholder=$'placeholder_snp\tplaceholder_peak\t1.0\t1.0\t0.0\t0\t0\tN\tN\t0\t0\t0\t501'


			cmd_make_in="seq $n_missing_tests | xargs -I {} echo '$non_significant_placeholder' >> '$in_file'"

			cmd_main="Rscript --verbose 'code/plot_qtl_results.R' \
					'$in_file' \
					'$out_dir' \
					'$MODE' \
					'$ALPHA'"


			case "$USE_CLUSTER" in

				false )

					eval "$cmd_make_in"
					eval "$cmd_main"
					;;


				true )

					njobs=$(( njobs + 1 ))
					echo "${cmd_make_in} && ${cmd_main}" >> "$jobs_file"
					;;

			esac

			;;

	esac

	## Job array

	if [[ "$USE_CLUSTER" == "true" ]]; then

		job_id="plot_fqtls_$(date -I)_${cell_type}"

		bsub << EOF
#!/usr/bin/env bash
#BSUB -R "rusage[mem=5G]"
#BSUB -q long
#BSUB -cwd ${PROJECT_DIR}
#BSUB -J "${job_id}"
#BSUB -o ${PROJECT_DIR}/code/bsub/logs/${job_id}.%I.out
#BSUB -e ${PROJECT_DIR}/code/bsub/logs/${job_id}.%I.err

set -euo pipefail

cd $PROJECT_DIR

source "$HOME/.bash_profile"
load-r-440


## Run

while IFS= read -r -d '' line; do

	eval "\$line"
	#sed -n "\${LSB_JOBINDEX}p" ${jobs_file} | bash

done < <(printf '%s\\0' "\$(cat ${jobs_file})")
EOF

	fi

done < <(find "${MATRIX_EQTL_OUTPUT_DIR}/chromatin-accessibility/${CT_MAP_ID}" -mindepth 1 -maxdepth 1 -type d -print0)
