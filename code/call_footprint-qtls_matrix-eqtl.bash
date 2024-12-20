#!/usr/bin/env bash
#
# Call footprintQTLs using matrix-eQTL on all cell-types
#
# Requires: 
#	- make_matrix-eqtl_input.ipynb output directory
# 
# Notes:
# 	- using the cluster option [-c] with a low $MULTPLE number yields a high nr of jobs.
# 	  This can be problematic for renv and give errors when loading. Always check for success of all jobs.
# 	- Correctly setup renv, renv.lock, .Rprofile, .Renviron
# 	  which contain the matrix-eQTL data for those categories.


### SETUP ###

set -eou pipefail


## Args

MODE="unset"
USE_CLUSTER=false

while getopts "m:c" opt; do

	case "$opt" in

		m ) MODE=$OPTARG ;;
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

source code/glob_vars.bash # MATRIX_EQTL_INPUT_DIR, MATRIX_EQTL_OUTPUT_DIR
source code/helpers/bash/matrix-eqtl.bash # format_check_meqtl_io

if [[ "$USE_CLUSTER" == "false" ]]; then

	source "${HOME}/.bash_profile"
	load-r-440

fi


## User variables

ALPHA=0.01 # Aygun-2023 caQTLS: 0.05
CIS_DIST=0 # Aygun-2023 caQTLS: within peak
MULTIPLE=5000 # Nr of test-batches to perform for every job
			 # For bulk-tests: All tests in 1 tests-batch
			 # For single-tests: 1 test in 1 tests-batch
			 # For peak-tests: ~2 tests in tests-batch


### SCRIPT ###


# Run function for all cell-types

while IFS= read -r -d '' cell_type_dir; do

	echo -e "\nProcessing file: ${cell_type_dir}"

	echo -e "\t- Multiples: ${MULTIPLE}"

	cell_type="$(basename "$cell_type_dir")"
	in_dir="${cell_type_dir}/${MODE}"
	out_dir="${MATRIX_EQTL_OUTPUT_DIR}/${cell_type}/${MODE}"
	tests_bed="${MATRIX_EQTL_INPUT_DIR}/${cell_type}/${MODE}/tests_snp_peak_pairs.bed"
	job_id="call_fqtls_$(date -I)_${cell_type}"


	## Clean old run files and create dir
	
	if [[ -d "${out_dir}/tests" ]]; then

		rm -rf "${out_dir}/tests" 

	fi

	if [[ -f "${out_dir}/qtl-testing_stats.tsv" ]]; then

		rm "${out_dir}/qtl-testing_stats.tsv" 

	fi

	mkdir -p "$out_dir"


	## Cluster relevant variables

	if [[ "$USE_CLUSTER" == "true" ]]; then

		njobs=0
		jobs_file="code/bsub/logs/call_fqtls_$(date -I)_${cell_type}.commands"

		if [[ -f "$jobs_file" ]]; then

			rm "$jobs_file"

		fi

	fi


	## Iterate over cell-types

	case "$MODE" in

		bulk-tests )

			SLICE_SIZE=2000 # Default in matrix-eQTL
			RANGE_LOWER=1 # Only 1 set of tests in bulk
			RANGE_UPPER=1 # Only 1 set of tests in bulk
			FORMAT_CHECK="TRUE"

			cmd_main="Rscript --verbose 'code/call_qtls_matrix-eqtl.R' \
					'$in_dir' \
					'$out_dir' \
					'qtls_all.tsv' \
					'$MODE' \
					'$ALPHA' \
					'$CIS_DIST' \
					'$SLICE_SIZE' \
					'$RANGE_LOWER' \
					'$RANGE_UPPER' \
					'$FORMAT_CHECK'"

			case "$USE_CLUSTER" in

				false )

					## Run
					eval "$cmd_main"

					;;

				true )
					njobs=$(( njobs + 1 ))
					echo $cmd_main >> "$jobs_file"
					;;

			esac

			;;


		single-tests )

			SLICE_SIZE=1 # Only 1 test at a time
			FORMAT_CHECK="FALSE"


			## Check input format

			format_check_meqtl_io \
				"$in_dir/genotype_NA_source.tsv" \
				"$in_dir/phenotype_source.tsv" \
				"$in_dir/covariates_source.tsv"


			## Iterate over user defined multiples of tests
			
			multiple=$MULTIPLE
			n_tests=$(wc -l "$tests_bed")
			
			RANGE_LOWER=1
			RANGE_UPPER=$multiple

			while true; do

				cmd_main="Rscript --verbose 'code/call_qtls_matrix-eqtl.R' \
							'$in_dir' \
							'$out_dir' \
							'snp-level_tests_${RANGE_LOWER}-${RANGE_UPPER}.tsv' \
							'$MODE' \
							'$ALPHA' \
							'$CIS_DIST' \
							'$SLICE_SIZE' \
							'$RANGE_LOWER' \
							'$RANGE_UPPER' \
							'$FORMAT_CHECK'"


				case "$USE_CLUSTER" in

					false )

						eval "$cmd_main"
						;;


					true )

						njobs=$(( njobs + 1 ))
						echo $cmd_main >> "$jobs_file"
						;;

				esac


				RANGE_LOWER=$(( RANGE_LOWER + multiple ))
				RANGE_UPPER=$(( RANGE_UPPER + multiple ))

				if [[ "$RANGE_LOWER" -gt "$n_tests" ]]; then

					break

				fi

				if [[ "$RANGE_UPPER" -gt "$n_tests" ]]; then

					RANGE_UPPER=$n_tests

				fi


			done

			;;


		peak-tests )

			SLICE_SIZE=2 # ~2 snps per peak
			FORMAT_CHECK="FALSE"

			
			## Check input format

			format_check_meqtl_io \
				"$in_dir/genotype_NA_source.tsv" \
				"$in_dir/phenotype_source.tsv" \
				"$in_dir/covariates_source.tsv"


			## Iterate over user defined multiples of tests
			
			multiple=$MULTIPLE
			n_tests=$(cut -f8 "$tests_bed" | sort | uniq | wc -l)
			
			RANGE_LOWER=1
			RANGE_UPPER=$multiple

			while true; do

				cmd_main="Rscript --verbose 'code/call_qtls_matrix-eqtl.R' \
						'$in_dir' \
						'$out_dir' \
						'tests_peak-level_${RANGE_LOWER}-${RANGE_UPPER}.tsv' \
						'$MODE' \
						'$ALPHA' \
						'$CIS_DIST' \
						'$SLICE_SIZE' \
						'$RANGE_LOWER' \
						'$RANGE_UPPER' \
						'$FORMAT_CHECK'"


				case "$USE_CLUSTER" in

					false )

						eval "$cmd_main"
						;;


					true )

						njobs=$(( njobs + 1 ))
						echo $cmd_main >> "$jobs_file"
						;;

				esac


				RANGE_LOWER=$(( RANGE_LOWER + multiple ))
				RANGE_UPPER=$(( RANGE_UPPER + multiple ))

				if [[ "$RANGE_LOWER" -gt "$n_tests" ]]; then

					break

				fi

				if [[ "$RANGE_UPPER" -gt "$n_tests" ]]; then

					RANGE_UPPER=$n_tests

				fi


			done

			;;

	esac


	## Job array


	if [[ "$USE_CLUSTER" == "true" ]]; then

			echo -e "\t- n_jobs = ${njobs}"

			bsub << EOF
#!/usr/bin/env bash
#BSUB -R "rusage[mem=2G]"
#BSUB -q medium
#BSUB -cwd ${PROJECT_DIR}
#BSUB -J "${job_id}[1-${njobs}]"
#BSUB -o ${PROJECT_DIR}/code/bsub/logs/${job_id}.%I.out
#BSUB -e ${PROJECT_DIR}/code/bsub/logs/${job_id}.%I.err

set -euo pipefail

cd $PROJECT_DIR

source "$HOME/.bash_profile"
load-r-440


## Run
sed -n "\${LSB_JOBINDEX}p" ${jobs_file} | bash
EOF

	fi
    
done < <(find "${MATRIX_EQTL_INPUT_DIR}" -mindepth 1 -maxdepth 1 -type d -print0)
