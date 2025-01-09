#!/usr/bin/env bash
#
# Call caQTLs using matrix-eQTL on all cell-types
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
MODE="bulk-tests"
USE_CLUSTER=false

while getopts "m:c" opt; do

	case "$opt" in

		#m ) MODE=$OPTARG ;;
		c ) USE_CLUSTER=true ;;
		? ) echo "Usage: $0 [-c] [-m mode {bulk-tests, pea-tests, single-tests}"; exit 1 ;;
	
	esac

done

case "$MODE" in

	single-tests | bulk-tests | peak-tests )

		:
		;;

	unset )

		echo "Usage: $0 [-c] [-m mode {bulk-tests, peak-tests, single-tests}"
		exit 1
		;;

	? )

		echo "Usage: $0 [-c] [-m mode {bulk-tests, peak-tests, single-tests}"
		exit 1
		;;

esac


## Environment

PROJECT_DIR=".."
PROJECT_DIR="$(realpath ${PROJECT_DIR})"
cd $PROJECT_DIR

source code/glob_vars.bash # MATRIX_EQTL_INPUT_DIR, MATRIX_EQTL_OUTPUT_DIR

if [[ "$USE_CLUSTER" == "false" ]]; then

	source "${HOME}/.bash_profile"
	load-r-440

fi


## User variables

ALPHA=0.01 # Aygun-2023 caQTLS: 0.05
CIS_DIST=1000 # Aygun-2023 caQTLS: +-1 kbp
MULTIPLE=1 # Nr of test-batches to perform for every job
			 # For bulk-tests: All tests in 1 tests-batch
			 # For single-tests: 1 test in 1 tests-batch
			 # For peak-tests: ~2 tests in tests-batch


### SCRIPT ###

while IFS= read -r -d '' cell_type_path; do

	cell_type="$(basename "$cell_type_path")"

	echo -e "$cell_type"


	echo -e "\t- Multiples: ${MULTIPLE}"

	in_dir="${cell_type_path}/${MODE}"
	out_dir="${MATRIX_EQTL_OUTPUT_DIR}/chromatin-accessibility/${CT_MAP_ID}/${cell_type}/${MODE}"
	# tests_bed="${cell_type_path}/${MODE}/tests_snp_peak_pairs.bed"
	job_id="call_caqtls_$(date -I)_${cell_type}"


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
		jobs_file="code/bsub/logs/call_caqtls_$(date -I)_${cell_type}.commands"

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

			echo "Not implemented yet"
			exit 1
			SLICE_SIZE=1 # Only 1 test at a time
			FORMAT_CHECK="FALSE"


			## Check input format

			format_check_meqtl_io \
				"$in_dir/genotype_NA_source.tsv" \
				"$in_dir/phenotype_source.tsv" \
				"$in_dir/covariates_source.tsv"


			## Iterate over user defined multiples of tests
			
			n_tests=$(wc -l "$tests_bed")
			multiple="$MULTIPLE"

			RANGE_LOWER=1

			if [[ "$n_tests" -ge "$MULTIPLE" ]]; then

				RANGE_UPPER=$multiple

			else

				RANGE_UPPER="$n_batches"

			fi
			

			while true; do

				cmd_main="Rscript --verbose 'code/call_qtls_matrix-eqtl.R' \
							'$in_dir' \
							'$out_dir' \
							'tests_snp-level_${RANGE_LOWER}-${RANGE_UPPER}.tsv' \
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

			echo "Not implemented yet"
			exit 1

			SLICE_SIZE=2 # ~2 snps per peak
			FORMAT_CHECK="FALSE"

			
			## Check input format

			format_check_meqtl_io \
				"$in_dir/genotype_NA_source.tsv" \
				"$in_dir/phenotype_source.tsv" \
				"$in_dir/covariates_source.tsv"


			## Iterate over user defined multiples of test batches
			
			n_batches=$(cut -f8 "$tests_bed" | sort | uniq | wc -l)
			multiple="$MULTIPLE"

			RANGE_LOWER=1

			if [[ "$n_batches" -gt "$multiple" ]]; then

				RANGE_UPPER=$multiple

			else

				RANGE_UPPER="$n_batches"

			fi
			

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

				if [[ "$RANGE_LOWER" -gt "$n_batches" ]]; then

					break

				fi

				if [[ "$RANGE_UPPER" -gt "$n_batches" ]]; then

					RANGE_UPPER=$n_batches

				fi


			done

			;;

	esac


	## Job array


	if [[ "$USE_CLUSTER" == "true" ]]; then

			echo -e "t\t\t- n_jobs = ${njobs}"

			bsub << EOF
#!/usr/bin/env bash
#BSUB -R "rusage[mem=20G]"
#BSUB -q long
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

done < <(find "${MATRIX_EQTL_INPUT_DIR}/chromatin-accessibility/${CT_MAP_ID}" -mindepth 1 -maxdepth 1 -type d -print0)
