#!/usr/bin/env bash
#
# Call footprintQTLs using matrix-eQTL on all cell-types
#
# Requires: 
#	- phenotype_file: PEAKS x DONORS
#	- peak_locations_file: PEAKS x (chr, start, end)
#	- phenotype_pcs_file: PCS x DONORS
#	- genotype_file: SNPS x DONORS
#	- snp_locations_file: SNPS x (chr, pos)
#	- genotype_pcs_file: PCS x DONORS
#	-> all tsvs
#	-> for details see:
#		https://github.com/andreyshabalin/MatrixEQTL/tree/master/data
#
# 
# Assumptions:
# 	- Correctly setup renv, renv.lock, .Rprofile, .Renviron
# 	- In input, the order of DONORS columns is the same between files.
# 	- Input directory that is scanned has only directories with the different cell-types/categories
# 	  which contain the matrix-eQTL data for those categories.


### SETUP ###

set -eou pipefail
set -x


# Args

MODE="unset"
USE_CLUSTER=false

while getopts "m:c" opt; do

	case "$opt" in

		m ) MODE=$OPTARG ;;
		c ) USE_CLUSTER=true ;;
		? ) echo "Usage: $0 [-c] [-m mode {bulk-tests,single-tests}"; exit 1 ;;
	
	esac

done

if [[ "$MODE" == "unset" ]]; then

	echo "Usage: $0 [-c] [-m mode {bulk-tests,single-tests}"
	exit 1

fi

# Environment

PROJECT_DIR=".."
PROJECT_DIR="$(realpath ${PROJECT_DIR})"
cd $PROJECT_DIR

source code/glob_vars.bash # MATRIX_EQTL_INPUT_DIR, MATRIX_EQTL_OUTPUT_DIR

if [[ "$USE_CLUSTER" == "false" ]]; then

	source "${HOME}/.bash_profile"
	load-r-440

	njobs=0
	jobs_file="code/bsub/logs/call_fqtls_$(date -I).commands"
	
fi


### SCRIPT ###

job_id="call_fqtls_$(date -I)"


# Run function for all cell-types

while IFS= read -r -d '' cell_type_dir; do

	echo "Processing file: ${cell_type_dir}."

	cell_type="$(basename "$cell_type_dir")"
	in_dir="${cell_type_dir}/${MODE}"
	out_dir="${MATRIX_EQTL_OUTPUT_DIR}/${cell_type}/${MODE}"


	# Clean old run files and create dir
	
	if [[ -d "${out_dir}" ]]; then

		rm -rf "${out_dir}" 

	fi

	mkdir -p "$out_dir"


	case "$MODE" in

		bulk-tests )

			phenotype_file="${in_dir}/footprints.tsv"
			peak_locations_file="${in_dir}/peak_locations.tsv"
			genotype_file="${in_dir}/genotype_NA.tsv"
			snp_locations_file="${in_dir}/snp_locations.tsv"
			covariates_file="${in_dir}/covariates.tsv"

			cmd="Rscript --verbose 'code/call_footprint-qtls_matrix-eqtl.R' \
				'$phenotype_file' \
				'$peak_locations_file' \
				'$phenotype_pcs_file' \
				'$genotype_file' \
				'$snp_locations_file' \
				'$genotype_pcs_file' \
				'$out_dir' \
				'$job_id'"

			case "$USE_CLUSTER" in

				false )
					eval "$cmd"
					;;

				true )
					njobs=$(( $njobs + 1 ))
					echo $cmd >> "$jobs_file"
					;;

			esac

			;;

		single-tests )

			while IFS='\t' read -r pair; do

				snp="$(echo "$pair" | cut -f 4)"
				peak="$(echo "$pair" | cut -f 8)"

				phenotype_file="${in_dir}/phenotypes/footprints%${peak}.tsv"
				peak_locations_file="${in_dir}/peak_locations/peak_location%${peak}.tsv"
				genotype_file="${in_dir}/genotypes/genotype_NA%${snp}.tsv"
				snp_locations_file="${in_dir}/snp_locations/snp_location%${snp}.tsv"
				covariates_file="${in_dir}/covariates/covariates%${peak}.tsv"

				cmd="Rscript --verbose 'code/call_footprint-qtls_matrix-eqtl.R' \
					'$phenotype_file' \
					'$peak_locations_file' \
					'$genotype_file' \
					'$snp_locations_file' \
					'$covariates_file' \
					'$out_dir' \
					'$job_id' \
					'FALSE' \
					'1'"

				case "$USE_CLUSTER" in

					false )
						eval "$cmd"
						;;

					true )
						njobs=$(( $njobs + 1 ))
						echo $cmd >> "$jobs_file"
						;;

				esac

			done < "${MATRIX_EQTL_INPUT_DIR}/${cell_type}/single-tests/tests_snp_peak_pairs.bed"

			;;

	esac
    
done < <(find "${MATRIX_EQTL_INPUT_DIR}" -mindepth 1 -maxdepth 1 -type d -print0)




# Job array

if [[ "$USE_CLUSTER" == "true" ]]; then

	bsub << EOF
#!/usr/bin/env bash
#BSUB -R "rusage[mem=50G]"
#BSUB -q medium
#BSUB -cwd ${PROJECT_DIR}
#BSUB -J "${job_id}[1-$njobs]"
#BSUB -o ${PROJECT_DIR}/code/bsub/logs/${job_id}.%I.out
#BSUB -e ${PROJECT_DIR}/code/bsub/logs/${job_id}.%I.err

set -euo pipefail

cd $PROJECT_DIR

source "$HOME/.bash_profile"
load-r-440

sed -n "\${LSB_JOBINDEX}p" ${jobs_file} | bash
eval "$cmd"
EOF

fi
