frag_file_to_bw_chrombpnet() {
	# Function: create_bw_from_frag_file
	# Description: Create bigwig coverage file.
	# Arguments:
	#	$1 - str - Path to ChromBPNet `reads_to_bigwig.py` script
	#	$2 - str - Path to fragment file location which will be used to create the profile.
	#		   .tsv or .tsv.gz
	#		   ```ChromBPNet: Both gzipped and plain text files are allowed. Sorting is not expected.```
	#	$3 - str - Path to the desired output bigwig file.
	#			   Do not include file extension.
	#	$4 - str - Path to reference genome fasta file file.
	#	$5 - str - Path to chromosoe sizes file.
	#			   Line format: "chr1\t10000000".
	#			   ! Chromosome naming must probably have the same format as in the frag file.
	#	$6 - str - Assay type ["ATAC" | "DNASE"].
	#	$7 - int - Plus shift of assay type. In case ChromBPNet doesn't automatically detect it.
	#				   Reevaluate with *.
	#				   Can be cleverly exploited to introduce custom shift.
	#				   E.g. og=+4 --> cust=+5 via plus_shift=+3 (+3-4=+1)
	#	$8 int - Minus shift of assay type. In case ChromBPNet doesn't automatically detect it.
	#				   Reevaluate with *.
	#				   Can be cleverly exploited to introduce custom shift.
	#				   E.g. og=-5 --> cust=-3 via minus_shift=-5 (-(-5)+(-4)=+1)
	#	$9 - str - opt - Use cluster? [true | false] (Default: false)
	#	$10 -str - opt - Project directory, needed if $9 == true
	# Returns:
	#	Bigwig file with coverage obtained from fragment file.
	#	Underlying base-pair indexing based on the indexing of the input fragment file.
	#	Usually 0-based half-open.
	#	Under the hood, ChromBPNet `reads_to_bigwig.py` uses a bedgraph conversion file.
	#	However, it still reads the fragment insertions from the fragment files so it's dependent on these.
	# Example:
	#		create_bw_from_frag_file \
	#			"~/.conda/envs/ian/lib/python3.9/site-packages/chrombpnet/helpers/preprocessing/reads_to_bigwig.py" \
	#			"$frag_files/frag-file-X.tsv.gz" \
	#			"grouped-frag-files/frag-file-X-coverages" \
	#			"data/GRCh38-p14/hg38.fa" \
	#			"ATAC"
	#
	# Notes:
	# 	- New*: ChromBPNet automatically detects the Tn5-shift and shifts it to +4/-4. If given manually, it will shift the fragments by the difference between the input shift and +4/-4. BE CAREFUL. E.g. manual args: ps=+4, ms=-5 => frag_start +0, frag_end +1; ps=0, ms=0 => first insertion +4, second insertion -5

	source "code/helpers/bash/utils.bash"


	## Args
	
	if [[ $# -lt 8 ]] || [[ $# -gt 10 ]]; then

		echo "Incorrect nr of arguments: <$#>"
		return 1

	fi

	local cbn_py_file="$1"
	local frag_file="$2"
	local out_bw="$3"
	local ref_genome_file="$4"
	local chrom_sizes="$5"
	local assay_type="$6"
	local plus_shift="$7"
	local minus_shift="$8"
	local use_cluster="${9:-false}"
	local project_dir="${10:-""}"

	local file_base
	file_base="$(basename "${frag_file%%.*}")"

	local out_dir="$(dirname "$out_bw")"
	local out_base="$(basename "$out_bw")"

	local randn="$(shuf -i 100000-999999 -n1)"
	local frag_file_tmp="${out_dir}/tmp/${out_base}.${randn}.gz"
	mkdir -p "${out_dir}/tmp"


	## Checks
	
	if [[ ! -f "$frag_file" ]]; then
	
		echo "Fragment file doesn't exist: ${frag_file}"
		return 1

	fi

	if [[ ! -f "$cbn_py_file" ]]; then

		echo "ChromBPNet python script not found! <$cbn_py_file> doesn't exist."
		return 1

	fi


	## Run
	
	case $use_cluster in

		false )

			strip_and_sort_fragment_file "$frag_file" "$frag_file_tmp"

			## Handle empty files: ignore
			
			case "$frag_file" in

				*.tsv )

					if [[ ! -s "$frag_file_tmp" ]]; then

						#touch "${out_bw}.bw"
						#rm "$frag_file_tmp"
						return 0

					fi
					;;

				*.gz )

					if [[ -z "$(bgzip -dc "$frag_file_tmp" | head -c 1 | tr '\0\n' __)" ]]; then

						#touch "${out_bw}.bw"
						#rm "$frag_file_tmp"
						return 0

					fi
					;;

				* )

					echo "Input file [.tsv | .tsv.gz] not <${frag_file_tmp}>"
					return 1
					;;

			esac


			## Build main command

			cmd_main="python \"$cbn_py_file\" \
				-ifrag \"$frag_file_tmp\" \
				-o \"$out_bw\" \
				-c \"$chrom_sizes\" \
				-g \"$ref_genome_file\" \
				-d \"$assay_type\""

			
			if [[ -n "$plus_shift" ]] && [[ -n "$minus_shift" ]]; then

				cmd_main="$cmd_main \
					-p $plus_shift \
					-m $minus_shift"


			elif [[ -n "$plus_shift" ]] || [[ -n "$minus_shift" ]]; then

				echo "If you are going to indicate the shifts, it must be both or none: ps=<$plus_shift> and ms=<$minus_shift>"
				return 1

			fi
			

			## Run
			
			eval "$cmd_main"
			mv "${out_bw}_unstranded.bw" "${out_bw}.bw"
			#rm "$frag_file_tmp"
			;;


		true )

			if [[ -z "$project_dir" ]]; then

				echo "When using the cluster indicate the project path! Argument 4"
				return 1

			fi

			local job_id="compute_coverage_$(date '+%Y-%m-%d')_${file_base}_low_mem"

			continue
			# Legacy code: delete if no problems
			cat <<EOF
#!/usr/bin/env bash
#BSUB -R "rusage[mem=10G]"
#BSUB -q long
#BSUB -cwd ${project_dir}
#BSUB -J ${job_id}
#BSUB -o ${project_dir}/code/bsub/logs/${job_id}.out
#BSUB -e ${project_dir}/code/bsub/logs/${job_id}.err

set -euo pipefail

cd ${project_dir}
source ~/.bash_profile
source "code/helpers/bash/utils.bash"

# python3.9
cmd_main="python \"$cbn_py_file\" \
	-ifrag \"$frag_file_tmp\" \
	-op \"$out_bw\" \
	-c \"$chrom_sizes\" \
	-g \"$ref_genome_file\" \
	-d \"$assay_type\""

strip_and_order_fragment_file "$frag_file" "$frag_file_tmp"

# Check if file is empty and return empty if so
case "$frag_file" in

	*.tsv )

		if [[ ! -s "$frag_file_tmp" ]]; then

			touch "${out_bw}.bw"
			rm "$frag_file_tmp"
			exit 0

		fi
		;;

	*.gz )

		if [[ -z "\$(bgzip -dc "$frag_file_tmp" | head -c 1 | tr '\0\n' __)" ]]; then

			touch "${out_bw}.bw"
			rm "$frag_file_tmp"
			exit 0

		fi
		;;

	* )
		echo "Input file [.tsv | .tsv.gz] not <${frag_file_tmp}>"
		exit 1
		;;

esac

if [[ -n "$plus_shift" ]] && [[ -n "$minus_shift" ]]; then

	cmd_main="\$cmd_main \
		-p $plus_shift \
		-m $minus_shift"

elif [[ -n "$plus_shift" ]] || [[ -n "$minus_shift" ]]; then

	echo "If you are going to indicate the shifts, it must be both not: ps=<$plus_shift> and ms=<$minus_shift>"
	exit 1
	
fi
			
eval "\$cmd_main"
mv "${out_bw}_unstranded.bw" "${out_bw}.bw"
rm "$frag_file_tmp"
EOF
			;;

		* )
			echo "Argument 9: use_cluster [true | false] not <$use_cluster>"
			return 1
			;;

	esac

}
