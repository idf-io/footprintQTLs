split_frag_file_by_groups() {
	# Function: Split fragments of a fragment file by specified groups.
        # Arguments:
        #       1 - str - Fragment file in tsv and bed format.
        #                 File name without extension will be considered the sample name.
        #       2 - str - Output directory
        #       3 - str - Path to files named after a grouping
        #                 and containing the new-line separated list of cells pertaining to that group.
        #                 File example: donorX_cell-typeY.txt
	#       4 - str, opt - Compress output? [True | False] (Defaults to False)
        #       5 - str, opt - Use cluster? [True | False] (Defaults to False)
	#       6 - str, opt - If using cluster indicate the project path. (Defaults to False)
	# Output:
	# 	The output directory is populated with the group-level fragment files.
	# 	In an intermediate step, output_directory/tmp is populated with group-sample level fragments.
	# Example:
	# 	split_frag_file_by_groups \
	# 		$HOME/docs/myProject/sample01.tsv.gz \
	# 		$HOME/docs/myProject/output \
	# 		$HOME/docs/myProject/cell-id_group_files \
	# 		True
	# Assumptions:
	# 	Fragment files:
	# 		fragment ids in format barcode_sample
	# 		named by their sample (case sensitive). Important because cell-ids are formatted using the sample: barcode_id
	# 	h5ad:
	# 		file with group/cell-type annotations in obs has the observations inde
	# Note:
	# 	Currently it's implemented such that group=donor_cell-type.
	#	Function replaces older files if existent.
	#	Function populates a group folder per frag_file/sample. If you execute the function for several samples with the same grouping scheme they will populate the same folders. Don't delete these folders by accident thinking they only contain the output of one sample.

	source "code/helpers/fragment_file_helpers.bash"
	
	local frag_file="$1"
	local out_dir="$2"
	local cellids_group_files_dir="$3"
	local compress="${4:-False}"
	local use_cluster="${5-False}"
	local project_dir="${6:-""}"

	local sample_id="$(basename $frag_file)"

	case "$sample_id" in 

		*.tsv.gz) local sample_id_short="${sample_id%.tsv.gz}" ;;
		*.tsv) local sample_id_short="${sample_id%.tsv}" ;;
		*)
			echo "File must be a tsv file! Optionally gzipped"
			return 1
			;;

	esac

	local date=$(date '+%Y-%m-%d')

	
	## Setup

	# Check if `cellids_group_files_dir` empty
	if ls "$cellids_group_files_dir"/* 1> /dev/null 2>&1; then

		:

	else

		echo "Directory <${cellids_group_files_dir}> is empty!"
		return 1

	fi


	## Main

	case $use_cluster in

		"False")

			shopt -s nullglob
			local cellids_group_files=("$cellids_group_files_dir"/*)
			#local tmp_path="/home/fichtner/projects/footprintQTL/data/intermediate-data/datasets/hca_brain-organoids/annotations_general/group_cell-ids/approach_2024-09-12"
			#local cellids_group_files=(
        			#"${tmp_path}/pelm_Glia.txt"
        			#"${tmp_path}/pelm_UL-EN.txt"
        			#"${tmp_path}/zoxy_Glia.txt"
			#)
		
			# For each cellids group file
			for cgf in "${cellids_group_files[@]}"; do 

        			local group="$(basename ${cgf%.*})"
				local out_dir_tmp="${out_dir}/tmp/${group}"

				mkdir -p $out_dir_tmp

		
				# Construct awk command parts
		
				# If input frag_file compressed
				case $frag_file in
		
					*.tsv.gz) local decompress_input_cmd="<(bgzip -dc \"${frag_file}\")" ;;
					*.tsv) local decompress_input_cmd="${frag_file}" ;;
					*)
						echo -e "Argument 1: fragment_file_path [*.tsv | *.tsv.gz (block zipped)] not <${frag_file}>."
						return 1
						;;
		
				esac
		
				# If output should be compressed
				case $compress in
		
					"True")
						local compress_cmd="| gzip"
						out_tmp_file="${out_dir_tmp}/${sample_id_short}%${group}.tsv.gz"
						;;

					"False")
						local compress_cmd=""
						out_tmp_file="${out_dir_tmp}/${sample_id_short}%${group}.tsv"
						;;

					*)
						echo "Argument 4: compress? [True | False] not <${compress}>."
						return 1
						;;
		
				esac


				# If file exists: delete
				if [[ -f "$out_tmp_file" ]]; then

					rm  "$out_tmp_file"

				fi

		
				# Main awk command
				local main_cmd="awk -v sample=\"$sample_id_short\" -F'\t' 'NR==FNR { cell_ids[\$1] ; next } { if ((\$4 \"_\" sample) in cell_ids) print \$0 }' OFS='\t' \"$cgf\" $decompress_input_cmd $compress_cmd > $out_tmp_file"
				# Execute main command
				#echo "$main_cmd"
				eval "$main_cmd"

				
			
			done
			;;
	
	
		"True")

			if [[ -z $project_dir ]]; then

				echo "When using the cluster indicate the project path! Argument 6"
				return 1

			fi

			local job_id="split_frags_bgroups_${date}_${sample_id_short}"

			#BSUB -W ${walltime}     # Walltime in hours:minutes
			
			cat <<EOF
#!/usr/bin/env bash
#BSUB -R "rusage[mem=50G]"
#BSUB -q verylong
#BSUB -cwd ${project_dir}
#BSUB -J ${job_id}
#BSUB -o ${project_dir}/code/bsub/logs/${job_id}.out
#BSUB -e ${project_dir}/code/bsub/logs/${job_id}.err

set -euo pipefail

cd ${project_dir}
source ~/.bash_profile
load-conda
conda activate ian

#CONDA_ENV_NAME="your_env_name"
#CONDA_PATH="\$HOME/miniconda3"  # Adjust this path as needed
#source "\${CONDA_PATH}/etc/profile.d/conda.sh"
#conda activate \$CONDA_ENV_NAME

# Main code to split the fragment file by groups
shopt -s nullglob
cellids_group_files=("$cellids_group_files_dir"/*)

tmp_path="/home/fichtner/projects/footprintQTL/data/intermediate-data/datasets/hca_brain-organoids/annotations_general/group_cell-ids/approach_2024-09-12"

#cellids_group_files=(
        #"\${tmp_path}/pelm_Glia.txt"
        #"\${tmp_path}/pelm_UL-EN.txt"
        #"\${tmp_path}/zoxy_Glia.txt"
#)


# For each cellids group file
for cgf in "\${cellids_group_files[@]}"; do 

	group="\$(basename \${cgf%.*})"
	out_dir_tmp="${out_dir}/tmp/\${group}"

	mkdir -p \$out_dir_tmp

	# Construct awk command parts
	# If input frag_file compressed
	case $frag_file in

		*.tsv.gz) decompress_input_cmd="<(bgzip -dc \"${frag_file}\")" ;;
		*.tsv) decompress_input_cmd="${frag_file}" ;;
		*)
			echo -e "Argument 1: fragment_file_path [*.tsv | *.tsv.gz (block zipped)] not <${frag_file}>."
			exit 1
			;;
	esac

	# If output should be compressed
	case $compress in
		"True")
			compress_cmd="| gzip"
			out_tmp_file="\${out_dir_tmp}/${sample_id_short}%\${group}.tsv.gz"
			;;
		"False")
			compress_cmd=""
			out_tmp_file="\${out_dir_tmp}/${sample_id_short}%\${group}.tsv"
			;;
		*)
			echo "Argument 4: compress? [True | False] not <${compress}>."
			exit 1
			;;
	esac


	# If file exists: delete
	if [[ -f "\$out_tmp_file" ]]; then

		rm  "\$out_tmp_file"
	
	fi

		

	# Main awk command
	main_cmd="awk -v sample=\"${sample_id_short}\" -F'\t' 'NR==FNR { cell_ids[\\\$1] ; next } { if ((\\\$4 \"_\" sample) in cell_ids) print \\\$0 }' OFS='\t' \"\$cgf\" \$decompress_input_cmd \$compress_cmd > \$out_tmp_file"
	
	# Execute main command
	eval "\$main_cmd"

done

EOF
			;;


		*)

			echo "Argument #4 [True, False] not <${use_cluster}>"
			return 1
			;;

	esac

}
