bgzip_and_index_tsv() {
	# Function: bgzip_and_index_tsv
	# Description: Block compresses and indexes the input tsv
	# Arguments:
	# 	1 - str - tsv file path
	# 	2 opt - str - Keep the original tsv file? [Default: False | True]
	# 	3 opt - str - New directory to move the compressed and tbi files. Default: cwd
	# Output:
	# 	1 - Compressed tsv file `.tsv.gz`
	# 	2 - Cognate indexing file `.tsv.gz.tbi`
	# Example:
	# 	bgzip_and_index_tsv $HOME/documents/fragments_file.tsv True

	local tsv="$1"
	local keep_tsv="${2:-"False"}"
	local out_dir="${3:-$(dirname "$1")}"
	
	local tsv_basename="$(basename $tsv)"

	mkdir -p $out_dir

	# Catch user errors
	if [[ "$#" -eq 0 || "$#" -gt 3 ]]; then

		echo "Too many args. 1: tsv file path, 2: Keep the original tsv? [def: False | True], 3 opt: new folder"
		exit 1

	fi


	# Compress tsv
	({ grep "^#" $tsv || true; }; { grep -v "^#" $tsv || true; } | sort -t"`printf '\t'`" -k1,1 -k2,2n) | bgzip -c > "${out_dir}/${tsv_basename}.gz"
    	tabix -p bed "${out_dir}/${tsv_basename}.gz"


	# Optionally remove original tsv
	if [[ $keep_tsv == "False" ]]; then

		rm $tsv

	elif [[ $keep_tsv == "True" ]]; then
		
		:

	else

		echo "Wrong argument nr 2: Keep false tsv? [False | True]"
		exit 1

	fi

}


check_jobs() {
	# Function: Given an array of job_ids, check if they are all finished.
	# Arguments:
	# 	job_ids - array - Job ids (str)
	# Output:
	# 	0 - If all jobs are absent and not in RUN or PEND mode.
	# 	1 - If any job is in RUN or PEND mode.
	# Example:
	# 	a=("0000001" "00000002" "00000003")
	# 	if check_jobs "${a[@]}"; then
	# 		echo "All jobs finished."
	# 	else
	# 		echo "Jobs not finished."
	# 	fi
	# Note:

	local cluster_jobs=("$@")
	local job # Ensure that there is no jobs variable from higher scope

	for job in "${cluster_jobs[@]}"; do

		if bjobs "$job" 2>/dev/null | grep -E "RUN|PEND|UNKWN" > /dev/null; then

			return 1

		fi

	done

	return 0

}

strip_and_sort_fragment_file() {
	# Function: Remove comment lines from fragment file and order it after
	# 				1) chromosome and 2) start position, then compress
	# Arguments:
	# 	1 - file_path - str - Handles zipped files automatically.
	# 	2 - new_file_path
	# Example:
	# 	strip_file fragments.tsv fragments_stripped.tsv
	# Notes:
	# 	- Compression done with bgzip
	# Assumptions:
	# 	Text file
	# 	Uniquely comment lines start with #
	# TODO:
	# 	- make compressing optional
	# 	- Include arg to exit 1 with empty file instead of creating an empty compressed file
	
	if [[ "$1" == *.tsv ]]; then

		{ grep -v "^#" "$1"  || true; } | sort -t"`printf '\t'`" -k1,1 -k2,2n | bgzip > "$2"

	elif [[ "$1" == *.tsv.gz ]]; then
		
		bgzip -dc $1 | { grep -v "^#" || true; } | sort -t"`printf '\t'`" -k1,1 -k2,2n | bgzip > "$2"

	else
		
			echo "Fragment file must be [.tsv | .tsv.gz]"
			exit 1

	fi

}



# Command line utilities

cmp_two_dirs() {
	# TLDR: Compare the files of two directories for equality.
	# Arguments:
	# 	Directory 1
	# 	Directory 2
	# Output:
	# 	List of files that differ
	# Example:
	# 	cmp_two_dirs dir1 dir2
	# Assumptions
	
	DIR1=$1
	DIR2=$2

	# Check if both directories exist
	if [ ! -d "$DIR1" ] || [ ! -d "$DIR2" ]; then

  		echo "One or both directories do not exist."
  		exit 1

	fi

	# Use diff to compare directories recursively, and cmp for binary files
	diff -qr "$DIR1" "$DIR2" | while read -r line; do

  		# If the output suggests the files differ, we compare using cmp
  		if [[ $line == *"differ"* ]]; then

    			FILE1=$(echo $line | cut -d ' ' -f 2)
    			FILE2=$(echo $line | cut -d ' ' -f 4)
    			cmp "$FILE1" "$FILE2" && echo "Binary files are identical." || echo "Binary files are different."

  		fi

	done
	
}






# Anndata

extract_homogeneous_col_from_anndata_obs() {

	local filename="$1"
	local col="$2"

	value="$(h5dump -d "/obs/${col}/categories" "$filename" |
				awk '/^   DATA {/,/}/' |
				grep -e '   ([0-9]*):' |  
				awk -F '"' '{ for(i=2; i<=NF; i+=2) print $i}')"

	echo "$value"

}
