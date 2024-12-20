format_check_meqtl_io() {
	# Check that all tsv files' headers contain the same sample names in equal order and uniqueness

	## Gather sample headers in array
	
	local tsvs=("$@")
	local sample_headers=()

	for tsv in "${tsvs[@]}"; do

		local sample_header="$(head -n 1 "$tsv" | cut -f2-)"
		sample_headers+=("$sample_header")


		## Check for sample uniqueness
		
		local n_samples_uniq="$(echo "$sample_header" | tr '\t' '\n' | sort | uniq | wc -l)"
		local n_samples="$(echo "$sample_header" | tr '\t' '\n' | wc -l)"

		if [[ "$n_samples_uniq" -ne "$n_samples" ]]; then

			echo "File <$tsv> has repeated header samples."
			return 1

		fi

	done


	## Check for equal samples and order

	local first_sample_header="${sample_headers[0]}"
	local counter=1

	for sample_header in "${sample_headers[@]}"; do


		if [[ "$sample_header" != "$first_sample_header" ]]; then

			echo "Samples in headers are not identical in all files."
			echo "Header of file 1 is different from file ${counter}"
			return 1

		fi

		local counter=$((counter + 1))

	done


	return 0 # If all correct

}
