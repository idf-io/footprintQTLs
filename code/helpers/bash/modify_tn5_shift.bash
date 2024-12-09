modify_tn5_shift() {

	# Change the shift of the Tn5 insertion sites of ATAC-seq data assuming +4/-5
	# Finally compress and index
	# Input: ATAC-seq fragments file, +4/-5 shift assumed , compressed with bgzip
	# Output:
	#  1. DATA: File .tsv.gz or directory of .tsv.gz files to modify
	#  2. OUT_DIR
	#  3. plus shift: new plus shift e.g. 1 would do start -3 since 

	DATA="$1"
	OUT_DIR="$2"
	plus_shift_og="$3"
	minus_shift_og="$4"
	
	plus_shift=$(( plus_shift_og - 4))
	minus_shift=$(( minus_shift_og + 5))

	if [ -d "$DATA" ]; then

		# Create shifted files
		find "$DATA" -maxdepth 1 -mindepth 1 -type f -name "*.tsv.gz" | while IFS= read -r f; do
			
			bgzip -dc "$f" | grep -v "^#" | awk -v ps="$plus_shift" -v ms="$minus_shift" -F'\t' 'BEGIN { OFS="\t" } { $2=$2+ps; $3=$3+ms; print }' - > "${OUT_DIR}/$(basename "${f%.tsv.gz}")_ps${plus_shift_og#+}_ms${minus_shift_og#-}.tsv"
		
		done
	
		# Compress and index
		if [ -n "${5-}" ]; then
	
			bgzip_and_index_tsvs "$OUT_DIR"
	
		fi

	elif [ -f "$DATA" ]; then

		# Create shifted file
		bgzip -dc "$DATA" | grep -v "^#" | awk -v ps="$plus_shift" -v ms="$minus_shift" -F'\t' 'BEGIN { OFS="\t" } { $2=$2+ps; $3=$3+ms; print }' - > "${OUT_DIR}/$(basename "${DATA%.tsv.gz}")_ps${plus_shift_og#+}_ms${minus_shift_og#-}.tsv"

	else

		echo "First argument is neither a file nor a directory"
		exit 1

	fi
	

}
