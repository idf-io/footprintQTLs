#!/usr/bin/env bash
set -euo pipefail

create_bw_from_frag_file_cbn() {
	# Function: create_bw_from_frag_file
	# Description: Create bigwig coverage file.
	# Arguments:
	#	$1 - str - Path to ChromBPNet `reads_to_bigwig.py` script
	#	$2 - str - Path to fragment file location which will be used to create the profile.
	#			   ```ChromBPNet: Both gzipped and plain text files are allowed. Sorting is not expected.```
	#	$3 - str - Path to the desired output bigwig file.
	#			   Do not include file extension.
	#	$4 - str - Path to reference genome fasta file file.
	#	$5 - str - Path to chromosoe sizes file.
	#			   Line format: "chr1\t10000000".
	#			   ! Chromosome naming must probably have the same format as in the frag file.
	#	$6 - str - Assay type ["ATAC" | "DNASE"].
	#	$7 opt - int - Plus shift of assay type. In case ChromBPNet doesn't automatically detect it.
	#				   Can be cleverly exploited to introduce custom shift.
	#				   E.g. og=+4 --> cust=+5 via plus_shift=+3 (+3-4=+1)
	#	$8 opt - int - Minus shift of assay type. In case ChromBPNet doesn't automatically detect it.
	#				   Can be cleverly exploited to introduce custom shift.
	#				   E.g. og=-5 --> cust=-3 via minus_shift=-5 (-(-5)+(-4)=+1)
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

	local cbn_py_file="$1"
	local frag_file="$2"
	local out_bw="$3"
	local ref_genome_file="$4"
	local chrom_sizes="$5"
	local assay_type="$6"

	if [[ $# -eq 6 ]]; then

		# python3.9
		python "$1" \
			-ifrag "$frag_file" \
			-op "$out_bw" \
			-c "$chrom_sizes" \
			-g "$ref_genome_file" \
			-d "$assay_type"
	
	elif [[ $# -eq 8 ]]; then
	
		local plus_shift="$7"
		local minus_shift="$8"

		python $1 \
			-ifrag $frag_file \
			-op $out_bw \
			-c $chrom_sizes \
			-g $ref_genome_file \
			-d $assay_type \
			-p $plus_shift \
			-m $minus_shift
			
	else

		echo -e "Wrong number of input parameters.\n Correct are 6 or 8: [cpn_binary, frag_file, out_bw, ref_genome_file, chrom_sizes, assay_type, (plus_shift, minus_shift)]"
		return 1

	fi

	# Rename ugly output file
	mv "${out_bw}_unstranded.bw" "${out_bw}.bw"

}


##### SCRIPT #####

source $HOME/.bash_profile
load-conda
conda activate ian

# Set paths
PROJECT_PATH=".."
PROJECT_PATH="$(realpath ${PROJECT_PATH})"
cd $PROJECT_PATH

TOOL="${1:-"CBN"}" # ChromBPNet

MODE="${2:-"borgs"}" # hca brain organoids

# ... for different modes
case $MODE in

	borgs )
		DATASET="hca_brain-organoids"
		CT_MAP_ID="approach_2024-09-12" # Only used to name files
		FRAG_FILES_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/fragment-files/grouped/${CT_MAP_ID}"

		case "$TOOL" in

			CBN )
				OUT_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/coverages/chrombpnet/${CT_MAP_ID}"
				;;

			SAFT )
				OUT_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/coverages/sc-atac-frag-tools/${CT_MAP_ID}"
				;;

			* ) echo "Choose coverage tool in argument 1 [CBN | SAFT]."; exit 1 ;;

		esac

		;;

	toy )
		echo "Not implemented"
		;;

	borgs_small )
		echo "Not implemented"
		#FRAG_FILES=(\
			#"melw_Glia.tsv" \
			#"melw_Midbrain-EN.tsv" \
			#"melw_UL-EN.tsv" \
			#"pelm_Glia.tsv" \
			#"pelm_Midbrain-EN.tsv" \
			#"pelm_UL-EN.tsv" \
			#"ualf_Glia.tsv" \
			#"ualf_Midbrain-EN.tsv" \
			#"ualf_UL-EN.tsv" \
			#"zoxy_Glia.tsv" \
			#"zoxy_Midbrain-EN.tsv" \
			#"zoxy_UL-EN.tsv")
		;;

	* )
		echo "Bug in code. Check MODE variable."
		exit 1
		;;

esac

# What to do if output directories already populated?
if [[ -n "$(find "$OUT_DIR" -mindepth 1 -print -quit 2>/dev/null)" ]]; then

	while true; do

		echo "Output directory <${OUT_DIR}> already contains files or folders."
		read -r -p "Delete contents? [Y,n]: " del

		case "$del" in

			[Yy]* ) rm -rf "${OUT_DIR:?}"; break ;;
			[Nn]* ) echo "I won't work with confounding folders or files."; exit 1 ;;
			* ) echo "Please answer [Y | n]." ;;

		esac

	done

fi

mkdir -p "$OUT_DIR"

find "$FRAG_FILES_DIR" -mindepth 1 -maxdepth 1 \( -type f -o -type l \) -print0 |
while IFS= read -r -d '' frag_file; do

	echo "Processing: ${frag_file}"

	case $TOOL in

		CBN )
			create_bw_from_frag_file_cbn \
				"${HOME}/.conda/envs/ian/lib/python3.9/site-packages/chrombpnet/helpers/preprocessing/reads_to_bigwig.py" \
				"${FRAG_FILES_DIR}/${frag_file}" \
				"${OUT_DIR}/${frag_file%.*}" \
				"data/GRCh38-p14/hg38.fa" \
				"data/GRCh38-p14/hg38.chrom.sizes" \
				"ATAC"
			;;

		SAFT )
			#TODO: implement
			echo "Not implemented yet"
			;;
			#scatac_fragment_tools bigwig \
				#-i "${FRAG_FILES_DIR}/${F}" \
				#-o "${OUT_DIR_SAFT}/${F%.*}.bw" \
				#-c "data/GRCh38-p14/hg38.chrom.sizes" \
				#-x \
				#-normalize False \
				#-scaling 1.0
				#–cut-sites False # Use 1 bp Tn5 cut sites (start and end of each fragment) instead of 
				# 		   whole fragment length for coverage calculation
			
		* )
			echo "Bug in code";
			exit 1
			;;

	esac

done

