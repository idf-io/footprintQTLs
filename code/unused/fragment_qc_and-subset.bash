#! /usr/bin/env bash

set -eu

PROJECT_PATH="/home/fichtner/projects/gemmo-tools/"
cd $PROJECT_PATH

DATA_PATH="/omics/groups/OE0540/internal/projects/HCA_organoid_2/cemm_sabrina-20Jul2022/"

declare -A metadata_map

while IFS=$'\t' read -r cell_id sample donor donor_id cell_type; do
	metadata_map["${cell_id}"]="$sample,$donor,$donor_id,$cell_type"
done < "data/intermediate-data/hca_brain_organoids/cell-types_map.tsv"

echo CreatedMap


while IFS=$'\t' read -r chr start end barcode counts; do

	if ! [[ $chr =~ ^# ]] \ # Skip comments
		&& [[ $(( end - start )) -gt 10 ]] \ # Min frag length
		&& [[ $(( end - start )) -lt 200 ]]\ # Max frag length
		[[ ]] # Pass cell-lvl QC
		then

		#printf "Barcode: %s, metadata: %s" "$barcode" "${metadata_map["${barcode}_sSL0146"]}"
		printf "%s\t%s\t%s\t%s\t%s\n" \
			"$chr" \
			"$start" \
			"$end" \
			"${barcode}_$(echo "${metadata_map["${barcode}_sSL0146"]}" | cut -d "," -f 2)" \
			"$counts" \
			>> "${PROJECT_PATH}data/intermediate-data/hca_brain_organoids/new.tsv"

	fi

done < <(zcat "${DATA_PATH}outputs/sSL0146_BrainO_R4_F_10xM_Multiome/sSL0146_BrainO_R4_F_10xM_Multiome/outs/atac_fragments.tsv.gz")


		#echo -e "$chr\t$start\t$end\t${barcode}_$(echo "${metadata_map["${barcode}_sSL0146"]}" | cut -d "," -f 2)\t$counts" >> "${project_path}data/intermediate-data/hca_brain_organoids/new.tsv"
		#
#find $DATA_PATH -maxdepth 1 -type d | while 

# Include delete file if exists
#SAMPLES=$((find . -maxdepth 1 -type f))
#SAMPLES=
#printf "cell_id: %s, values: %s\n" "$cell_id" "$(echo "${metadata_map["$cell_id"]}" | cut -d "," -f 2)"
#echo "${metadata_map["TTTGGTGCATAGCAGG-1_sSL0146"]}" | cut -d "," -f 2
