#! /usr/bin/env bash

set -euo pipefail

PROJECT_PATH=/home/fichtner/projects/gemmo-tools
cd $PROJECT_PATH

DATA_PATH=/omics/groups/OE0540/internal/projects/HCA_organoid_2/cemm_sabrina-20Jul2022
DATA_PATH=data/datasets/hca_brain-organoids_small
OUTPUT_FOLDER=data/datasets/hca_brain-organoids_grouped_scTools
SAMPLE_FRAG_MAP=data/datasets/hca_brain-organoids_data/sample_to_fragment.tsv
GROUP_BARCODE_MAP=data/datasets/hca_brain-organoids_data/group_to_cell_barcode.tsv


if [[ ! -d data/datasets/hca_brain-organoids_data ]]; then
	mkdir data/datasets/hca_brain-organoids_data
fi


# Create sample->fragment map file
if [[ -f $SAMPLE_FRAG_MAP ]]; then
	rm $SAMPLE_FRAG_MAP
fi

SAMPLES_LONG=($(ls $DATA_PATH/outputs))
SAMPLES_LONG=(\
	"sSL0136_BrainO_R4_F_10xM_Multiome" \
	"sSL0146_BrainO_R4_F_10xM_Multiome" \
	"sSL0170_BrainO_R4_F_10xM_Multiome")
SAMPLES=()

printf 'sample\tpath_to_fragment_file\n' >> $SAMPLE_FRAG_MAP

for S in "${SAMPLES_LONG[@]}"; do
	if [[ $S =~ (^sSL[0-9]{4}[A-Z]?)_BrainO_(.*)_10xM_Multiome$ ]]; then
		sample_id=${BASH_REMATCH[1]}
		sample_code=${BASH_REMATCH[2]}
		SAMPLES+=($sample_id)
		printf '%s\t%s/outputs/%s_BrainO_%s_10xM_Multiome/%s_BrainO_%s_10xM_Multiome/outs/atac_fragments.tsv.gz\n' "$sample_id" "${DATA_PATH}" "${sample_id}" "${sample_code}" "${sample_id}" "${sample_code}" >> $SAMPLE_FRAG_MAP
	fi
done


# Create group->barcode map file

if [[ -f $GROUP_BARCODE_MAP ]]; then
	rm $GROUP_BARCODE_MAP
fi

python code/create_group2barcode_map.py




# Create grouped frags files
 scatac_fragment_tools split \
    -f $SAMPLE_FRAG_MAP \
    -b $GROUP_BARCODE_MAP \
    -c $PROJECT_PATH/data/GRCh38-p14/chrom_lengths.chrom.sizes \
    -o $OUTPUT_FOLDER
    -v \
    -n 32 \ # All available
    # –clear_temp True \
    # -cell_type_column cell_group # Didn't work as expected
