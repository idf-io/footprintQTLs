#!/usr/bin/env bash

set -eoux

DATA_PATH=/omics/groups/OE0540/internal/projects/HCA_organoid_2/cemm_sabrina-20Jul2022
SAMPLES=($(find "$DATA_PATH" -maxdepth 1 -mindepth 1 -type d -printf "%P\n"))

#sSL0146_BrainO_R4_F_10xM_Multiome/
#sSL0146_BrainO_R4_F_10xM_Multiome/outs/atac_fragments.tsv.gz

folders_string=$(IFS="\t"; echo "${SAMPLES[*]}")
echo -e "folders_string" >> custom_region_fragments_CT1.tsv


for d in "${SAMPLES[@]}"; do
	tabix $d chr5:140011317-140012000 >> custom_region_fragments_CT1.tsv
done


