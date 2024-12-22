#!/usr/bin/env bash
set -euo pipefail

cd "${HOME}/projects/footprintQTL"


## Fragment file

FRAG_FILE_DIR_OLD="/omics/groups/OE0540/internal/projects/HCA_organoid_2/cemm_sabrina-20Jul2022/outputs"
FRAG_FILE_DIR_NEW="data/datasets/hca_brain-organoids/atac-seq/fragment-files"

mkdir -p \
	$FRAG_FILE_DIR_NEW


while IFS= read -r -d '' sample_ref_path; do 

	sample_ref="$(basename "$sample_ref_path")"

	echo $sample_ref_path
	echo $sample_ref
	base_path="${FRAG_FILE_DIR_OLD}/${sample_ref}/${sample_ref}/outs/atac_fragments"

	if [[ "$sample_ref" =~ sSL[0-9]{4}[a-zA-Z]* ]]; then

		sample=$BASH_REMATCH

		source_gz="${base_path}.tsv.gz"
		target_gz="${FRAG_FILE_DIR_NEW}/${sample}.tsv.gz"

		source_tbi="${base_path}.tsv.gz.tbi"
		target_tbi="${FRAG_FILE_DIR_NEW}/${sample}.tsv.gz.tbi"

		ln -sf $source_gz $target_gz

		ln -sf $source_tbi $target_tbi

	else

		echo "File <${sample_ref}> didn't have an atac fragments file, SKIPPED"

	fi

done < <(find $FRAG_FILE_DIR_OLD -maxdepth 1 -type d ! -name "." ! -regex "^${FRAG_FILE_DIR_OLD}$" -print0)
