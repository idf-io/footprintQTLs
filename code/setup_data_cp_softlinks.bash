#!/usr/bin/env bash
#
# Setup data structure from private sources by copying and making softlinks.

PROJECT_DIR=".."
PROJECT_DIR="$(realpath "$PROJECT_DIR")"
cd $PROJECT_DIR

source "code/glob_vars.bash"

BORGS_ROOT="/omics/groups/OE0540/internal/projects/HCA_organoid_2/cemm_sabrina-20Jul2022/"

cp "${BORGS_ROOT}/eQTL_mapping/eSNPs_significant_all_celltypes_HVGs.tsv" "$PRECOMPUTED_EQTLS_TSV"

#rm \
	#"$PRECOMPUTED_EQTLS_TSV"
