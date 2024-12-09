#!/usr/bin/env bash
#
# Calculate how many variants are within all- and cell-type-specific peaks.
#
#BSUB -R "rusage[mem=10G]"
#BSUB -q short
#BSUB -cwd "/omics/groups/OE0540/internal_temp/users/fichtner/projects/footprintQTL"
#BSUB -J "calculate_n_variants_within_peaks"
#BSUB -o "/omics/groups/OE0540/internal_temp/users/fichtner/projects/footprintQTL/code/bsub/logs/%I.out"
#BSUB -e "/omics/groups/OE0540/internal_temp/users/fichtner/projects/footprintQTL/code/bsub/logs/%I.err"

# Time: 1' in bworker

### Setup ###
set -euo pipefail

if [[ "$(basename $(pwd))" == "calculations" ]]; then

	cd "../.."

fi

module load bedtools/2.29.2

source "code/glob_vars.bash"
#  1/2. SNP_LOCS_BED, ATAC_CHROM_ACCESS_METADATA_DIR 
#  3/4. SNP_LOCS_BED, ATAC_CHROM_ACCESS_METADATA_DIR, MEQTL_CA_QTLS_DIR, CT_MAP_ID
#  5. SNP_LOCS_BED, FOOTPRINTS_METADATA_TSV

source "${HOME}/.bash_profile"
load-micromamba
micromamba activate "$MAIN_ENV"



### Calculations ###

# 1. SNPs within all peaks

echo -e -n "\n------\n\n1. SNPs within all peaks.\nNr: "
bedtools intersect -a "$SNP_LOCS_BED" -b "${ATAC_CHROM_ACCESS_METADATA_DIR}/peaks_all.bed" -wa -wb |
	tee "${ATAC_CHROM_ACCESS_METADATA_DIR}/snps_overlapping_peaks_all.bed" |
	wc -l


# 2. SNPs within all peaks +- 1kbp window at each flank

echo -e -n "\n------\n\n1. SNPs within all peaks +- 1kbp around each flank.\nNr: "
bedtools window -a "$SNP_LOCS_BED" -b "${ATAC_CHROM_ACCESS_METADATA_DIR}/peaks_all.bed" -w 1000 |
	tee "${ATAC_CHROM_ACCESS_METADATA_DIR}/snps_overlapping_peaks_all_pm_1kbp-window.bed" |
	wc -l


# 3. SNPs within caQTL peaks (cell-type-level)

echo -e "\n------\n\n2. SNPs wihtin cell-type-specific peaks used for caQTL calling:"

while IFS= read -r -d '' cell_type_dir; do

	cell_type="$(basename "$cell_type_dir")"

	echo -e "\nProcessing file: ${cell_type_dir}/peaks.bed"

	out_bed="${ATAC_CHROM_ACCESS_METADATA_DIR}/selected_peaks_ca-qtls/${CT_MAP_ID}/${cell_type}/snps_overlapping_peaks.bed"
	mkdir -p "$(dirname "$out_bed")"

	echo -n "Number of snps within peaks: "
	bedtools intersect -a "$SNP_LOCS_BED" \
		-b <(awk -F'\t' -v OFS='\t' 'NR > 1 {print $2, $3, $4, $1}' "${cell_type_dir}/peak_locations_all.tsv" | sort -k1,1 -k2,2n) \
		-wa -wb |
		tee "$out_bed" |
		wc -l

done < <(find "$CA_QTLS_MEQTL_INPUT_DIR" -mindepth 1 -maxdepth 1 -type d -print0)


# 4. snps within caQTL peaks +-1kbp window (cell-type-level)

echo -e "\n------\n\n3. SNPs within cell-type-specific peaks used for caQTL calling +- 1kbp around each flank:"

while IFS= read -r -d '' cell_type_dir; do

	cell_type="$(basename "$cell_type_dir")"

	echo -e "\nProcessing file: ${cell_type_dir}/peaks.bed"

	out_bed="${ATAC_CHROM_ACCESS_METADATA_DIR}/selected_peaks_ca-qtls/${CT_MAP_ID}/${cell_type}/snps_overlapping_peaks_pm_1kbp-window.bed"
	mkdir -p "$(dirname "$out_bed")"

	echo -n "Number of snps within peaks +-1kbp: "
	bedtools window -a "$SNP_LOCS_BED" \
		-b <(awk -F'\t' -v OFS='\t' 'NR > 1 {print $2, $3, $4, $1}' "${cell_type_dir}/peak_locations_all.tsv" | sort -k1,1 -k2,2n) \
		-w 1000 |
		tee "$out_bed" |
		wc -l

done < <(find "$CA_QTLS_MEQTL_INPUT_DIR" -mindepth 1 -maxdepth 1 -type d -print0)


# 4. snps within footprint peaks (cell-type-level)

echo -e "\n------\n\n4. SNPs within cell-type-specific peaks used in footprints:"

while IFS= read -r -d '' cell_type_dir; do

	cell_type="$(basename "$cell_type_dir")"

	echo -e "\nProcessing file: ${cell_type_dir}/peaks.bed"

	echo -n "Number of snps within peaks: "
	bedtools intersect -a "$SNP_LOCS_BED" -b "${cell_type_dir}/peaks.bed" -wa -wb | tee "${cell_type_dir}/snps_overlapping_peaks.bed" | wc -l

done < <(find "$FOOTPRINTS_METADATA_DIR" -mindepth 1 -maxdepth 1 -type d -print0)
