MAIN_ENV="main06"

DATASET="hca_brain-organoids"
CT_MAP_JSON="config/cell-type_groupings/${DATASET}/approach_2024-09-12.json"
CT_MAP_ID="$(basename "${CT_MAP_JSON%.json}")"
COVERAGE_TOOL="chrombpnet"

GROUPED_FRAG_FILES_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/fragment-files/grouped/${CT_MAP_ID}" # Only used in bash
GROUPED_BIGWIG_FILES_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/coverages/${COVERAGE_TOOL}/grouped/${CT_MAP_ID}" # Only used in bash


PRECOMPUTED_EQTLS_TSV="data/datasets/${DATASET}/rna-seq/gene-expression/eqtls/eqtls_fdr_all_cell-types_hvgs.tsv"


ATAC_PEAKS_H5AD_NEW="data/datasets/hca_brain-organoids/atac-seq/chromatin_accessibility/peak-matrix_rna-qc-cells_norm-reads-in-tss.h5ad"
ATAC_CHROM_ACCESS_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/chromatin-accessibility"
ATAC_CHROM_ACCESS_METADATA_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/chromatin-accessibility/metadata"
CA_QTLS_MEQTL_INPUT_DIR="data/intermediate-data/datasets/${DATASET}/matrix-eqtl/chromatin-accessibility/${CT_MAP_ID}"

FOOTPRINT_APPROACH="js_divergence"
FOOTPRINTS_DIR="results/datasets/${DATASET}/atac-seq/footprints/${FOOTPRINT_APPROACH}/${CT_MAP_ID}" # bash only
FOOTPRINTS_EDA="results/datasets/${DATASET}/atac-seq/footprints/${FOOTPRINT_APPROACH}/${CT_MAP_ID}/eda"
FOOTPRINTS_METADATA_DIR="results/datasets/${DATASET}/atac-seq/footprints/${FOOTPRINT_APPROACH}/${CT_MAP_ID}/metadata"


GENOTYPES_TSV="data/datasets/${DATASET}/genotype/genotype_NA.tsv"
GENOTYPES_VCF="data/datasets/${DATASET}/genotype/genotype.vcf"

SNP_LOCS_BED="data/intermediate-data/datasets/${DATASET}/genotype/snp_locations.bed"


MATRIX_EQTL_INPUT_DIR="data/intermediate-data/datasets/${DATASET}/matrix-eqtl/footprints/${FOOTPRINT_APPROACH}/${CT_MAP_ID}"
MATRIX_EQTL_OUTPUT_DIR="results/datasets/${DATASET}/matrix-eqtl/footprints/${FOOTPRINT_APPROACH}/${CT_MAP_ID}"

REF_GENOME_FASTA="data/GRCh38-p14/hg38.fa"
CHROM_SIZES="data/GRCh38-p14/hg38.chrom.sizes"
