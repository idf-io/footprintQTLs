MAIN_ENV="main06"

DATASET="hca_brain-organoids"
CT_MAP_JSON="config/cell-type_groupings/${DATASET}/approach_2024-09-12.json"
CT_MAP_ID="$(basename "${CT_MAP_JSON%.json}")"
COVERAGE_TOOL="chrombpnet"
ALGORITHMS=("js_divergence" "counts")

GENOTYPES_TSV="data/datasets/${DATASET}/genotype/genotype_NA.tsv"
GENOTYPES_VCF="data/datasets/${DATASET}/genotype/genotype.vcf"
SNP_LOCS_BED="data/intermediate-data/datasets/${DATASET}/genotype/snp_locations.bed"
REF_GENOME_FASTA="data/GRCh38-p14/hg38.fa"
CHROM_SIZES="data/GRCh38-p14/hg38.chrom.sizes"
PRECOMPUTED_EQTLS_TSV="data/datasets/${DATASET}/rna-seq/gene-expression/eqtls/eqtls_fdr_all_cell-types_hvgs.tsv"

ATAC_PEAKS_H5AD_OLD="/omics/groups/OE0540/internal/projects/HCA_organoid_2/cemm_sabrina-20Jul2022/peakmatrix_RNA_QC_cells_reads_in_tss_norm.h5ad"
ATAC_CHROM_ACCESS_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/chromatin-accessibility"


ATAC_PEAKS_H5AD_NEW="data/datasets/hca_brain-organoids/atac-seq/chromatin_accessibility/peak-matrix_rna-qc-cells_norm-reads-in-tss.h5ad"
ATAC_PEAKS_PROCESSED_H5AD="data/intermediate-data/datasets/${DATASET}/atac-seq/chromatin-accessibility/adata/peak_matrix.h5ad"
ATAC_CHROM_ACCESS_METADATA_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/chromatin-accessibility/metadata"
SELECT_PEAKS_TSV_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/chromatin-accessibility/metadata/select_peaks/${CT_MAP_ID}"

GROUPED_FRAG_FILES_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/fragment-files/grouped/${CT_MAP_ID}" 
GROUPED_BIGWIG_FILES_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/coverages/${COVERAGE_TOOL}/grouped/${CT_MAP_ID}" 

FOOTPRINTS_DIR="results/datasets/${DATASET}/atac-seq/footprints"
#FOOTPRINT_APPROACH="js_divergence"
#FOOTPRINTS_EDA="results/datasets/${DATASET}/atac-seq/footprints/${FOOTPRINT_APPROACH}/${CT_MAP_ID}/eda"
#FOOTPRINTS_METADATA_DIR="results/datasets/${DATASET}/atac-seq/footprints/${FOOTPRINT_APPROACH}/${CT_MAP_ID}/metadata"

MATRIX_EQTL_INPUT_DIR="data/intermediate-data/datasets/${DATASET}/matrix-eqtl"
MATRIX_EQTL_INPUT_FOOTPRINTS_DIR="data/intermediate-data/datasets/${DATASET}/matrix-eqtl/footprints"
MATRIX_EQTL_OUTPUT_DIR="results/datasets/${DATASET}/matrix-eqtl"
MATRIX_EQTL_OUTPUT_FOOTPRINTS_DIR="results/datasets/${DATASET}/matrix-eqtl/footprints"




CA_QTLS_MEQTL_INPUT_DIR="data/intermediate-data/datasets/${DATASET}/matrix-eqtl/chromatin-accessibility/${CT_MAP_ID}"
