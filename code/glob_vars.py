import os

MAIN_ENV = 'main06'

DATASET = 'hca_brain-organoids'
CT_MAP_JSON = f'config/cell-type_groupings/{DATASET}/approach_2024-09-12.json'
CT_MAP_ID = os.path.splitext(os.path.basename(CT_MAP_JSON))[0]
COVERAGE_TOOL = 'chrombpnet'
ALGORITHMS = ('js_divergence', 'counts')


RNA_H5AD_OLD = '/omics/groups/OE0540/internal/projects/HCA_organoid_2/cemm_sabrina-20Jul2022/outputs_allsamples/sabrina_allsamples_rna_final_after_atac.h5ad' # Only python
RNA_H5AD_NEW = f'data/datasets/{DATASET}/rna-seq/gene_expression/gene-expression_qc_after-atac.h5ad' # Only python
PRECOMPUTED_EQTLS_TSV = f'data/datasets/{DATASET}/rna-seq/gene-expression/eqtls/eqtls_fdr_all_cell-types_hvgs.tsv'
EQTLS_DIR = f'data/intermediate-data/datasets/{DATASET}/rna-seq/gene-expression/eqtls'


GROUPED_FRAG_FILES_DIR = f'data/intermediate-data/datasets/{DATASET}/atac-seq/fragment-files/grouped/{CT_MAP_ID}'
GROUPED_BIGWIG_FILES_DIR = f'data/intermediate-data/datasets/{DATASET}/atac-seq/coverages/{COVERAGE_TOOL}/grouped/{CT_MAP_ID}'


ATAC_PEAKS_H5AD_OLD = f'/omics/groups/OE0540/internal/projects/HCA_organoid_2/cemm_sabrina-20Jul2022/peakmatrix_RNA_QC_cells_reads_in_tss_norm.h5ad' # Only python
ATAC_PEAKS_H5AD_NEW = f'data/datasets/{DATASET}/atac-seq/chromatin_accessibility/peak-matrix_rna-qc-cells_norm-reads-in-tss.h5ad'
ATAC_PEAKS_PROCESSED_H5AD = f'data/intermediate-data/datasets/{DATASET}/atac-seq/chromatin-accessibility/adata/peak_matrix.h5ad'
ATAC_CHROM_ACCESS_DIR = f'data/intermediate-data/datasets/{DATASET}/atac-seq/chromatin-accessibility'
ATAC_CHROM_ACCESS_METADATA_DIR = f'data/intermediate-data/datasets/{DATASET}/atac-seq/chromatin-accessibility/metadata'
SELECT_PEAKS_TSV_DIR = f'data/intermediate-data/datasets/{DATASET}/atac-seq/chromatin-accessibility/metadata/select_peaks/{CT_MAP_ID}'
CA_QTLS_MEQTL_INPUT_DIR = f'data/intermediate-data/datasets/{DATASET}/matrix-eqtl/chromatin-accessibility/{CT_MAP_ID}'





FOOTPRINT_APPROACH = 'js_divergence'
FOOTPRINTS_DIR = f'results/datasets/{DATASET}/atac-seq/footprints/'
#FOOTPRINTS_EDA = f'results/datasets/{DATASET}/atac-seq/footprints/{FOOTPRINT_APPROACH}/{CT_MAP_ID}/eda'
#FOOTPRINTS_METADATA_DIR = f'results/datasets/{DATASET}/atac-seq/footprints/{FOOTPRINT_APPROACH}/{CT_MAP_ID}/metadata'


GENOTYPES_TSV = f'data/datasets/{DATASET}/genotype/genotype_NA.tsv'
GENOTYPES_PROCESSED_TSV = f'data/intermediate-data/datasets/{DATASET}/genotype/genotype_NA.tsv'
GENOTYPES_VCF = f'data/datasets/{DATASET}/genotype/genotype.vcf'
GENOTYPE_PCS_TSV = f'data/datasets/{DATASET}/genotype/genotype_pcs_T.tsv'

GENOTYPE_DIR = f'data/intermediate-data/datasets/{DATASET}/genotype'
SNP_LOCS_BED = f'data/intermediate-data/datasets/{DATASET}/genotype/snp_locations.bed'




MATRIX_EQTL_INPUT_DIR = f'data/intermediate-data/datasets/{DATASET}/matrix-eqtl/footprints'
MATRIX_EQTL_OUTPUT_DIR = f'results/datasets/{DATASET}/matrix-eqtl/footprints/{FOOTPRINT_APPROACH}/{CT_MAP_ID}'

CHROM_SIZES="data/GRCh38-p14/hg38.chrom.sizes"
