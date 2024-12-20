library(tools)

DATASET = 'hca_brain-organoids'
CT_MAP_JSON = paste0('config/cell-type_groupings/', DATASET, '/approach_2024-09-12.json')
CT_MAP_ID = file_path_sans_ext(basename(CT_MAP_JSON))

FOOTPRINT_APPROACH = 'js_divergence'
#FOOTPRINTS_DIR = f'results/datasets/{DATASET}/atac-seq/footprints/{FOOTPRINT_APPROACH}/{CT_MAP_ID}'
#FOOTPRINTS_EDA = f'results/datasets/{DATASET}/atac-seq/footprints/{FOOTPRINT_APPROACH}/{CT_MAP_ID}/eda'

GENOTYPES_TSV = sprintf('data/datasets/%s/genotype/genotype_NA.tsv', DATASET)
#GENOTYPE_PCS_TSV = f'data/datasets/{DATASET}/genotype/genotype_pcs_T.tsv' # python only

MATRIX_EQTL_INPUT_DIR = paste0('data/intermediate-data/datasets/', DATASET, '/matrix-eqtl/footprints/', FOOTPRINT_APPROACH, '/', CT_MAP_ID)
MATRIX_EQTL_OUTPUT_DIR = paste0('results/datasets/', DATASET, '/matrix-eqtl/footprints/', FOOTPRINT_APPROACH, '/', CT_MAP_ID)
