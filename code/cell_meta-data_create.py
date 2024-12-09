import anndata as ad
import os
from pathlib import Path

PROJECT_PATH = '/home/fichtner/projects/gemmo-tools'
DATA_PATH = '/omics/groups/OE0540/internal/projects/HCA_organoid_2/cemm_sabrina-20Jul2022/'
RNA_AD = 'outputs_allsamples/sabrina_allsamples_rna_final_after_atac.h5ad'
OUTPUT_FILE = 'data/intermediate-data/hca_brain_organoids/cell-types_map.tsv'

os.chdir(PROJECT_PATH)

borgs_rna = ad.read_h5ad(os.path.join(DATA_PATH + RNA_AD))
cell_metadata = borgs_rna.obs[['sample', 'donor', 'donor_id', 'celltype_predicted_vertesy']].copy()

if not os.path.isdir(Path(OUTPUT_FILE).parent):
    os.mkdir(Path(OUTPUT_FILE).parent)

with open(OUTPUT_FILE, 'w') as out_tsv:

    for row in cell_metadata.itertuples(index=True, name='Pandas'):
        out_tsv.write(f'{row.Index}\t{row.sample}\t{row.donor}\t{row.donor_id}\t{row.celltype_predicted_vertesy}\n')
