import anndata as ad
import os
from pathlib import Path
import re

PROJECT_PATH = '/home/fichtner/projects/gemmo-tools'
DATA_PATH = '/omics/groups/OE0540/internal/projects/HCA_organoid_2/cemm_sabrina-20Jul2022/'
RNA_AD = 'outputs_allsamples/sabrina_allsamples_rna_final_after_atac.h5ad'
OUTPUT_FILE = 'data/datasets/hca_brain-organoids_data/group_to_cell_barcode.tsv'

os.chdir(PROJECT_PATH)


def regex_get(text, pattern=None):
    
    pattern_map = {'sample_name': r'sSL\d*',
                   'barcode': r'^\D*-\d*'}

    if pattern in pattern_map.keys():
        pattern = pattern_map[pattern]
    
    match = re.search(pattern, text)
    return match.group()


borgs_rna = ad.read_h5ad(os.path.join(DATA_PATH + RNA_AD), backed='r')
cells_coldata = borgs_rna.obs[['sample', 'donor', 'celltype_predicted_vertesy']].copy()
cells_coldata['barcode'] = [regex_get(i, 'barcode') for i in cells_coldata.index.tolist()]

if not os.path.isdir(Path(OUTPUT_FILE).parent):
    os.mkdir(Path(OUTPUT_FILE).parent)
    
if os.path.isfile(OUTPUT_FILE):
    os.remove(OUTPUT_FILE)

with open(OUTPUT_FILE, 'a') as out_tsv:
    
    out_tsv.write('sample\tcell_type\tcell_barcode\n')
    
    for row in cells_coldata.itertuples(index=True, name='Pandas'):
        
        out_tsv.write('{}\t{}\t{}\n'.format(row.sample, row.donor + '_' + row.celltype_predicted_vertesy, row.barcode))