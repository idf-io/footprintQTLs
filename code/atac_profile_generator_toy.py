import os
import pandas as pd
import pysam
import anndata as ad
import re
from scipy.sparse import csr_matrix, coo_matrix
import numpy as np

PROJECT_PATH = '/home/fichtner/projects/gemmo-tools'
DATA_PATH = '/omics/groups/OE0540/internal/projects/HCA_organoid_2/cemm_sabrina-20Jul2022/'
DATA_PATH_ALT = '/home/fichtner/data/datasets/hca-brain-organoids_small'
#SAMPLES = [f for f in os.listdir(DATA_PATH + 'outputs') if os.path.isdir(os.path.join(DATA_PATH + 'outputs', f))]
SAMPLES = ['sSL0146_BrainO_R4_F_10xM_Multiome']
RNA_AD = 'outputs_allsamples/sabrina_allsamples_rna_final_after_atac.h5ad'
OUT_FILE = 'CD14_filtered_fragments_10k.tsv'

os.chdir(PROJECT_PATH)



def regex_get(text, pattern=None):
    
    pattern_map = {'sample_name': r'sSL\d*',
                   'barcode': r'^\D*-\d*'}

    if pattern in pattern_map.keys():
        pattern = pattern_map[pattern]
    
    match = re.search(pattern, text)
    return match.group()



borgs_rna = ad.read_h5ad(DATA_PATH + RNA_AD)

cells_coldata = borgs_rna.obs[['sample', 'donor', 'donor_id', 'celltype_predicted_vertesy']].copy()
cells_coldata.rename(columns={'celltype_predicted_vertesy': 'cell_type'}, inplace=True)
cells_coldata['barcode'] = [regex_get(i, 'barcode') for i in cells_coldata.index.tolist()]

qc_cells = cells_coldata.index.tolist()

cell_types = tuple(set(cells_coldata['cell_type']))
cell_types = sorted(cell_types, key=str.upper)

samples = tuple(set(cells_coldata['sample']))
samples = sorted(samples, key=str.upper)

donors = tuple(set(cells_coldata['donor']))
donors = sorted(donors, key=str.upper)

# groups = donors
groups = ['bima']
cells_coldata


regions = [['cd14', 'chr5', 140011317, 140012000]] # CD14, 1-based
print(*regions[0])

min_frag_size = 10
max_frag_size = 200


if os.path.exists(OUT_FILE):
    print('Output file already exists')
    os.remove(OUT_FILE)


frags = pd.DataFrame(columns=['chr', 'start', 'end', 'barcode', 'counts', 'sample', 'donor', 'donor_id', 'cell_type'])

for s in SAMPLES:

    sample_path = os.path.join(DATA_PATH_ALT, 'outputs', s, s, 'outs/atac_fragments_10k.tsv.gz')
    sample_name = regex_get(s, 'sample_name')

    tbx = pysam.TabixFile(sample_path, parser=pysam.asTuple())

    for row in tbx.fetch(*regions[0][1:0]):

        cell_id = '{}_{}'.format(row[3], sample_name)


        if cell_id in qc_cells and (int(row[2]) - int(row[1])) > min_frag_size and (int(row[2]) - int(row[1])) < max_frag_size:

            entry = {'chr': row[0],
                      'start': int(row[1]),
                      'end': int(row[2]),
                      'barcode': row[3],
                      'counts': row[4],
                      'sample': cells_coldata.loc[cell_id, 'sample'],
                      'donor': cells_coldata.loc[cell_id, 'donor'],
                      'donor_id': cells_coldata.loc[cell_id, 'donor_id'], 
                      'cell_type': cells_coldata.loc[cell_id, 'cell_type']}
            frags = frags._append(entry, ignore_index=True)
            
            
print(frags[0:15])
