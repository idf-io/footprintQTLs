import os
import gzip
import anndata as ad
import re

PROJECT_PATH = '/home/fichtner/projects/gemmo-tools'
os.chdir(PROJECT_PATH)

DATA_PATH = 'data/datasets/hca_brain-organoids' # softlink to /omics/groups/OE0540/internal/projects/HCA_organoid_2/cemm_sabrina-20Jul2022/
DATA_PATH = 'data/datasets/hca_brain-organoids_small'
SAMPLE_CODES = [f for f in os.listdir(DATA_PATH + '/outputs')]
# SAMPLES = ['sSL0146_BrainO_R4_F_10xM_Multiome', 'sSL0136_BrainO_R3_B_10xM_Multiome', 'sSL0170_BrainO_R4_F_10xM_Multiome']
RNA_AD = 'outputs_allsamples/sabrina_allsamples_rna_final_after_atac.h5ad'
OUT_FOLDER = os.path.join(PROJECT_PATH, 'data/datasets/hca_brain-organoids_grouped/')



def ct_format(cell_type):
    return cell_type.replace(' ', '-').replace('.', '')

def regex_get(text, pattern=None):
    
    pattern_map = {'sample_name': r'sSL\d*[A-Za-z]*',
                   'barcode': r'^\D*-\d*'}

    if pattern in pattern_map.keys():
        pattern = pattern_map[pattern]
    
    match = re.search(pattern, text)
    return match.group()


def close_multiple_files(file_objects):
    """
    Close multiple files.
    
    :param file_objects: Dictionary of file objects with file names as keys
    """
    for file in file_objects.values():
        file.close()

        
def open_multiple_files(file_names, mode='w'):
    """
    Open multiple files simultaneously as specified by a list of strings.
    
    :param file_names: List of file names (without extensions)
    :param mode: Mode in which to open the files (default is 'w' for write mode)
    :return: Dictionary of file objects with file names as keys
    """
    
    # Check for format
    for name in file_names:
        if (' ' in name) or ('.' in name):
            print('Invalid characters in `file_names`. E.g. space or point.')
            return {}
            
   
    file_objects = {}
    
    for name in file_names:
        
        file_path = os.path.join(OUT_FOLDER, f'{name}_atac_fragments.tsv')
        
        # Ensure file doesn't exists, else abort
        try:
            if os.path.isfile(file_path):
                raise ValueError(f'File: `{file_path}` already exists. ABORTING OPERATION.')
                                                 
        except ValueError as e:
            print(e)
            close_multiple_files(file_objects)
            return {}

        file_objects[name] = open(file_path, mode)

                                 
    return file_objects

# borgs_rna = ad.read_h5ad(DATA_PATH + RNA_AD, backed='r')
borgs_rna = ad.read_h5ad('data/datasets/hca_brain-organoids/' + RNA_AD, backed='r')


# Metadata of cells that passed QC
cells_coldata = borgs_rna.obs[['sample', 'donor', 'celltype_predicted_vertesy']].copy()
cells_coldata['cell_type'] = cells_coldata['celltype_predicted_vertesy'].apply(ct_format)
cells_coldata = cells_coldata.drop('celltype_predicted_vertesy', axis=1)
cells_coldata['barcode'] = [regex_get(i, 'barcode') for i in cells_coldata.index.tolist()]

ann_idxs = cells_coldata.index.tolist()

cell_types = tuple(set(cells_coldata['cell_type']))
cell_types = sorted(cell_types, key=str.upper)

samples = tuple(set(cells_coldata['sample']))
samples = sorted(samples, key=str.upper)

donors = tuple(set(cells_coldata['donor']))
donors = sorted(donors, key=str.upper)

# Note: (indexes ~ cell id) = barcode_sample --> are unique but barcodes aren't


# Open output files to write
try:
    
    # Define groups
    groups = [f'{d}_{ct_format(ct)}'
              for d in donors
              for ct in cell_types]

    out_files = open_multiple_files(groups, 'a')


    # Parse each sample file
    for C in SAMPLE_CODES:

        sample = regex_get(C, pattern='sample_name')
        file_path = os.path.join(DATA_PATH, "outputs", C, C, "outs/atac_fragments.tsv.gz")
        
        with gzip.open(file_path, 'rt') as fragments_file:
            
            for line in fragments_file:
                
                if not line.strip().startswith('#'):
                
                    cell_data = line.strip().split('\t')
                    ann_index = f'{cell_data[3]}_{sample}'
                    # print(f'current cell: {ann_index}')

                    if ann_index in ann_idxs:

                        donor, cell_type = cells_coldata.loc[ann_index][['donor', 'cell_type']]
                        # print(f'Passed QC: donor={donor}, cell_type={cell_type}')
                        # print(f'Writing to: {donor}_{cell_type}')
                        out_files[f'{donor}_{cell_type}'].write(line)
                

finally:
    if 'out_files' in locals():
        close_multiple_files(out_files)