import os
import re
import anndata as ad

def regex_get(text, pattern=None):
    
    pattern_map = {'sample_id': r'sSL\d*[A-Za-z]?',
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


def get_anndata_obs_essentials(anndata_path,
                               cols: list = ['sample', 'donor', 'donor_id', 'cell_type'],
                               cell_type_col: str = 'cell_type'):
    """ 
    Extracts and formats obs from anndata.

    Args:
        anndata_path (h5ad): Path to anndata file.
        cols (list[str]): List of columns to extract from the anndata.
                          If column not found it will be skipped unless extractabable from the index name such as sample and barcode.
        cell_type_col (str, opt): Column containing the annotation feature AKA the cell-type feature.
                                  It will be formatted as categorical and saved as `cell_type`.
    
    Return:
        list:
            obs_ess: Formatted obs pandas dataframe.
            ann_idxs: Indexes of obs.
            cell_types (opt): If `cell_type` argument given. List of unique cell-types.
            samples (opt): If `samples` in anndata.obs.columns or extractable from the index. List of unique donors.
            donors (opt): If `donors` in anndata.obs.columns or extractable from the index. List of unique donors.
        
    Example:
        >>> get_anndata_obs_essentials("path_to_ad.h5ad", cell_type_col="cts_predicted")
        [obs_ess, ann_idxs, cell_types, samples, donors]
        
    Raise:
        
    Note:
        Newer than `get_anndata_coldata`
    """
   
    adata = ad.read_h5ad(anndata_path, backed='r').obs
    
    cols_exist = cols.copy()
    
    for i in cols_exist:
        
        if i not in adata.columns.tolist():
            
            print(f"get_anndata_obs_essentials:" \
                  f"Column `{i}` not found in anndata.obs. SKIPPED")
            cols_exist.remove(i)
        
    obs_ess = adata[cols_exist].copy()
    
    # Format 'cell_type' column
    if cell_type_col:
        
        obs_ess['cell_type'] = obs_ess[cell_type_col].apply(ct_format).astype('category')
        
        if cell_type_col != 'cell_type':
            obs_ess = obs_ess.drop(cell_type_col, axis=1)
    
    # Add 'barcode' and 'sample' columns
    if ('barcode' in cols) and (not 'barcode' in obs_ess.columns.tolist()):
        
        print(f"get_anndata_obs_essentials:")
        print("\t`barcode` column not found, attempting to extract from index")
        obs_ess['barcode'] = [regex_get(i, 'barcode') for i in obs_ess.index.tolist()]
        print("\tExtracted.")
        
    if ('sample' in cols) and (not 'sample' in obs_ess.columns.tolist()):
        
        print(f"get_anndata_obs_essentials:")
        print("\t`sample` column not found, attempting to extract from index")
        obs_ess['sample'] = [regex_get(i, "sample_id") for i in obs_ess.index.tolist()]
        print("\tExtracted.")

    # Some extracted variables
    ann_idxs = obs_ess.index.tolist()

    # Return
    
    out_list = [obs_ess, ann_idxs]
    
    if 'cell_type' in obs_ess.columns.tolist():
        
        cell_types = tuple(set(obs_ess['cell_type']))
        cell_types = sorted(cell_types, key=str.upper)
        
        out_list.append(cell_types)
        
            
    if 'sample' in obs_ess.columns.tolist():

        samples = tuple(set(obs_ess['sample']))
        samples = sorted(samples, key=str.upper)
        
        out_list.append(samples)


    if 'donor' in obs_ess.columns.tolist():
        
        donors = tuple(set(obs_ess['donor']))
        donors = sorted(donors, key=str.upper)
        
        out_list.append(donors)

    return out_list


selected_cts = ['Glia', 'UL-EN', 'Midbrain-EN']
selected_donors = ['pelm', 'zoxy', 'ualf', 'melw']
