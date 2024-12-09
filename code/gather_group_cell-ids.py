def gather_group_cellids(anndata_path,
                         group_col_id: str,
                         out_dir: str,
                         ct_map_path: str=None):
    """
    
    Extract a categorical annotation (e.g. cell-type) from an anndata object and create a file containing all cell-ids pertinent to the donor-category. Optionally combine several categories after a scheme.
    
    Args:
        anndata (h5ad): Anndata object with cells x features and a column in the .obs with the group annotation.
        ann_group_col (str): Column in the anndata.obs df with the group annotation.
        out_dir (str): Path to directory where the group files containing the cell ids will be created.
                       Naming format for new files: <group>.txt
        dataset_id (str): Query dataset ID.
        ct_map_path (str, opt): Path to json with scheme to group annotation groups by.
                      #ID of scheme to group annotation groups by.
                      #It will search in the `/config/cell-type_groupings/` folder for a json file named <ct_map_id>.json.
                      Json structure: key = new group name --> array of old group names to be gathered.
                      If not provided, the groupings in the anndata.obs['<group_ann_col>'] will be used.

    
    Returns:
        nothing: creates group files containing cell-ids in target directory.
        
    Raises:
    
    Examples:
        >>> gather_group_cellids(pbmcs.h5ad, "cell_type", "/home/user/documents/project/out", "only-neurons.json")
        
    Note:

    """

    import os
    import json
    import sys
    from helpers.helpers import ct_format, get_anndata_obs_essentials

    if not os.path.isdir('out_dir'):
        os.makedirs('out_dir')
    
    ## Gather relevant data
        
    # formatted ad.obs
    cells_metadata, _, cell_types, _, donors = get_anndata_obs_essentials(anndata_path, cell_type_col='cell_type')
    
    # Get cell-type map
    if ct_map_path:
        
        ct_map_id = '_' + os.path.basename(ct_map_path).split('.')[0]

        with open(ct_map_path, 'r') as json_file:
            ct_map= json.load(json_file)

        reverse_ct_map = {ct_format(ct): ct_format(group)
                          for group, cts in ct_map.items() for ct in cts}

        new_cts = [ct_format(ct) for ct in list(ct_map.keys())]
        
    else:
        
        ct_map_id = ''
        ct_map = None
        new_cts = [ct_format(ct) for ct in cell_types]
        reverse_ct_map = {i: i for i in new_cts}


    # Extend ad.obs
    group_col_id = f'cell_type{ct_map_id}'
    cells_metadata[group_col_id] = cells_metadata['cell_type'].map(reverse_ct_map)
    
    
    ## Gather donor-group specific cells and output
    
    for d in donors:
        
        for cat in new_cts:

            donor_cat = f'{d}_{cat}'
            out_file = f'{out_dir}/{donor_cat}.txt'
            
            # Avoid repeated run concatenation issues
            if os.path.isfile(out_file):
                os.remove(out_file)
            
            d_group_cell_ids = cells_metadata[(cells_metadata['donor'] == d) & (cells_metadata[group_col_id] == cat)].index.tolist()
            
            with open(out_file, 'w') as tsv:
                
                for cid in d_group_cell_ids:
                    
                    tsv.write(f'{cid}\n')
                    
if __name__ == "__main__":
    
    import os
    import sys
    
    PROJ_ROOT = ".."
    os.chdir(PROJ_ROOT)

    #DATASET_ID = "hca_brain-organoids"
    DATASET_ID = "toy_atac-seq"

    #CT_MAP_PATH = f"config/cell-type_groupings/{DATASET_ID}/{CT_MAP_KEY}.json"
    CT_MAP_PATH = f"/home/fichtner/data/toy_datasets/atac-seq_fragments_simple/ct_map.json"
    #ATAC_PEAKS_PATH = f"data/datasets/{DATASET_ID}/atac-seq/chromatin_accessibility/peak-matrix_rna-qc-cells_norm-reads-in-tss.h5ad"
    ATAC_PEAKS_PATH = f"/home/fichtner/data/toy_datasets/atac-seq_fragments_simple/peaks_raw.h5ad"
    #OUT_DIR = f"data/intermediate-data/datasets/{DATASET_ID}/annotations_general/group_cell-ids/{CT_MAP_KEY}"
    OUT_DIR = f"/home/fichtner/data/toy_datasets/atac-seq_fragments_simple/gci"
    
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)

    gather_group_cellids(ATAC_PEAKS_PATH,
                         group_col_id="cell_type",
                         out_dir=OUT_DIR,
                         ct_map_path=CT_MAP_PATH)
