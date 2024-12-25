"""
Pre-annotate the anndata object with
    - donor, donor_id
 	- nr. of cells per donor
 	- nr. of insertions per donor
 	- nr. of fragments per donor
 	- locations of peaks
 	- metadata from reference anndata

 Requires:
 	- anndata
 	- path with donor fragment files, to extract:
        - donor, donor_id
 		- n_insertions per donor
 		- n_fragmetns per donor
 		- n_cells per donor: 1st option to extract this data
 		- path format: "%d_myDS.tsv.gz" where %d represents the donors that will be matched t the anndata obs
 	- reference anndata with other .obs and .vars annotations, to extract:
 		- n_cells per donor, 2nd option to extract this data
        - var and obs metadata

 Assumptions:
 	- anndata object contains an .obs column with the ct that each observation(donor) belongs to.
 	- in requirement 2: donor-ct grouped fragment files, all files which match *ct* (1 lvl deep) are used to match cell-type from footprint anndata and grouped fragment-file anndata
"""

import os
import sys
import argparse
import re
import subprocess
from typing import Dict, Tuple
import anndata as ad
import pandas as pd

sys.path.append(os.path.dirname(__file__))
from utils import list_files_and_links, ct_format
from pandas_utils import aggregate_df_by_columns
from anndata_utils import port_obs_adata, port_var_adata
from fragment_files_utils import annotate_adata_counts_from_fragment_files


def parse_args():

    parser=argparse.ArgumentParser("Description: Annotate anndata with n_cells and n_insertions per donor(-cell_type)")

    parser.add_argument('-a', '--anndata', type=str, required=True,
                        help="Path to anndata file to annotate [.h5ad].")
    parser.add_argument('-o', '--output', type=str, required=True,
                        help="Path to output anndata file.")
    parser.add_argument('-f', '--fragment-files', type=str, required=True,
                        help="""Path to fragment files used to create footprints anndata.
                                Use %s as a placeholder of the donor; e.g. <2024_%d_myDS.tsv.gz>. 
                                The matched donors within the directory will be mapped to the donors in the anndata. 
                                Used to count the nr. of insertions, fragments and optionally cells per donor.""")
    parser.add_argument('-r', '--reference-adata', type=str, required=True, 
                        help="""Path to reference pre-pseudobulked anndata object [.h5ad] which 
                                will be used to (1) port some .vars and .obs annotations 
                                and optionally (2) count the cell nr. per donor.""")
    parser.add_argument('-u', '--use-ref-anndata', action='store_true', 
                        help="Use reference anndata to extract cell counts per donor.")
    #parser.add_argument('-F', '--filter', type=str, required=False, default=False, 
                        #help="""Indicate anndata.obs column and key pairs to optionally filter the reference anndata 
                                #(separate key:value by colon and key:value pairs with comma
                                #e.g. cell_type:DL-EN, QC:pass""")
    parser.add_argument('-k', '--kwargs', type=str, required=False, 
                        help="""Kwargs when -u is used. 
                                'grouping_col': categorical column name to group and count cells according to. 
                                'filter_col': .obs column to filter before counting obs in a group. 
                                'filter_key': .obs[filter_col] key to filter for. 
                                Example: <grouping_col:donor,filter_col:cell_type,filter_key:DL-EN""")

    args = parser.parse_args()

    return args


def annotate_ncells_adata_ref_adata(adata: str, adata_ref, kwargs: Dict[str, str] = None):
    """
    Annotate the anndata (e.g. donors) with the nr of cells per donor of the pre-pseudobulk single-cell dataset.
    Use a reference pre-pseudobulked dataset.

    Input:
        - adata (anndata)
        - adata_ref_path (str): Anndata with pre-pseudobulked obervations.
        - kwargs (dict):
            - 'grouping_col': categorical column in reference anndata.obs  to group and count cells according to.
            - 'filter_col': optional, if you want to filter the adata, indicate an .obs column and ...
            - 'filter_key': optional, key to filter for.
            - 'obs_map_cal': which column to use as the obs map.
                E.g. donor if you are mapping donors to their cell nr. Can also indicate 'index' for the .obs_names.

    Output:
        - adata with columns adata.obs['n_cells']
    """

    adata_ref = adata_ref.copy()

    # Format kwargs

    kwarg_options = ['grouping_col', 'filter_col', 'filter_key', 'obs_map_col']

    for k in kwargs.keys():
        assert k in kwarg_options

    assert len(kwargs) in [2,4]


    # Get ncells per group

    if 'filter_col' in kwarg_options and 'filter_key' in kwarg_options:

        adata_ref = adata_ref[adata_ref.obs[kwargs['filter_col']] == kwargs['filter_key']]
        assert adata_ref.obs[kwargs['filter_col']].nunique() == 1

    donor_ncells_map = adata_ref.obs.groupby(kwargs['grouping_col']).size().to_dict()


    # Annotate

    adata_annotated = adata.copy() # Just in case anndata objects are mutable (change within function)

    if kwargs['obs_map_col'] == 'index':

        print(len(donor_ncells_map))
        print(donor_ncells_map)

        print(len(adata_annotated.obs))
        print(adata_annotated.obs)
        print(adata_annotated.obs.index)

        print(len(adata_annotated.obs.index.map(donor_ncells_map)))
        print(adata_annotated.obs.index.map(donor_ncells_map))

        adata_annotated.obs['n_cells'] = adata_annotated.obs.index.map(donor_ncells_map).fillna(0).astype(int)

    else:   

        adata_annotated.obs['n_cells'] = adata_annotated.obs[kwargs['obs_map_col']].map(donor_ncells_map).fillna(0).astype(int)


    return adata_annotated



def annotate_ncells_adata_frag_files(adata: str, fragment_files_obs_map: dict):
    """
    Annotate anndata with nr of cells per observation (e.g. donor) using obs-mapped fragment files.

    Input:
        - adata (anndata)
        - obs_fragment_files_map (dict): observation -> fragment file map to annotate the anndata by. Fragment files are bed formatted with at leas chr, start, end, cell_id columns.

    Output:
        - adata with columns adata.obs['n_cells']

    Assumptions:
        - cell_ids in fragment file column 4 must be cell-type unique
        - When using a fragment file: the file was the sole source of fragments for creating the footprints anndata.
    """

    # Checks
    for f in obs_fragment_files_map.values():
        assert f.endswith('.tsv') or f.endswith('.tsv.gz')


    n_cells_map = {} # obs: n_cells

    for obs, frag_file in obs_fragment_files_map.items():

        # Construct bash command

        if frag_file.endswith('.gz'):

            command_count_lines = 'gunzip -c \'{frag_file}\' | grep -v \'^#\' | cut -f4 | sort | uniq | wc -l'

        else:

            command_count_lines = 'grep -v \'^#\' \'{frag_file}\' | cut -f4 | sort | uniq | wc -l'


        # Count lines
        result = subprocess.run(command_count_lines, shell=True, capture_output=True, text=True, check=True)

        if result.stderr:
            raise RuntimeError(f'Command returned an error: {n_lines.stderr}')

        n_cells = int(result.stdout)

        n_cells_map[obs] = n_cells


    # Annotate anndata using maps
    adata_annotated = adata.copy() # Just in case anndata objects are mutable (change within function)
    adata_annotated.obs['n_cells'] = adata_annotated.obs.index.map(n_cells_map).fillna(0).astype(int)

    return adata_annotated



def annotate_peak_locs_adata(adata: str):
    """
    Annotate anndata by constructing .vars['peak_name', 'chr', 'start', 'end'] based on the adata peak names (var_names).

    Assumptions:
        - `adata.var_names` in format 'contig:start:end:length'
    """

    adata_annotated = adata.copy()
    adata_annotated.var['peak_name'] = adata_annotated.var.index
    adata_annotated.var['chr'] = adata_annotated.var['peak_name'].str.split(':').str[0]
    adata_annotated.var['start'] = adata_annotated.var['peak_name'].str.split(':').str[1]
    adata_annotated.var['end'] = adata_annotated.var['peak_name'].str.split(':').str[2]
    adata_annotated.var['length'] = adata_annotated.var['peak_name'].str.split(':').str[3]

    return adata_annotated
    



def main():

    args = parse_args()

    # checks
    assert args.anndata.endswith('.h5ad')

    # Make donor:frag-file map
    frag_files_dir = os.path.dirname(args.fragment_files)
    files = list_files_and_links(frag_files_dir, extension=".tsv.gz")
    frag_files_pattern1 = os.path.basename(args.fragment_files)
    frag_files_pattern2 = os.path.basename(args.fragment_files).replace('%d', '(.*)')


    donor_frag_file_map = {}

    for file in files:

        match = re.search(frag_files_pattern2, file)

        if match:

            donor = match.group(1)
            donor_frag_file_map[donor] = os.path.join(frag_files_dir, frag_files_pattern1.replace('%d', donor))


    # Format kwargs

    kwargs = {}

    for pair in args.kwargs.split(','):

        key, value = pair.split(':')

        key = key.strip()
        value = value.strip()

        kwargs[key] = value


    # Annotate adata
    adata = ad.read_h5ad(args.anndata)
    adata_ref = ad.read_h5ad(args.reference_adata)


    donor_id_map = adata_ref.obs[['donor', 'donor_id']].set_index('donor')['donor_id'].to_dict()
    adata.obs['donor_id'] = adata.obs['donor'].map(donor_id_map)


    adata_annotated = annotate_peak_locs_adata(adata)
    adata_annotated = annotate_adata_counts_from_fragment_files(adata_annotated, donor_frag_file_map, ['fragments', 'insertions'])


    if args.use_ref_anndata:

        adata_annotated = annotate_ncells_adata_ref_adata(adata_annotated, adata_ref, kwargs)

    else:

        adata_annotated = annotate_ncells_adata_frag_files(adata_annotated, donor_frag_file_map)


    adata_annotated = port_obs_adata(adata_annotated, adata_ref, neg_filter_ref=['leiden'], kwargs=kwargs)
    adata_annotated = port_var_adata(adata_annotated, adata_ref, key='peak_name', neg_filter_ref=['chr', 'start', 'end'])


    # Save
    adata_annotated.write(args.output, compression='gzip')
    adata.file.close()



if __name__ == '__main__':

    main()
