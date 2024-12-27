"""
Process the pre-annotated footprint anndatas:
    - calculate descriptive statistics: peak-level mean, std, var, min, max, 25%, 50%, 75%
    - calculate PCA, UMAP and clustering with scanpy

Requires:
    - footprint anndatas
        - adata.var['n_cells'] annotation
"""

import os
import sys
import argparse
from typing import Tuple
import warnings
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import anndata as ad
import scanpy as sc

sys.path.append(os.path.dirname(__file__))
from anndata_utils import check_anndata, subset_common_donors
from utils import create_dir

sys.path.append(os.getcwd() + '/code')
from glob_vars import GENOTYPES_TSV, GENOTYPE_PCS_TSV



def parse_args():

    parser = argparse.ArgumentParser("Description: Process anndata by calculating the PCA, UMAP and clustering.")

    parser.add_argument('-a', '--anndata', type=str, required=True,
                        help='Path to anndata [.h5ad] to annotate.')
    parser.add_argument('-o', '--output-anndata', type=str, required=True,
                        help='Path to output anndata')
    parser.add_argument('-m', '--min-cells', type=int, required=True, default=0,
                        help='Minimum nr of cells per donor to keep. Requires anndata to have .obs[n_cells]')
    parser.add_argument('-p', '--min-peak-counts', type=int, required=False, default=0,
                        help='Minimum nr of counts per peak. Extracted from .var[peak_counts_total]')

    return parser.parse_args()


def plot_peak_count_histogram(adata, out_file: str):
    '''
    Plot histogram of total number of counts for a peak (across donors, peak-level)

    Notes:
        - Uses .var['peak_counts_total'] for count extractino
        - Uses .var['length'] to determine the range of the restricted range plot
    '''

    data = adata.var['peak_counts_total']
    create_dir(out_file)


    ## Full range

    range = (np.min(data), np.max(data))

    plt.figure(figsize=(11, 6))
    plt.hist(data, range=range, bins=np.linspace(range[0], range[1] , int(np.ceil(np.subtract(*range[::-1]))) + 1))
    plt.xlabel('Insertions', fontsize=14)
    plt.ylabel('#', fontsize=14)

    plt.savefig(out_file, dpi=300, bbox_inches='tight')


    ## Restricted range

    if 'length' not in adata.var.columns:

        warnings.warn('adata.var did not contain column <length>. Skipping restricted range plot.')

    if adata.var['length'].nunique() != 1:

        warnings.warn('adata.var[length] had heterogeneous values. Skipping restricted range plot.')


    length = int(adata.var['length'][0])
    range = (np.min(data), np.subtract(length, 0.000000001))

    plt.figure(figsize=(11, 6))
    plt.hist(data, range=range, bins=np.linspace(range[0], range[1] , int(np.ceil(np.subtract(*range[::-1]))) + 1))
    plt.xlabel('Insertions', fontsize=14)
    plt.ylabel('#', fontsize=14)

    out_file_name = os.path.basename(out_file)
    out_file_restricted = os.path.join(os.path.dirname(out_file), f'{os.path.splitext(out_file_name)[0]}_resricted-range{os.path.splitext(out_file_name)[1]}')
    plt.savefig(out_file_restricted, dpi=300, bbox_inches='tight')


def descriptive_stats_adata(adata, dim: Tuple[str]):
    """
    Calculate var(peak)-level mean, std, var, min, max, 25%, 50%, 75%

    Input:
        - adata
        - dim (tuple) {'var', 'obs'}: Choose for which dimension(s) to compute the descriptive statistics.
    """

    adata_processed = adata.copy()

    if 'var' in dim:

        stats = adata_processed.to_df().describe().T
        stats['var'] = stats['std'] ** 2
        stats = stats[['mean', 'std', 'var', 'min', 'max', '25%', '50%', '75%']]

        assert adata_processed.var.index.equals(stats.index), "Index contents or order aren't equal!"

        adata_processed.var = adata_processed.var.join(stats, how='left')

    if 'obs' in dim:

        raise ValueError("Not implemented yet")


    return adata_processed



def highly_variable_peaks_adata(adata, top_percent: int = 0.1, n_top_peaks: int = None):
    """
    Calculate the highly variable peaks.

    TODOS:
        - [ ] implement absolute values
    """
    print('HVPs: Absolute values not yet implemented')


    # Checks
    assert (top_percent or n_top_peaks ) and not (top_percent and n_top_peaks)


    adata_processed = adata.copy()


    # Prep

    if not 'std' in adata.var.columns.to_list():

        adata_processed = descriptive_stats_adata(adata_processed, dim=['var'])

    var = adata.var


    # Calculate hvps and their ranks
    if top_percent:

        top_peaks = var['std'].nlargest(int(len(var) * top_percent)).index

    elif n_top_peaks:
    
        top_peaks = var['std'].sort_values(ascending=False)[0:n_top_peaks].index

    rank = var['std'].rank(ascending=False, method='first').astype(int)


    # Annotate
    adata_processed.var['highly_variable_std'] = var.index.isin(set(top_peaks))
    adata_processed.var['std_rank'] = var.index.map(rank.to_dict())


    return adata_processed



def process_adata(adata):

    adata_processed = adata.copy()

    # Checks
    assert adata.shape[0] > 1, 'Anndata has 0 observations!'
    assert adata.shape[1] > 1, 'Anndata has 0 features!'


    # Processing
    adata_processed = highly_variable_peaks_adata(adata_processed)
    adata_processed = sc.tl.pca(adata_processed, copy=True)
    adata_processed = sc.pp.neighbors(adata_processed, copy=True)
    adata_processed = sc.tl.leiden(adata_processed, flavor='igraph', n_iterations=2, copy=True)
    adata_processed = sc.tl.umap(adata_processed, copy=True)

    return adata_processed


def main():

    args = parse_args()
    adata = ad.read_h5ad(args.anndata)

    ## Checks

    assert args.anndata.endswith('.h5ad')
    check_anndata(adata, min_obs=20, light=True)


    ## Pre-process

    # Subset
    adata_processed = adata[adata.obs['n_cells'] > args.min_cells].copy() # Filter min n_cells per donor,
                                                                          # should be done at footprint creation idealy
    adata_processed = subset_common_donors(adata_processed, GENOTYPES_TSV, GENOTYPE_PCS_TSV) # Remove donors not
    adata_processed = adata_processed[:, adata_processed.var['peak_counts_total'] > args.min_peak_counts].copy()

    # Peak counts histogram
    peak_counts_png = os.path.join(os.path.dirname(args.output_anndata), 'metadata', 'peak_counts_total_histogram.png')
    plot_peak_count_histogram(adata_processed, peak_counts_png)
                                                      
    adata_processed = descriptive_stats_adata(adata_processed, dim=['var'])

    # Checkpoint save
    adata_processed.write(args.output_anndata, compression='gzip')

    # Process
    check_anndata(adata_processed, min_obs=0, light=False) # Don't permit NaNs in X
    adata_processed = process_adata(adata_processed)
        
    # Save
    adata_processed.write(args.output_anndata, compression='gzip')


if __name__ == '__main__':

    main()
    os._exit(0) # For some reason the script doesn't end by itself
