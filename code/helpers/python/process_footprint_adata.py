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
import pandas as pd
import anndata as ad
import scanpy as sc

sys.path.append(os.path.dirname(__file__))
from anndata_utils import check_anndata, subset_common_donors

sys.path.append(os.getcwd() + '/code')
from glob_vars import GENOTYPES_TSV, GENOTYPE_PCS_TSV



def parse_args():

    parser = argparse.ArgumentParser("Description: Process anndata by calculating the PCA, UMAP and clustering.")

    parser.add_argument('-a', '--anndata', type=str, required=True,
                        help='Path to anndata [.h5ad] to annotate.')
    parser.add_argument('-o', '--output-anndata', type=str, required=True,
                        help='Path to output anndata')
    parser.add_argument('-m', '--min-cells', type=int, required=True,
                        help='Minimum nr of cells per donor to keep. Requires anndata to have .obs[n_cells]')

    return parser.parse_args()



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
    adata_processed = descriptive_stats_adata(adata_processed, dim=['var'])
    adata_processed = highly_variable_peaks_adata(adata_processed)
    adata_processed = sc.tl.pca(adata_processed, copy=True)
    adata_processed = sc.pp.neighbors(adata_processed, copy=True)
    adata_processed = sc.tl.leiden(adata_processed, flavor='igraph', n_iterations=2, copy=True)
    adata_processed = sc.tl.umap(adata_processed, copy=True)

    return adata_processed



def main():

    args = parse_args()
    adata = ad.read_h5ad(args.anndata)

    # Checks
    assert args.anndata.endswith('.h5ad')
    check_anndata(adata, min_obs=20)


    # Process adata
    adata_processed = adata[adata.obs['n_cells'] > args.min_cells].copy() # Filter min n_cells per donor,
                                                                          # should be done at footprint creation idealy
    adata_processed = subset_common_donors(adata_processed, GENOTYPES_TSV, GENOTYPE_PCS_TSV) # Remove donors not
                                                                                             # found in genotype data
    adata_processed = process_adata(adata_processed)
        
    # Save
    adata_processed.write(args.output_anndata, compression='gzip')

    sys.exit(0) # For some reason the script doesn't end by itself


if __name__ == '__main__':
    
    main()
