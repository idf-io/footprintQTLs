from typing import Tuple, List, Dict, Callable
import numpy as np
import pandas as pd
import anndata as ad

import os
import sys
sys.path.append(os.path.dirname(__file__))
from pandas_utils import aggregate_df_by_columns



def check_anndata(adata, min_obs: int = 20, obs_criteria_kwargs: List[Dict[str, Callable]] = []):
    """
    Check anndata integrity for criteria:
	- Unique obs and vars
	- No NaNa
        - min nr obs
        - custom .obs criteria
    """

    # Data integrity
    assert adata.obs.index.nunique() == len(adata.obs)               # Obs unique
    assert adata.var.index.nunique() == len(adata.var)               # Var unique
    assert not np.any(np.isnan(adata.X.data))                        # Non NaNs

    # obs criteria
    for c in obs_criteria_kwargs:

        assert len(c.keys()) == 2
        assert 'col' in c.keys()
        assert 'func' in c.keys()
        
        col = c['col']
        func = c['func']

        assert func(adata.obs[col]), f'obs_criteria_kwargs check failed: col={col}, func={func}'

    # Min
    if min_obs:
        assert adata.shape[0] > min_obs, "Not enough obs/donors!" # Min donors


    assert not any(adata.obs['donor_id'].isna())                     # obs['donor_id'] complete
    assert not adata.obs['donor_id'].nunique == len(adata.obs)       # obs['donor_id'] unique


def subset_common_donors(adata, genotypes_tsv: str, genotype_pcs_tsv: str):
    """
    Subset to common donors in all data sources (AKA also in  genotype data)

    Notes: adata donor information extracted from .obs['donor_id']
    """

    adata_subset = adata.copy()

    # adata donors
    donor_ids_adata = set(adata_subset.obs['donor_id'].to_list())

    # genotype donors
    genotype = pd.read_csv(genotypes_tsv, sep='\t', header=0, index_col=0)
    donor_ids_genotype = set(genotype.columns.to_list())

    # genotype pcs donors
    genotype_pcs = pd.read_csv(genotype_pcs_tsv, sep="\t", index_col=0)
    donor_ids_genotype_pcs = set(genotype_pcs.index.to_list())

    # common donors
    donor_ids = donor_ids_adata & donor_ids_genotype & donor_ids_genotype_pcs


    # Subset adata
    adata_subset = adata_subset[adata_subset.obs['donor_id'].isin(donor_ids), :].copy()

    del genotype
    del genotype_pcs

    return adata_subset



def port_obs_adata(adata, adata_ref, suffix: str = '_ref_ad', pos_filter_ref: Tuple[str] = None, neg_filter_ref: Tuple[str] = None, kwargs: Dict[str, str] = None):
    """
    Expand anndata .obs columns with the columns of a second anndata of equal or subset of obs dimensions.

    Input:
        - adata (anndata): anndata to expand its obs
        - adata_ref (anndata): reference anndata to expand the first anndata .obs columns with
        - suffix (str): if both anndatas have a common column, suffix the column from the ref_anndata
        - pos_filter_ref (list [str], optional): filter reference var by selecting columns.
        - neg_filter_ref (list [str], optional): filter reference var by discarding columns.
        - grouping_col (str): Group and aggregate the reference anndata.obs by a column.
            Aggregation strategy: mean for numerical, for categorical only keep if homogeneous values per group.
            Useful when reference anndata is e.g. pre-pseudobulked.
        - kwargs (dict):
            - 'grouping_col': categorical column in reference anndata.obs to group obs/cells according to.
            - 'filter_col': optional, if you want to filter the adata, indicate an .obs column and ...
            - 'filter_key': optional, key to filter for.
            - 'obs_map_cal': which column to use as the obs map.
                E.g. donor if you are averging over single-cells. Can also indicate 'index' for the .obs_names.

    Output:
        - adata with .obs columns expanded
    """

    obs_ref = adata_ref.obs

    # Format kwargs

    kwarg_options = ['grouping_col', 'filter_col', 'filter_key', 'obs_map_col']

    for k in kwargs.keys():
        assert k in kwarg_options

    assert len(kwargs) in [0,2,4]


    # Filter cols

    if pos_filter_ref:
        obs_ref = obs_ref[pos_filter_ref]

    if neg_filter_ref:
        obs_ref = obs_ref.drop(columns=neg_filter_ref)


    # Optionally: Filter by col

    if 'filter_col' in kwarg_options and 'filter_key' in kwarg_options:

        obs_ref = obs_ref[obs_ref[kwargs['filter_col']] == kwargs['filter_key']]
        assert obs_ref[kwargs['filter_col']].nunique() == 1


    # Optionally: group by col

    if 'grouping_col' in kwarg_options and 'obs_map_col' in kwarg_options:

        obs_ref = aggregate_df_by_columns(obs_ref, kwargs['grouping_col'])

        if kwargs['obs_map_col'] != 'index':

            obs_ref = obs_ref.set_index(kwargs['obs_map_col'], drop=True)


    # Annotate
    obs = adata.obs
    obs_exp = obs.join(obs_ref, rsuffix=suffix, how='left')

    adata_annotated = adata.copy()
    adata_annotated.obs = obs_exp

    return adata_annotated


    
def port_var_adata(adata, adata_ref, suffix = '_ref_ad', key: str = None, pos_filter_ref: Tuple[str] = None, neg_filter_ref: Tuple[str] = None):
    """
    Expand anndata .var columns with the columns of a second anndata of equal or subset of var dimensions.

    Input:
        - adata (anndata): anndata to expand its var
        - adata_ref (anndata): reference anndata to expand the first anndata .var columns with
        - key (str, optional): if the reference anndata has the var index to join on in a column instead of index.
        - pos_filter_ref (list [str], optional): filter reference var by selecting columns.
        - neg_filter_ref (list [str], optional): filter reference var by discarding columns.

    Output:
        - adata with .var columns expanded
    """

    var = adata.var
    var_ref = adata_ref.var

    if pos_filter_ref:
        var_ref = var_ref[pos_filter_ref]

    if neg_filter_ref:
        var_ref = var_ref.drop(columns=neg_filter_ref)

    if key:
        var_ref = var_ref.set_index(key, drop=True)

    var_exp = var.join(var_ref, rsuffix = suffix, how='left')

    adata_annotated = adata.copy()
    adata_annotated.var = var_exp

    return adata_annotated
