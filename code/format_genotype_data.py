"""
Format all genotype relevant data prior to the pipeline.

Genotype matrix: correct snp label format and remove indels

Input:
    - GENOTYPES_TSV

Output:
    - GENOTYPES_PROCESSED_TSV
"""

from typing import Callable
import pandas as pd
from glob_vars import GENOTYPES_TSV, GENOTYPES_PROCESSED_TSV


### Functions ###

def index_filter_no_indels(label: str):
    
    ref_allele, alt_allele = label.split('_')[2:4]

    return len(ref_allele) == 1 and len(alt_allele) == 1


### Script ###

# Genotype matrix

gt = pd.read_csv(GENOTYPES_TSV, sep='\t', header=0, index_col=0)


# Format index
gt = gt.rename(index={i: 'chr' + i for i in gt.index})

# Remove indels
gt = gt[gt.index.to_series().apply(index_filter_no_indels)]

# Remove duplicates
gt = gt.loc[gt.index.drop_duplicates(), :]
gt = gt.groupby(level=0).first() # Some artifact makes the previous line still have duplicates


# Save
gt.to_csv(GENOTYPES_PROCESSED_TSV, sep='\t')
