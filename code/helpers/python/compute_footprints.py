import os
import argparse
from typing import List, Tuple, Callable

import numpy as np
from scipy.sparse import csr_matrix
import anndata as ad

from scipy.spatial.distance import jensenshannon as js_distance
import pyBigWig

import utils

def parse_args():

    parser=argparse.ArgumentParser(description="Calculate the footprints of coverage files pertaining a an experimental condition.")
    parser.add_argument('-c', '--covs-dir', type=str, required=True, help='Coverage files directory. Should only contain coverage files in bigwig format [.bw] based on bedgraph format.')
    parser.add_argument('-o', '--out-adata', type=str, required=True, help='Output footprint file [.adata] to store the footprint matrix in.')
    parser.add_argument('-p', '--peaks-file', type=str, required=True, help='tsv file with peak ids as first column.')
    parser.add_argument('-H', '--header', action='store_true', help='Skip the header line of the peaks tsv file.')
    parser.add_argument('-a', '--algorithm', choices=['js_div', 'js_divergence', 'jensen_shannon_divergence'], default='jensen_shannon_divergence', type=str, help='Footprint calculation algorithm.')
    parser.add_argument('-A', '--obs-annotation', type=str, required=False, default=None, help='Format: "str1,str2". Annotate the observations with an .obs column named after the first string, and populate the column with the second string.')

    args = parser.parse_args()

    return args


def extract_donor_name(file_name: str):
    
    return file_name.split('_')[0]


def compute_footprints(
        coverage_files_dir: str,
        out_adata: str,
        peak_ids: List[str],
        algorithm: str = 'js_divergence',
        func_for_donor_naming: Callable[[List[str]], List[str]] = lambda x: os.path.splitext(x)[0],
        obs_annotation: Tuple[str, str] = None
        ):
    """
    Computes the footprint for all coverage files within a directory for a set of peaks/regions.

    Input:
        - coverage_files_dir: Directory of coverage files to take into account. Only considers coverage files within the directory which are in bigwig (.bw) format and they must be created from bedgraph files.
        - out_adata (str): Path of output anndata file.
        - peak_ids (List(str)): Peak ids to calculate the footprints on. Format: contig:start:end:length
        - algorithm (str, opt): Currently not implemented
        - func_for_donor_naming (func, opt): How to extract donor name from file names within coverage_files_dir
        - obs_annotation (Tuple(str, str)): Annotate the observations with an .obs column named after the first list element and populate the column with the secondn list element.

    TODOs:
        - 01: use multi-dimensional arrays in very grand scheme for js_distance calculation
        - scipy bug: https://github.com/scipy/scipy/pull/20786
            --> If working with low insertion counds, errors can occur. In this case and if this wasn't corrected in scipy:>1.13.1, solve by installing from the above mentioned branch or manually fixing the code.
    """
    js_approach = 'norm_total_no_all0' # absolute, norm_total_no_all0, norm_total_pseudocount, norm_min_max
        # absolute: no normalization. Coverages of all zeros, give them a small uniform nr.
        # norm_total_no_all0: normalize to total counts. Coverages of all zeros give them a small uniform nr.
        # norm_total_pseudocount:
        # norm_min_max: https://github.com/kundajelab/basepairmodels/blob/cf8e346e9df1bad9e55bd459041976b41207e6e5/basepairmodels/cli/metrics.py#L129

    condition = os.path.basename(coverage_files_dir)
    file_names = utils.list_files_and_links(coverage_files_dir, '.bw')
    donor_map = {func_for_donor_naming(fn): fn
              for fn in file_names} # {donor: donor_file}
    donor_map = {key: donor_map[key] for key in sorted(donor_map.keys())}
    donors_sorted = sorted(list(donor_map.keys())) # Reference indexes for np.array in case dict order fails: older python versions
    peak_ids = sorted(peak_ids, key=utils.sort_key_peaks)

    # Check peaks unique
    if len(set(peak_ids)) != len(peak_ids):
        raise ValueError(f'Peaks are not unique!')

    peaks = {}  # peak: ['chr', 'start', 'end']

    # Get peaks data
    for index, p in enumerate(peak_ids):

        utils.check_peak_format(p)

        p_split = p.split(':')
        peaks[p] = [str(p_split[0]), int(p_split[1]) - 1, int(p_split[2])]  # Indexing basis correction: 1-base fully-closed to 0-based half-open

    # Init vars
    arr = np.empty([len(donor_map), len(peaks.keys())])
    covs = {p: {} for p in peaks.keys()}

    # Populate
    for d, bw_name in donor_map.items():

        bw_file = f'{coverage_files_dir}/{bw_name}'

        # Populated bw file
        if not os.stat(bw_file).st_size == 0:

            bw = pyBigWig.open(bw_file)

            for p in peak_ids:

                # If peak in contigs
                if peaks[p][0] in set(bw.chroms().keys()):

                    cov = bw.values(*peaks[p][0:3])
                    cov = np.nan_to_num(cov, 0)

                    if js_approach =='norm_total_no_all0':

                        if np.sum(cov) == 0:

                            cov = np.full(len(cov), 1e-9)

                    covs[p][d] = cov
    
                else:

                    covs[p][d] = np.full(peaks[p][2] - peaks[p][1], 1e-9)

            bw.close()

        # Empty bw file
        else:

            for p in peak_ids:

                covs[p][d] = np.full(peaks[p][2] - peaks[p][1], 1e-9)


    for p_idx, p in enumerate(peak_ids):

        cov_arr = np.vstack(list(covs[p].values()))
        cov_mean = np.mean(cov_arr, axis=0)

        for d_idx, d in enumerate(donors_sorted):


            if js_approach == 'norm_total_pseudocount':

                pseudocount = 0.001

                p = (covs[p][d] + pseudocount) / (np.sum(covs[p][d]) + pseudocount)
                q = (cov_mean + pseudocount) / (np.sum(cov_mean) + pseudocount)

                js_dist = js_distance(p, q, base = 2)
                js_div = js_dist ** 2

                

            elif js_approach == 'norm_total_no_all0':

                js_dist = js_distance(cov_mean, covs[p][d], base = 2)
                #print(f'JSd {p}: d({d})={covs[p][d]} m={cov_mean} = {js_dist}')
                js_div = js_dist ** 2
                arr[d_idx, p_idx] = js_div
                #print(f'JSD {p} = {js_div}')

            else:

                raise ValueError(f'JS approach <{js_approach}> not implemented yet')

    # Create anndata
    adata = ad.AnnData(csr_matrix(arr))
    adata.obs_names = donors_sorted
    adata.var_names = peak_ids

    # Annotate
    adata.obs['donor'] = adata.obs.index

    if obs_annotation:
        adata.obs[obs_annotation[0]] = obs_annotation[1]

    # Save anndata
    utils.create_dir(out_adata, parent=True)
    adata.write(out_adata, compression='gzip')


def main():

    args = parse_args()

    peak_ids = utils.read_rownames(args.peaks_file, header=args.header)

    if args.obs_annotation:

        obs_annotation = args.obs_annotation.split(',')[0:2]

    else:

        obs_annotation = None
    

    compute_footprints(coverage_files_dir = args.covs_dir,
                       out_adata = args.out_adata,
                       peak_ids = peak_ids,
                       algorithm = args.algorithm,
                       func_for_donor_naming = extract_donor_name,
                       obs_annotation = obs_annotation)


if __name__ == '__main__':

    main()
