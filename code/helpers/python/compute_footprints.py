import os
import argparse
from typing import List, Tuple, Callable, Dict

import numpy as np
from scipy.sparse import csr_matrix
import anndata as ad

from scipy.spatial.distance import jensenshannon as js_distance
import pyBigWig

import utils


def parse_args():

    parser=argparse.ArgumentParser(description="Calculate the footprints of coverage files.")
    parser.add_argument('-c', '--covs-dir', type=str, required=True, help='Coverage files directory. Should only contain coverage files in bigwig format [.bw] based on bedgraph format.')
    parser.add_argument('-o', '--out-adata', type=str, required=True, help='Output footprint file [.adata] to store the footprint matrix in.')
    parser.add_argument('-p', '--peaks-tsv', type=str, required=True, help='tsv file with peak ids.')
    parser.add_argument('-a', '--algorithm', choices=['counts', 'js_divergence'], default='js_divergence', type=str, help='Footprint calculation algorithm.')
    parser.add_argument('-A', '--obs-annotation', type=str, required=False, default=None, help='Format: "str1,str2". Annotate the observations with an .obs column named after the first string, and populate the column with the second string.')

    args = parser.parse_args()

    return args


def compute_footprints(
        coverage_files_dir: str,
        out_adata: str,
        peaks: Dict[str, list],
        algorithm: str = 'js_divergence',
        func_for_donor_naming: Callable[[List[str]], List[str]] = lambda x: os.path.splitext(x)[0],
        obs_annotation: Tuple[str, str] = None
        ):
    """
    Computes the footprint for all coverage files within a directory for a set of peaks/regions.

    Input:
        - coverage_files_dir: Directory of coverage files to take into account. Only considers coverage files within the directory which are in bigwig (.bw) format and they must be created from bedgraph files.
        - out_adata (str): Path of output anndata file.
        - peaks (Dict(str, list)): Peaks to calculate the footprints on. Format: dict{peak_name: [contig,start,end]}
        - algorithm (str, opt): Algorithm to compute footprints with {counts, js_divergence, wasserstein_dist}.
        - func_for_donor_naming (func, opt): How to extract donor name from coverage file names.
        - obs_annotation (Tuple(str, str)): Annotate the observations with an .obs column named after the first list element and populate the column with the secondn list element.

    Output:
        - anndata with computed footprints (donors x peaks/regions)

    TODOs:
        - 01: use multi-dimensional arrays in very grand scheme for js_distance calculation
        - scipy bug: https://github.com/scipy/scipy/pull/20786
            --> If working with low insertion counds, errors can occur. In this case and if this wasn't corrected in scipy:>1.13.1, solve by installing from the above mentioned branch or manually fixing the code.
    """

    ## Get variables

    file_names = utils.list_files_and_links(coverage_files_dir, '.bw')
    donor_map = {func_for_donor_naming(file): file
                 for file in file_names} # {donor: donor_file}
    donor_map = {key: donor_map[key] for key in sorted(donor_map.keys())}

    # Reference indexes for np.array in case dict order fails: older python versions
    donors_sorted = sorted(list(donor_map.keys())) 
    peak_ids = peaks.keys()

    js_approach = 'all_zeros_to_uniform' # {all_zeros_to_uniform, norm_min_max}
        # norm_total_no_all0: normalize to total counts. Coverages of all zeros give them a small uniform nr.
        # norm_min_max: https://github.com/kundajelab/basepairmodels/blob/cf8e346e9df1bad9e55bd459041976b41207e6e5/basepairmodels/cli/metrics.py#L129
        # norm_total_pseudocount


    ## Checks

    if len(set(list(peaks.keys()))) != len(set(list(peaks.keys()))):
        raise ValueError(f'Peak ids are not unique!')

    algorithm_options = ('counts', 'js_divergence')
    assert algorithm in algorithm_options

    assert utils.are_sublists_unique(peaks.values()), 'Peak locations not unique!'

    js_approach_options = ('all_zeros_to_uniform', 'norm_min_max') # norm_total_pseudocount
    assert js_approach in js_approach_options, f'JS approach <{js_approach}> not implemented yet'


    ## Get coverages

    covs = {peak_id: {} for peak_id in peaks.keys()} # {peak_id: {donor: coverage}}

    for donor, donor_bw in donor_map.items():

        bw_file = f'{coverage_files_dir}/{donor_bw}'
        assert os.stat(bw_file).st_size != 0, f'Bigwig file empty! <{bw_file}>'

        bw = pyBigWig.open(bw_file)

        for peak_id in peak_ids:

            # If peak contig in bw contigs
            if peaks[peak_id][0] in set(bw.chroms().keys()):

                cov = bw.values(*peaks[peak_id][0:3])
                cov = np.nan_to_num(cov, 0)

                if algorithm == 'js_divergence' and js_approach == 'all_zeros_to_uniform':

                    if np.sum(cov) == 0:

                        cov = np.full(len(cov), 1e-9)

                covs[peak_id][donor] = cov

            # No counts in contig
            else:
                    covs[peak_id][donor] = np.full(peaks[peak_id][2] - peaks[peak_id][1], 1e-9)

        bw.close()


    ## Compute profiles

    arr_metric = np.empty([len(donor_map), len(peaks.keys())])
    arr_counts = np.empty([len(donor_map), len(peaks.keys())])

    for peak_idx, peak_id in enumerate(peak_ids):


        # Compute metric

        if algorithm == 'js_divergence':

            cov_arr = np.vstack(list(covs[peak_id].values()))
            cov_mean = np.mean(cov_arr, axis=0)

            for donor_idx, donor in enumerate(donors_sorted):

                if js_approach == 'all_zeros_to_uniform':

                    js_dist = js_distance(cov_mean, covs[peak_id][donor], base = 2)
                    js_div = js_dist ** 2
                    arr_metric[donor_idx, peak_idx] = js_div

                elif js_approach == 'norm_total_pseudocount':

                    raise ValueError('Not implemented: norm_total_pseudocount')
                    # Later if necessary, must adatpt
                    pseudocount = 0.001

                    p = (covs[p][d] + pseudocount) / (np.sum(covs[p][d]) + pseudocount)
                    q = (cov_mean + pseudocount) / (np.sum(cov_mean) + pseudocount)

                    js_dist = js_distance(p, q, base = 2)
                    js_div = js_dist ** 2

                    arr_metric[donor_idx, peak_idx] = js_div


        elif algorithm == 'wasserstein_dist':

            raise ValueError('Not implemented: wasserstein_dist')


        # Always compute counts

        for donor_idx, donor in enumerate(donors_sorted):

            arr_counts[donor_idx, peak_idx] = np.sum(covs[peak_id][donor])


            if algorithm == 'counts':

                arr_metric = arr_counts


                

    ## Create anndata

    if algorithm == 'js_divergence':

        # Scipy bug in `scipy.spatial.distance.jensenshannon`:
        #   returns NaNs when distributions are very closely similar but not exactly the same
        #   - Issue: https://github.com/scipy/scipy/issues/20083
        #   - Possible solution: https://github.com/scipy/scipy/pull/20786
        arr_metrix = np.nan_to_num(arr_metric, 0)


    adata = ad.AnnData(csr_matrix(arr_metric))
    adata.obs_names = donors_sorted


    peak_ids_classic = []

    for peak_idx, peak_id in enumerate(peak_ids):

        fields = peaks[peak_id]
        peak_ids_classic.append('{}:{}:{}:{}:*:{}'.format(fields[0],
                                                          fields[1],
                                                          fields[2],
                                                          fields[2] - fields[1],
                                                          peak_idx))

    adata.var_names = peak_ids_classic
    adata.var['centre_snp'] = peak_ids
    adata.var['peak_counts_total'] = np.sum(arr_counts, axis=0).astype(int)


    # Annotate

    adata.obs['donor'] = adata.obs.index

    if obs_annotation:
        adata.obs[obs_annotation[0]] = obs_annotation[1]


    # Save anndata

    utils.create_dir(out_adata)
    adata.write(out_adata, compression='gzip')


def main():

    ## Args

    args = parse_args()

    if args.obs_annotation:

        obs_annotation = args.obs_annotation.split(',')[0:2]

    else:

        obs_annotation = None


    ## Get and format peaks

    peaks = {}

    with open(args.peaks_tsv, 'r') as f:

        for line in f:

            fields = line.strip().split('\t')

            peak_name = fields[3]
            contig = fields[0]
            start = int(fields[1]) # Indexing basis: 0-based half-open
            end = int(fields[2])

            peaks[peak_name] = (contig, start, end)
    

    ## Compute footprints

    compute_footprints(coverage_files_dir = args.covs_dir,
                       out_adata = args.out_adata,
                       peaks = peaks,
                       algorithm = args.algorithm,
                       func_for_donor_naming = lambda x: os.path.splitext(x)[0],
                       obs_annotation = obs_annotation)


if __name__ == '__main__':

    main()
