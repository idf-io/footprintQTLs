import gzip
from typing import Tuple, Dict


def get_counts_from_fragment_file(frag_file: str, types: Tuple[str] = ('fragments')):
    """
    Parse fragment file (bed format) and extract any combination of the following count types:
        - fragments
        - insertions (fragments * 2)
        - ids

    Notes:
        - the ids counts can be useful if they represent a cell-id (must be unique for the cell)

    Assumptions:
        - for ids, the 4th column ID must be unique for every category it represents
    """

    types = tuple(types)

    # Checks

    frag_file_extensions = ('.tsv', '.bed', '.tsv.gz', '.bed.gz')
    assert any(frag_file.endswith(ext) for ext in frag_file_extensions)

    types_options = ('fragments', 'insertions', 'ids')
    assert all([c in types_options for c in types])


    # Read file

    if frag_file.endswith('.gz'):

        with gzip.open(frag_file, 'rt') as f:

            lines = f.readlines()

    else:

        with open(frag_file, 'r') as f:
        
            lines = f.readlines()

    
    # Parse

    frag_count = 0
    cell_ids = []


    for l in lines:

        l = l.strip()

        if l.startswith('#'):

            continue


        frag_count += 1
        cell_ids.append(l.split('\t')[3])
        
    insertion_count = frag_count * 2
    id_count = len(set(cell_ids))


    types_count_map = {'fragments': frag_count, 'insertions': insertion_count, 'ids': id_count}

    return (types_count_map[t] for t in types)




def annotate_adata_counts_from_fragment_files(adata: str, donor_frag_file_map: Dict[str, str], types: Tuple[str] = ('fragments')):
    """
    Annotate anndata with any combination of counts: {fragments, insertions or ids} per donor.

    Input:
        - adata (anndata)
        - donor_frag_file_map (dict): donors -> fragment file path map to annotate the anndata by.

    Output:
        - adata with columns adata.obs['n_fragments', 'n_insertions', 'n_ids']

    If completely broken see legacy function below.
    """


    count_maps = {'fragments': {}, 'insertions': {}, 'ids': {}} # {'type': {'donor': count}}

    for donor, frag_file in donor_frag_file_map.items():


        n_frags, n_ins = get_counts_from_fragment_file(frag_file, types=['fragments', 'insertions'])

        count_maps['fragments'][donor] = n_frags
        count_maps['insertions'][donor] = n_ins
        

    # Annotate anndata using maps
    adata_annotated = adata.copy()

    for t in types:

        adata_annotated.obs['n_' + t] = adata_annotated.obs.index.map(count_maps[t]).fillna(0).astype(int)


    return adata_annotated



def _annotate_nfrags_nins_adata(adata: str, obs_fragment_files_map: dict):
    """
    ### LEGACY ###
    -> Improved by annotate_nfrags_nins_adata


    Annotate anndata with nr. of fragments and insertions per observation from obs(e.g. donor)-mapped fragment files.

    Input:
        - adata (anndata)
        - obs_fragment_files_map (dict): observation -> fragment file map to annotate the anndata by.

    Output:
        - adata with columns adata.obs['n_frags', 'n_insertions']

    Assumptions:
        - When using a fragment file: the file was the sole source of fragments for creating the footprints anndata.
    """

    # Checks
    for f in obs_fragment_files_map.values():
        assert f.endswith('.tsv') or f.endswith('.tsv.gz')


    n_frags_map = {} # obs: n_frags
    n_insertions_map = {} # obs: n_insertions

    for obs, frag_file in obs_fragment_files_map.items():

        # Construct bash command

        if frag_file.endswith('.gz'):

            command_count_lines = f'gunzip -c \'{frag_file}\' | grep -v \'^#\' | wc -l'

        else:

            command_count_lines = f'grep -v \'^#\' \'{frag_file}\' | wc -l'


        # Count lines
        result = subprocess.run(command_count_lines, shell=True, capture_output=True, text=True, check=True)

        if result.stderr:
            raise RuntimeError(f'Command returned an error: {result.stderr}')

        n_lines = int(result.stdout)

        n_frags_map[obs] = n_lines
        n_insertions_map[obs] = n_lines * 2

    
    # Annotate anndata using maps
    adata_annotated = adata.copy()
    adata_annotated.obs['n_frags'] = adata_annotated.obs.index.map(n_frags_map).fillna(0).astype(int)
    adata_annotated.obs['n_insertions'] = adata_annotated.obs.index.map(n_insertions_map).fillna(0).astype(int)

    return adata_annotated
