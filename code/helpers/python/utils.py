import os
import io
import re
from typing import List, Dict, Callable, Tuple
import numpy as np
import pandas as pd
import anndata as ad



def create_dir(path: str, no_extension: bool = False, path_type: str = None):
    """
    If input is a directory, creates the directory, 
    if the input is a file name, then the parent directory of the file is created.
    Optionally indicate to get the parent of the string.

    Input:
	- path (str): Directory or file path.
	- no_extension (bool): Use True if the file has no extension
	- path_type (str) ['file' | 'dir'] optional: Explicitly indicate the query type.

    Assumptions:
	- Directory name doesn't include '.'
    """

    path = path.strip()
    path = os.path.abspath(path)


    def make_dir_from_dir(path: str):

    	os.makedirs(path, exist_ok=True)


    def make_dir_from_file(path: str):

    	os.makedirs(os.path.dirname(path), exist_ok=True)


    if path_type == 'file' or (path_type == None and (no_extension == True or '.' in os.path.basename(path))):

    	make_dir_from_file(path)


    elif path_type == 'dir' or path_type is None:
    
    	make_dir_from_dir(path)



def list_files_and_links(
        directory: str = '.',
        extension: str = None):
    """
    List only files and symbolic links in the given directory
    
    Args:
        directory (str): Path to directory to scan. Defaults to current directory.
        extension (str): Optional filter for file extension.
    
    Returns:
        list: Names of files and symbolic links in the directory
    """
    
    search_results = []

    if extension:

        for entry in os.listdir(directory):

            file_path = os.path.join(directory, entry)

            if entry.endswith(extension) and (os.path.isfile(file_path) or os.path.islink(file_path)):

                search_results.append(entry)

    else:

        for entry in os.listdir(directory):

            file_path = os.path.join(directory, entry)

            if os.path.isfile(file_path) or os.path.islink(file_path):

                search_results.append(entry)

    return search_results


def list_dirs(directory='.'):
    """
    List directories in query directory. Depth=1
    """

    return [e for e in os.listdir(directory) if os.path.isdir(os.path.join(directory, e))]


def read_rownames(file: str,
                  header: bool = False,
                  sep: str = '\t'):
    """
    Extract the values of the first column of a line-based file.
    """

    row_ids = []

    with open(file, 'r') as lines_file:

        if header:

            next(lines_file)

        for line in lines_file:

            row_id = line.strip().split(sep)[0]
            row_ids.append(row_id)

    return row_ids


def middle_indexes(length: int, radial_range: int = 1):
    """
    Get indexes of the middle. Optionally indicate how many elements should result.
    """
    c = length // 2
    s = radial_range // 2

    if radial_range % 2 != 0:
        c_u = c + s
    else:
        c_u = c + s - 1

    c_d = c - s

    return [c_d, c_u+1]



### SORT KEYS ###

def sort_key_batches(batch):
    """
    Sort donor names of format: sSL0123A
    """

    match = re.search(r'sSL([0-9]+)([a-zaA-Z]*)', batch)
    
    return (match.group(1), match.group(2))


def sort_key_peaks(peak):
    """
    Sort peak names of format chrX:start:end and contigX:start:end
    """

    contig, start, end  = peak.split(':')[0:3]
    
    if contig.startswith('chr'):
        # If contig format: `chr13`
        contig_nr = int(contig[3:]) if contig[3:].isdigit() else float('inf')

        primary_sort = (0, contig_nr)

    else:

        primary_sort = (1, contig)

    return (*primary_sort, int(start), int(end))


### FORMAT ###

def ct_format(cell_type):
    return cell_type.replace(' ', '-').replace('.', '')


def ct_format_alt(cell_type):
    return cell_type.replace(' ', '_').replace('.', '')


def regex_get(text, pattern=None):
    
    pattern_map = {'sample_id': r'sSL\d*[A-Za-z]?',
                   'barcode': r'^\D*-\d*'}

    if pattern in pattern_map.keys():
        pattern = pattern_map[pattern]
    
    match = re.search(pattern, text)
    return match.group()




### CHECKS ###

def check_peak_format(peak):
    """
    Check if the format of the peak name is known. Optionally output the indexing basis.
    """

    pattern = r'^[a-zA-Z0-9]{3,}:(\d+):(\d+):(\d+)'
    match = re.search(pattern, peak)

    if match:

        diff = int(match.group(2)) - int(match.group(1))

        if diff == int(match.group(3)):

            return "0-based half-open"

        elif diff == int(match.group(3)) + 1:

            return "1-based fully-closed"

    else:
        raise ValueError(f'Peak format not recognized. Expected: {pattern}, not {peak}')



### Parsers ###

def parse_vcf(vcf_file: str):
    """
    Parse a VCF file from either a file path or content string.
    
    Input:
        vcf_file (str): VCF file path.
        
    Returns:
        pandas.DataFrame: Parsed VCF data

    Assumptions:
        VCF file first line starts with '#CHROM' and comment lines are preceded with '#'
    """

    with open(vcf_file, 'r') as f:
        content = f.read()
    
    lines = content.strip().split('\n')
    
    # Find the header line (starts with #CHROM)
    header_line = None

    for i, line in enumerate(lines):

        if line.startswith('#CHROM'):

            header_line = line
            data_start_index = i + 1

            break
    
    if not header_line:
        raise ValueError("Could not find header line starting with #CHROM")


    # Create a string buffer with the header and data
    header_line = header_line.lstrip('#')
    data_str = '\n'.join([header_line] + lines[data_start_index:])
    
    # Read the tsv data using pandas
    df = pd.read_csv(io.StringIO(data_str), 
                     sep='\t',
                     comment='#',
                     dtype={'POS': 'Int64', 'QUAL': 'float64'})
    
    return df


### Bed tools port ###

def bedtools_intersect(
    bed_a: str, 
    bed_b: str, 
    flags: Tuple[str] = None, 
    output_file: str = None
) -> str:
    """
    Python port for `bedtools intersect`.
    
    Args:
        bed_a (str): Path to the first BED file
        bed_b (str): Path to the second BED file
        flags (optional, Tuple[str]): Flags for bedtools intersect
        output_file (optional, str): Output file path
    
    Returns:
        str: Command output or contents of output file
    """

    # Checks
    assert all(os.path.exists(file) for file in [bed_a, bed_b])
    
    if flags:
        assert all(f.startswith('-') for f in flags), "All flags must start with '-'"
    
    # Construct command
    cmd = ['bedtools', 'intersect', '-a', bed_a, '-b', bed_b]
    if flags:
        cmd.extend(flags)
    
    print(f'Performing: shell$ {" ".join(cmd)}')
    
    try:

        # If output file is specified, write directly to it
        if output_file:

            with open(output_file, 'w') as f:

                subprocess.run(cmd, check=True, stdout=f, text=True)
            
            # Read and return file contents
            with open(output_file, 'r') as f:
                
                return f.read()
        
        # If no output file, capture and return output
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return result.stdout
    
    except subprocess.CalledProcessError as e:
        print(f'Command failed with return code: {e.returncode}')
        print(f'Error output: {e.stderr}')
        raise
