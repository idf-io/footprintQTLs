# General todos

- [ ] In format input anndatas, already make the cell_type column formatted. This requires to change operation when ct_format_alt is used downstream.

# Notes
- If ?: I haven't checked yet


# Scripts

glob_vars.bash:
    - :: Contains environmental variables that bash scripts may import
    
glob_vars.py
    - :: Contains environmental variables that python scripts may import

setup_data_cp_softlinks.bash
    - :: Makes links and copies from original data files and changes their name. Useful for data source tracing.
    - I:
        - PRECOMPUTED_EQTLS_TSV
    - O: ?
    - N: ?
    - T: ?
    - Todos: ?

format_genotype_data.bash(.py)
    - :: Correct snp naming and filter out indels of genotype matrix
    - R:
    - I:
        - GENOTYPES_TSV
    - O:
        - GENOTYPES_PROCESSED_TSV
    - L: ?
    - T: #genotype
    - Todos: If it becomes bigger, move the filtering to `process_genotype_data.py


make_extra_files_genotype.ipynb
    - :: Make some extra genotype-linked files
        - snp locations
    - R: jupyterhub: <1'
    - I:
        - GENOTYPES_TSV
    - O:
        - SNP_LOCS_BED
    - A:
        - Manu's genotype matrix: formatted indexes as 1_12323_A_T
    - L:
        - make_matrix-eqtl_input_*
        - calculations/*
    - T: #positions #locations #genotype



format_gene-expression_adata.ipynb ? For uniformity, not made yet
process_gene-expression_adata.ipynb ? For uniformity, not made yet


format_peaks_adata.ipynb
    - :: Formats the peaks anndata which will be used for the pipeline
    - R:
        - bsub: 75G, 10'
    - I:
        - ATAC_PEAKS_H5AD_OLD
        - RNA_H5AD_OLD
    - O: 
        - ATAC_PEAKS_H5AD_NEW
        - RNA_H5AD_NEW
    - A: RNA- and ATAC-seq anndata contents e.g. X are already computed and cell_types were annotated.
    - L: 
        - format_anndata.bash # BSUB launcher
        - process_peaks_adata.ipynb
    - T: #format #ct-map #json


process_peaks_adata.ipynb
    - :: Process the peaks anndata
        - to create cell-type pseudo-bulked anndatas (remove 'Discard' cell-types)
    - R:
        - bsub: 350G, 25'
    - I:
        - ATAC_PEAKS_H5AD_NEW
        - GENOTYPES_TSV
        - GENOTYPE_PCS_TSV
    - O:
        - ATAC_PEAKS_PROCESSED_H5AD
        - ATAC_CHROM_ACCESS_DIR/adata/
            - peak_matrix_cell-type-donors-pseudobulk.h5ad
            - CT_MAP_ID/<CT>/
                    - peak_matrix_cells_<CT>.h5ad
                            - peak_matrix_donors-pseudobulk_<CT>.h5ad
    - D:
        format_peaks_adata.ipynb
    - L: ?
        - process_peaks_adata.bash # BSUB launcher
    - T: #format #ct-map #json #donor-artifact #subset #discard


make_extra_files_peaks-accessibility.ipynb
    - :: Make some extra files from peaks accessibility adata (all cells)
        - Peak beds
    - R: jupyterhub: <1'
    - I:
        - FOOTPRINTS_DIR/ footprints_<CT>_processed.h5ad
    - O:
        - FOOTPRINTS_METADATA_DIR/ <CT>/ peaks.bed
    - A:
        - Peaks format: 1-based fully open
    - D:
        - compute_footprints.bash
    - L: ?
    - T: #footprints #locations #phenotype #bed


### caQTLs

make_matrix-eqtl_input_footprint-qtls_bulk-tests.bash
    - :: Make phenotype, genotype, peak location, snp location and covariates(chromatin acc. PCs and genotype PCs) tsvs
    - R: ?
    - I: ?
    - O: ?
    - A: ?
    - L: ?
    - D: ?
    - N:
        - Already exported different subsets and pseudobulks in process_peaks_adata.bash
        - Refactoring from "matrix-eQTL_input_refactored_annotation.ipynb"
    - T: ?
    - Todos: ?

call_ca-qtls_matrix-eqtl.R
call_ca-qtls_matrix-eqtl.bash


make_extra_files_peaks-accessibility_cell-type-level.ipynb
    - :: Make some extra files from peaks accessibility adata (pseudo-bulk: mean-aggregated across cells of same donor)
        - Peak beds
    - R: jupyterhub:
    - I:
        - FOOTPRINTS_DIR/ footprints_<CT>_processed.h5ad
    - O:
        - FOOTPRINTS_METADATA_DIR/ <CT>/ peaks.bed
    - A:
        - Peaks format: 1-based fully open
    - D:
        - compute_footprints.bash
    - L: ?
    - T: #footprints #locations #phenotype #bed


### Group fragment files

build_dir_tree.bash
    - :: Build the directories and softlinks necessary to create the template directory and file structure for project.
    - F: 
    - L:
    - D:
    - V: Newest

gather_group_cell-ids.py
    - :: Create group files containing the cells of the group (grouping based on annotation and grouping scheme)
    - F: Given an anndata.h5py with an annotation column and a group joining scheme in json format, create the files above.
    - L: split_fragments_by-groups.bash
    - D: code/helpers/fragment_file_helpers.bash

split_fragments_by-groups.bash
    - :: Takes fragment files from a folder and splits them according to a grouping config. Creates the cellid_group_files.
    - F: Given a fragment file and a folder with group files containing the cells pertaining to the group, split the file according to the groups.
    - L: join_fragments_by-groups.bash
    - D: 
        - code/helpers/fragment_file_helpers.bash
        - gather_group_cell-ids.py
    - EF:
        - check_jobs
        - (wait for cluster jobs)
    - N:
        - Fragment ids formatted as barcode_sample. h5ad must contain correctly labeled cells and group annotations. Fragment file names must be sample names.
        - !! CLuster and gathering part unfinished.
    - Todos: Join with script below and modularize wait for cluster
    - Tags: #jobs #cluster #cell-id #barcode

join_fragments_by-groups.bash
    - :: Takes fragment files in the tmp folder from `split_fragments_by-groups.bash` and joins them to donor_ct level .tsv.gz files
    - I:
        - SAMPLES_DIR="data/datasets/${DATASET}/atac-seq/fragment-files" --> donor_ct/sample%donor_ct.tsv.gz
        - OUT_DIR="data/intermediate-data/datasets/${DATASET}/atac-seq/fragment-files/grouped/${CT_MAP_KEY}/"
        + DATASET, CT_MAP_{JSON,KEY}
    - O:
        - OUTDIR/donor_ct.tsv.gz
    - L: compute_atac-coverage_bw.bash
    - D:
        - split_fragments_by-groups.bash
        - helpers/bash/utils.bash
        - helpers/bash/join_frag_files_by_groups.bash
    - N:
    - Tags: #jobs #cluster #gz #join


### footprints(QTLs)

compute_atac-coverage_bw.bash
    - :: Computes bw coverages from tsv.gz fragment files using ChromBPNet or scATAC-seq fragment tools.
    - F: frag_file_to_bw_chrombpnet
    - L: compute_footprints.bash
    - D: strip_and_order_fragment_file
    - N: The coverages are performed on fragment-file-level so these must represent the correct groups already (e.g. cell-type-donor)
    - T:
        #jobs #cluster #bw #chrombpnet


compute_footprints.bash
    - :: Computes the footprint/shape metric for all samples/donors with respect to the average profile and outputs an anndata object per cell type.
    - F: compute_footprints.py
        + DATASET, CT_MAP_JSON, CT_MAP_KEY
    - L: [...call fQTLs...]
    - D:
    - N: 
        - The footprints are computed using bw(bg) files on a sample/donor-level, a folder with only these bw files is required.
        - Discar: `Discard` group
    - T: #jobs #cluster #footprint #jensen-shannon #shape


pre-annotate_footprint_adata.py
    - :: Annotate the anndata object with
         donor, donor_id, ncells, nfrags, nins per donor and obs and var metadata from reference anndata
    - R:
        - t: 8' compute
        - m: 45G
    - I:
        - FOOTPRINTS_DIR/ footprints_<CT>_raw.h5ad
        - GROUPED_FRAG_FILES_DIR/ <donor>_<CT>.tsv.gz
        - ATAC_PEAKS_H5AD_NEW
    - O:
        - FOOTPRINTS_DIR/ footprints_<CT>_pre-annotated.h5ad
    - A:
        - Group fragment files generated above contain 
          all cells and fragments used for coverage and footprint computation.
        - Reference anndata correctly depicts information about the current. 
          AKA vars of same id can be assumed to be the same, obs atac-seq data (pseudobulked) 
          was used to create the current anndata object. 
    - L: process_footprint_adata.ipynb
    - D: compute_footprints.bash
    - T: #cluster #annotation #simple-launcher #latest-launcher #simple-bsub
    - Todo:
        - Check if first A is true.


process_footprint_adata.bash
    - :: Process the pre-annotated footprint anndatas:
        - min_cells_donor, subset donors to common in all data sources, descr. stats, PCA, clustering, UMAP
    - R:
        - bsub: 4G, 20'
    - I:
        - FOOTPRINTS_DIR/ footprints_<CT>_pre-annotated.h5ad
        - GENOTYPES_TSV
        - GENOTYPE_PCS_TSV
    - O:
        - FOOTPRINTS_DIR/ footprints_<CT>_processed.h5ad
    - L:
        - footprints_eda_notebook.ipynb
        - footprints_eda_all-cts_notebook.ipynb
    - D: pre-annotate_footprint_adata.bash
    - N: Currently executes but doesn't end. Probably some danling open file.
    - T: #job-array #filter-ncells #threshold #donors-artifact #subset #donors-artifact-redundant?


footprints_eda.bash
    - :: Exploratory data analysis on anndata of footprints (donors x peaks) at the cell-type level.
    - R:
        - bsub: 1G, 14'
    - I:
        - FOOTPRINTS_DIR/ footprints_<CT>_processed.h5ad
        - GROUPED_BIGWIG_FILES_DIR/ <CT>/<donor>_<CT>.bw
    - O:
        - Figs: descr. stats, PCA, UMAP, profile comparison
        - FOOTPRINTS_EDA
                          /eda_<CT>.ipynb
                          /<CT>/
                                  footprints_<CT>.pdf
    - D: process_footprint_adata.ipynb
    - T: #eda #profiles #pdf #ipynb


footprints_eda_all-cts.bash
    - :: Exploratory data analysis on anndata of footprints (donors x peaks) across all cell-types.
    - R:
        - bsub: 1G, 1'
        - notebook: 1'
    - I:
        - FOOTPRINTS_DIR/ footprints_<CT>_processed.h5ad
    - O:
        - Figs: descr. stats, HVP set analysis
        - FOOTPRINTS_EDA
                          /eda_all-cell-types.ipynb
    - D: process_footprint_adata.ipynb
    - T: #eda #ipynb #upsetplot
    - Todos:
        - Variance analysis: cell-types vs donors


make_extra_files_footprints_cell-type-level.ipynb
    - :: Make some extra files from footprints adata
        - Peak beds
    - R: jupyterhub: <1'
    - I:
        - FOOTPRINTS_DIR/ footprints_<CT>_processed.h5ad
    - O:
        - FOOTPRINTS_METADATA_DIR/ <CT>/ peaks.bed
    - A:
        - Peaks format: 1-based fully open
    - D:
        - compute_footprints.bash
    - L:
        - calculations/*
    - T: #footprints #locations #phenotype #bed


make_matrix-eqtl_input_footprint-qtls.bash (-->ipynb)
    - :: Make phenotype, genotype, peak location, snp location and covariates(footprint PCs and genotype PCs) tsvs
    - R:
        - bsub:
            - mode=single-tests: 6G, 100-190'
            - mode=bulk-tests: 215G, 10'
    - I:
        - CT_MAP_JSON
        - FOOTPRINTS_DIR/ footprints_<CT>_processed.h5ad
        - GENOTYPES_TSV
        - GENOTYPES_VCF
        - GENOTYPE_PCS_TSV
    - O:
        - MATRIX_EQTL_INPUT_DIR/ <CT>
                                /bulk-tests/
                                    make_matrix-eqtl_input_footprint-qtls.ipynb
                                    footprints.tsv
                                    peak_locations.tsv
                                    footprint_pcs.tsv
                                    genotype_NA.tsv
                                    snp_locations.tsv
                                    genotype_pcs.tsv
                                single-tests/
                                    make_matrix-eqtl_input_footprint-qtls.ipynb
                                    phenotypes/
                                    peak_locations/
                                    footprint_pcs/
                                    genotypes/
                                    snp_locations/
                                    covariates/
        - MATRIX_EQTL_OUTPUT_DIR/qtl_input_stats.tsv
    - A: donors/obs in anndata are also found in genotype databases (processing step in process_footprint_adata.ipyb)
    - D: process_footprint_adata.bash
    - L: call_footprint-qtls_matrix-eqtl.bash
    - T: #qtls #io #memory-control #garbage-collector
    - Todos:
        - [ ] Make individual steps functions if reused in other scenarios
    - Note:
        - Major parameter: mode={bulk-tests, single-tests, peak-tests}

call_footprint-qtls_matrix-eqtl.bash
    - :: Call QTLs on the footprint data
    - R: bsub(multiple=500): 0.5M, 25' bsub(multiple=5000): 0.5M, 10'
    - I:
        - MATRIX_EQTL_INPUT_DIR
    - O:
        - MATRIX_EQTL_OUTPUT_DIR
    - A:
    - L: 
    - D: helpers/R/matrix-eqtl.R
    - N:
    - T: #qtl #job-array #alpha #cis-dist
    - Todos:
    - Note: renv has problems with many parallel cluster jobs

plot_qtl_results_footprints.bash --> .R
    - :: Collect, FDR and plot results of QTL testing
    - R: bsub: 0.5M, 6'
    - I: MATRIX_EQTL_OUTPUT_DIR
    - O: MATRIX_EQTL_OUTPUT_DIR/
                                qtls_all.tsv
                                qtls_all_fdr.tsv
                                qtls_stats.tsv
    - L:
    - D: call_footprint-qtls_matrix-eqtl.bash 
    - N:
    - T: #fdr
    - Todos:




process_single-test_qtl-calling_footprint-qtls.bash
    - :: Merge the tests and compute statistics and plots from the output directory of qtl-calling of single-tests.
    - R:
    - I:
    - O:
    - A:
    - L:
    - D:
    - N:
    - T:
    - Todos:





template.ipynb
    - :: Tldr
    - R:
    - I:
    - O:
    - A:
    - L:
    - D:
    - N:
    - T:
    - Todos:





# QTL mapping pipeline

matrix-eQTL_input_refactored_annotation.ipynb




matrix-eqtl_args.R
matrix-eQTL_input_refactored.ipynb
matrix-eQTL_input_refactored2.bash
matrix-eqtl_args.cp.R
caQTLS_analysis.ipynb
caqtl_stats.R
matrix-eQTL_input.ipynb

matrix-eqtl_restricted.R
matrix-eQTL_input_refactored.R


matrix-eQTL_input_refactored2/
