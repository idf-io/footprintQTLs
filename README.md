## Workflow

```mermaid
flowchart TD

subgraph Setup
S1[format_genotype_data.bash]
S2[make_extra_files_genotype.ipynb]
S3[format_peaks_adata.ipynb]
S4[process_peaks_adata.bash]
S3 --> S4
S6[format_fragment-files.bash]
end

X1[gather_group_cell-ids.py]
X2[split_fragments_by-groups.bash]
X3[join_fragments_by-groups.bash]

subgraph Group fragment file
S6 --> X1
X1 --> X2
X2 --> X3
end

K1[select_peaks_and_create_ca-qtl_adata.bash]
K2[make_matrix-eqtl_input_chromatin-accessibility-qtls.bash]
K3[call_chromatin-accessibility-qtls_matrix-eqtl.bash]
K4[plot_qtl_results_chromatin-accessibility.bash]

subgraph caQTLs
S1 & S2 & S4 --> K1
K1 --> K2
K2 --> K3
K3 --> K4
end


A[compute_atac-coverage_bw.bash]
B[compute_footprints.bash]
C[pre-annotate_footprint_adata.bash]
D[process_footprint_adata.bash]
E[make extra_files_footprints_cell-type-level]
D1[footprints_eda.bash]
D2[footprints_eda_all-cts.bash]
F[make_matrix-eqtl_input_footprint-qtls.bash]
G[call_footprint-qtls_matrix-eqtl.bash]
H[plot_qtl_results_footprints.bash]

K2 --> A2[generate_phenotype-modality_and_define-tests.ipynb]
A2 --> B & F

subgraph fQTLs
X3 --> A
A --> B
B --> C
C --> D
D ---> E
D --> D1 & D2
D & E --> F
S1 & S2 --> F
F --> G
G --> H
end

Z[gather-and-plot_qtl_metatable.ipynb]
H & K4 --> Z
```

**Legend**
- Squared blocks: scripts
- Rounded corner square blocks: function files
- Hexagonal blocks: non-existent files
- Left indented squared blocks: unfinished file
- Dotted line: path not used yet
- X-ended arrows: Not intended downstream (e.g. temporary file)
- o-ended arrows: Newer/updated version of the file


## Environment setup

### R

For reproducibility, `renv.lock` and `DESCRIPTION` files were provided for project-level R-package dependencies setup using `renv`.

Execute the following script to setup the `renv` environment:

```
renv::init(bioconductor = TRUE)
# renv will ask which dependencies to use -> `explicit` from DESCRIPTION file

renv::restore()
```

[renv's](https://rstudio.github.io/renv/index.html) documentation and articles are very useful for setup and package installation debugging.

For setup on the DKFZ odcf cluster (dated: 2024-11-27) follow additional configuration steps in `.setup-odcf-env.bash` BEFORE installing packages with `renv` as above. This scripot also creates the file `load-odcf-env.bash` to load specific modules. Source when working with the project.
