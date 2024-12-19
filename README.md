## Environment setup


### Workflow

```mermaid
flowchart TB
subgraph fQTLs
E[make_matrix-eqtl_input_footprint-qtls.bash] --> D[process_footprint_adata.bash]
F[call_footprint-qtls_matrix-eqtl.bash] --> E
G[plot_qtl_results_footprints.bash] --> F
end
```


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
