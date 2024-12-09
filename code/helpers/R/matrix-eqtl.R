input_format_check <- function(gt.tsv, snp.loc.tsv, phe.tsv, peak.loc.tsv, cov.tsv) {
    # Check the tabular format of the input including sample name (columns) order.
    #
    # Todos:
    #    - Check for missing values and if they are the same as the user specified and expected

    phe.header <- colnames(fread(phe.tsv, nrows=1, header=TRUE))
    gt.header <- colnames(fread(gt.tsv, nrows = 1, header=TRUE))[-1]
    cov.header <- colnames(fread(cov.tsv, nrows=1, header=TRUE))[-1]

    # Ensure identical columns and column order (sample labels)

    
    if (all(gt.header == phe.header)) {
        
        print("CORRECT: Same columns and their order.")
        
    } else if (all(sort(gt.header) == sort(phe.header))) {
        
        stop("FORMAT ERROR: Genotype and phenotype files have same headers but different order!")
        
    } else {
        
	cat(gt.header)
	cat(phe.header)
        stop("FORMAT ERROR: Genotype and phenotype file headers don't match!")
        
    }
    

    if (all(gt.header == cov.tsv)) {
        
        print("CORRECT: Same columns and their order.")
        
    } else if (all(sort(gt.header) == sort(cov.tsv))) {
        
        stop("FORMAT ERROR: Genotype and covariate files have same headers but different order!")
        
    } else {
        
        stop("FORMAT ERROR: Genotype and covariate file headers don't match!")
        
    }
    

    if (all(phe.header == cov.tsv)) {
        
        print("CORRECT: Same columns and their order.")
        
    } else if (all(sort(phe.header) == sort(cov.tsv))) {
        
        stop("FORMAT ERROR: Phenotype and covariate files have same headers but different order!")
        
    } else {
        
        stop("FORMAT ERROR: Phenotype and covariate file headers don't match!")
        
    }
    

    # Ensure identical rows (snp & peak labels)
    
    gt.idxs <- fread(gt.tsv, select = 1, header = TRUE)[[1]]
    snp.loc.idxs <- fread(snp.loc.tsv, select = 1, header = TRUE)[[1]]
    
    if (!all(gt.idxs %in% snp.loc.idxs)) {
        
        stop("FORMAT ERROR: Genotype and snp location file indexes (snp labels) don't match!")
        
    }
    
    
    phe.idxs <- fread(phe.tsv, select = 1, header = TRUE)[[1]]
    peak.loc.idxs <- fread(peak.loc.tsv, select = 1, header = TRUE)[[1]]
    
    if (! all(phe.idxs %in% peak.loc.idxs)) {
        
        stop("FORMAT ERROR: Phenotype and peak file indexes (peak labels) don't match!")
        
    }

}



merge_covariates <- function(tsv1, tsv2, tsv.out) {
    # Merge covariate tsvs
    #
    # Note: LEGACY: not used I think.

    gt.pcs <- fread(tsv1, sep = '\t', header = TRUE) %>% 
        as.data.frame.matrix()
    
    rownames(gt.pcs) <- gt.pcs$id
    gt.pcs <- gt.pcs %>% select(-id)
    

    phe.pcs <- fread(PHE.PCS.FILE, sep = '\t', header = TRUE) %>% 
        as.data.frame.matrix()
    
    rownames(phe.pcs) <- phe.pcs$id
    phe.pcs <- phe.pcs %>% select(-id)
    
    
    if (all(colnames(gt.pcs) == colnames(phe.pcs))) {
        
        print("Colnames match.")
        
    } else if (all(sort(colnames(gt.pcs)) == sort(colnames(phe.pcs)))) {
        
        stop("Same headers/columns, different order!")
        
    } else {
        
        stop("Header don't match! Check order and equivalent sets.")
        
    }
    
    cvrts <- rbind(gt.pcs[phe.header], phe.pcs[phe.header]) %>%
        rownames_to_column(var = 'id')
    
    CVRT.FILE <- tsv.out

    return(CVRT.FILE)
    #fwrite(cvrts, CVRT.FILE, sep="\t")

}
