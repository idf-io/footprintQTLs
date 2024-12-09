# Error handling function
handle_error <- function(e) {
    cat("Error: ", conditionMessage(e), "\n")
    cat("Please check your inputs and try again.\n")
    quit(status = 1)
}

tryCatch({
    
    PROJ.ROOT <- "/home/fichtner/projects/footprintQTL"
    #cat("Number of arguments:", length(args), "\n")
    #print(args)
    
    ## Libraries
    
    if (!requireNamespace("renv", quietly = TRUE)) {
        install.packages("renv")
    }
    
    setwd(paste0(PROJ.ROOT, '/code'))
    library(renv)
    renv::activate()
    setwd(PROJ.ROOT)
    
    library(MatrixEQTL)
    library(data.table)

    library(dplyr)

    library(tibble)
    library(stringr)

    source("code/helpers/R/utils.R")
    
    
    ##  User parameters
    
    cell.type.r <- "Neural progenitors"
    job.id <- "today"
    ct.map.id <- "approach_2024-09-12"
    #ct.map.id <- "original_clean"
    dataset.id <- "hca_brain-organoids" # ! Change DATA.DIR when 'hca_brain-organoids_`processed`' is removed
    
    MODEL <-  modelLINEAR
    ALPHA <-  1e-2
    ERROR.COV <-  numeric()
    CIS.DIST <-  0 # As Aygun?
    #CIS.DIST <-  1000 # As Aygun? 1k for normal caQTLS!
    
    ## Args
    
    #args <- commandArgs(trailingOnly = TRUE)
    
    #if (length(args) == 4) {
        
        #cell.type.r <- as.character(args[1])
        #job.id <- as.character(args[2])
        #ct.map.id <- as.character(args[3])
        #dataset <- as.character(args[4])
        
    #} else {
        
        #print("Wrong number of arguments!")
        #stop()
        
    #}
    cat("Args:")
    cat(cell.type.r)
    cat(job.id)
    cat(ct.map.id )
    cat(dataset.id)
    
    cell.type = ct.format(cell.type.r)
    cell.type.alt = ct.format.alt(cell.type.r)
    
    
    ## Dir variables
    
    DATA.DIR <-  paste0("data/datasets/", dataset.id, "/")
    DATA.INTERM.DIR <- paste0("data/intermediate-data/datasets/", dataset.id, "/matrix-eQTL/footprints/", ct.map.id)
    RESULTS.DIR <-  paste0("results/datasets/", dataset.id, "/matrix-eQTL/footprints", ct.map.id, "/")
    
    STATS.FILE <- paste0(RESULTS.DIR, "caQTL_stats.tsv")
    OUT.FILE.SIMPLE <- paste0(RESULTS.DIR, cell.type, "/caqtls_peaks_simple_", cell.type, ".tsv")
    OUT.FILE <- paste0(RESULTS.DIR, cell.type, "/caqtls_peaks_", cell.type, ".tsv")
    
    SNP.GT.FILE <-  paste0(DATA.INTERM.DIR, "/", cell.type, "/genotype_NA.tsv")
    SNP.LOC.FILE <-  paste0(DATA.INTERM.DIR, "/", cell.type, "/snp_locations.tsv")
    PEAKS.CA.FILE <- paste0(DATA.INTERM.DIR, "/", cell.type, "/footprints.tsv")
    PEAKS.LOC.FILE <- paste0(DATA.INTERM.DIR, "/", cell.type, "/peak_locations.tsv")
    
    GT.PCS.FILE <- paste0(DATA.INTERM.DIR, "/", cell.type, "/genotype_pcs.tsv")
    CA.PCS.FILE <- paste0(DATA.INTERM.DIR, "/", cell.type, "/footprint_pcs.tsv")


    #DATA.DIR <-  "data/datasets/hca_brain-organoids_processed/"
    #DATA.INTERM.DIR <- paste0("data/intermediate-data/datasets/", dataset.id, "/matrix-eQTL_io/")
    #RESULTS.DIR <-  paste0("results/datasets/", dataset.id, "/matrix-eQTL/", ct.map.id, "/")
    
    #STATS.FILE <- paste0(RESULTS.DIR, "caQTL_stats.tsv")
    #OUT.FILE.SIMPLE <- paste0(RESULTS.DIR, cell.type, "/caqtls_peaks_simple_", cell.type, ".tsv")
    #OUT.FILE <- paste0(RESULTS.DIR, cell.type, "/caqtls_peaks_", cell.type, ".tsv")
    
    #SNP.GT.FILE <-  paste0(DATA.INTERM.DIR, "genotype/", ct.map.id, "/genotype_NA_", cell.type , ".tsv")
    #SNP.LOC.FILE <-  paste0(DATA.INTERM.DIR, "genotype/snp_locations/", ct.map.id, "/", cell.type , ".tsv")
    #PEAKS.CA.FILE <- paste0(DATA.INTERM.DIR, "chromatin_accessibility/peak_ca/", ct.map.id, "/", cell.type , ".tsv")
    #PEAKS.LOC.FILE <- paste0(DATA.INTERM.DIR, "chromatin_accessibility/peak_locations/", ct.map.id, "/", cell.type , ".tsv")
    
    #GT.PCS.FILE <- paste0(DATA.INTERM.DIR, "genotype/genotype_pcs/", ct.map.id, "/", cell.type , ".tsv")
    #CA.PCS.FILE <- paste0(DATA.INTERM.DIR, "chromatin_accessibility/ca_pcs/", ct.map.id, "/", cell.type, ".tsv")

    
    files <- c(DATA.DIR,
              DATA.INTERM.DIR,
              RESULTS.DIR,
              #
              STATS.FILE,
              OUT.FILE.SIMPLE,
              OUT.FILE,
              #
              SNP.GT.FILE,
              SNP.LOC.FILE,
              PEAKS.CA.FILE,
              PEAKS.LOC.FILE,
              #
              GT.PCS.FILE,
              CA.PCS.FILE)
    
    for (d in files) {
        
        create_parent_dir(d)
        
    }
    
    
    ## Get colnames for sanity check
    snp.gt.header <- colnames(fread(SNP.GT.FILE, nrows = 1, header=TRUE))[-1]
    peaks.ca.header <- colnames(fread(PEAKS.CA.FILE, nrows=1, header=TRUE))[-1]
    
    
    ## Make covariates file
    
    gt.pcs <- fread(GT.PCS.FILE, sep = '\t', header = TRUE) %>% 
        as.data.frame.matrix()
    
    rownames(gt.pcs) <- gt.pcs$id
    
    gt.pcs <- gt.pcs %>% select(-id)
    
    ca.pcs <- fread(CA.PCS.FILE, sep = '\t', header = TRUE) %>% 
        as.data.frame.matrix()
    
    rownames(ca.pcs) <- ca.pcs$id
    
    ca.pcs <- ca.pcs %>% select(-id)
    
    
    if (all(colnames(gt.pcs) == colnames(ca.pcs))) {
        
        print("Good.")
        
    } else if (all(sort(colnames(gt.pcs)) == sort(colnames(ca.pcs)))) {
        
        stop("Same headers, different order!")
        
    } else {
        
        stop("Headers don't match! Check order and equivalent sets.")
        
    }
    
    cvrts <- rbind(gt.pcs[snp.gt.header], ca.pcs[snp.gt.header]) %>%
        rownames_to_column(var = 'id')
    
    fwrite(cvrts, paste0(DATA.INTERM.DIR, "/", cell.type, "/covariates.tsv"), sep="\t")
    
    CVRT.FILE <- paste0(DATA.INTERM.DIR, "/", cell.type, "/covariates.tsv")
    
    
    ## Sanity check: Ensure correct rows(indexes) and columns
    
    cvrt.header <- colnames(fread(CVRT.FILE, nrows=1, header=TRUE))[-1]
    
    if (all(snp.gt.header == peaks.ca.header)) {
        
        print("Good.")
        
    } else if (all(sort(snp.gt.header) == sort(peaks.ca.header))) {
        
        stop("Same headers, different order!")
        
    } else {
        
        print("Headers don't match! Check order and equivalent sets.")
        
    }
    
    if (all(snp.gt.header == cvrt.header)) {
        
        print("Good.")
        
    } else if (all(sort(snp.gt.header) == sort(cvrt.header))) {
        
        stop("Same headers, different order!")
        
    } else {
        
        stop("Headers don't match! Check order and equivalent sets.")
        
    }
    
    if (all(peaks.ca.header == cvrt.header)) {
        
        print("Good.")
        
    } else if (all(sort(peaks.ca.header) == sort(cvrt.header))) {
        
        stop("Same headers, different order!")
        
    } else {
        
        stop("Headers don't match! Check order and equivalent sets.")
        
    }
    
    
    snps.gt.idxs <- fread(SNP.GT.FILE, select = 1, header = TRUE)[[1]]
    snps.loc.idxs <- fread(SNP.LOC.FILE, select = 1, header = TRUE)[[1]]
    
    if (!all(snps.gt.idxs %in% snps.loc.idxs)) {
        
        stop("Indexes don't match! Check order and equivalent sets.")
        
    }
    
    
    peaks.ca.idxs <- fread(PEAKS.CA.FILE, select = 1, header = TRUE)[[1]]
    peaks.loc.idxs <- fread(PEAKS.LOC.FILE, select = 1, header = TRUE)[[1]]
    
    if (! all(peaks.ca.idxs %in% peaks.loc.idxs)) {
        
        stop("Indexes don't match! Check order and equivalent sets.")
        
    }
    
    
    
    ## Calculate QTLs
    
    snps.gt = SlicedData$new()
    snps.gt$fileDelimiter = "\t"      # the TAB character
    snps.gt$fileOmitCharacters = "NaN" # denote missing values
    snps.gt$fileSkipRows = 1          # one row of column labels
    snps.gt$fileSkipColumns = 1       # one column of row labels
    snps.gt$fileSliceSize = 2000      # read file in pieces of 2,000 rows
    snps.gt$LoadFile( SNP.GT.FILE )
    
    snps.pos <- read.table(SNP.LOC.FILE, sep = '\t', header = TRUE)
    
    peaks.ca = SlicedData$new()
    peaks.ca$fileDelimiter = "\t"      # the TAB character
    peaks.ca$fileOmitCharacters = "NaN" # denote missing values
    peaks.ca$fileSkipRows = 1          # one row of column labels
    peaks.ca$fileSkipColumns = 1       # one column of row labels
    peaks.ca$fileSliceSize = 2000      # read file in pieces of 2,000 rows
    peaks.ca$LoadFile( PEAKS.CA.FILE )
    
    peaks.pos <- read.table(PEAKS.LOC.FILE, sep = '\t', header = TRUE)
    
    cvrt = SlicedData$new()
    cvrt$fileDelimiter = "\t"      # the TAB character
    cvrt$fileOmitCharacters = "NaN" # denote missing values
    cvrt$fileSkipRows = 1          # one row of column labels
    cvrt$fileSkipColumns = 1       # one column of row labels
    cvrt$fileSliceSize = 2000      # read file in pieces of 2,000 rows
    cvrt$LoadFile( CVRT.FILE )
    
    
    
    me = Matrix_eQTL_main(
        snps = snps.gt,
        snpspos = snps.pos,
        gene = peaks.ca,
        genepos = peaks.pos,
        cvrt = cvrt,
        output_file_name.cis = OUT.FILE.SIMPLE,
        pvOutputThreshold = 0, # Skip trans-eQTL mappinh
        pvOutputThreshold.cis = ALPHA,
        cisDist = CIS.DIST,
        useModel = MODEL,
        errorCovariance = ERROR.COV,
        verbose = TRUE,
        pvalue.hist = TRUE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE)
    
    unlink(OUT.FILE.SIMPLE)
    
    cat("1")
    ## Results:
    
    cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
    cat('Detected QTLs:', '\n');
    show(me$cis$eqtls)
    
    # Times
    # Midbrain EN - DRD4 - 11 016 cis-tests = 33s
    # Added CA PCs as covariate: cis-tests = 44s
    
    ## Plot the histogram of all p-values
    png(paste0(RESULTS.DIR, cell.type, "_hist_", cell.type, ".png"))
    
    plot(me)
    
    dev.off()

    ## Manhattan plot
    library(qqman)
    
    caqtls <- me$cis$eqtls
    caqtls <-  caqtls %>%
        mutate(snp_chr = sapply(str_split(snps,'_'), function(z) as.numeric(str_sub(z[1], 4)))) %>% 
        mutate(snp_pos = sapply(str_split(snps,'_'), function(z) as.numeric(z[2]))) %>% 
        mutate(snp_ref = sapply(str_split(snps,'_'), function(z) z[3])) %>% 
        mutate(snp_all = sapply(str_split(snps,'_'), function(z) z[4])) %>% 
        mutate(peak_chr = sapply(str_split(gene,':'), function(z) as.numeric(str_sub(z[1], 4)))) %>%
        mutate(peak_start = sapply(str_split(gene,':'), function(z) as.numeric(z[2]) - 1)) %>% 
        mutate(peak_end = sapply(str_split(gene,':'), function(z) as.numeric(z[3]))) %>% 
        mutate(peak_len = sapply(str_split(gene,':'), function(z) as.numeric(z[4]))) %>%
        arrange(desc(FDR)) 
        
        
    
    
    # snps per chromosome
    # as.data.frame(table(gwasResults$CHR))
    
    ## Plots
    
    # Manhattan FDR
    png(paste0(RESULTS.DIR, cell.type, "manhattan-plot_fdr_", cell.type, ".png"))
    manhattan(caqtls, chr = 'snp_chr', bp = 'snp_pos', p = 'FDR', snp = 'snps')
    dev.off()
    
    # Manhattan nominal
    png(paste0(RESULTS.DIR, cell.type, "manhattan-plot_nominal_", cell.type, ".png"))
    manhattan(caqtls, chr = 'snp_chr', bp = 'snp_pos', p = 'pvalue', snp = 'snps')
    dev.off()
    
    # Q-Q nominal
    png(paste0(RESULTS.DIR, cell.type, "qq_plot_nominal_", cell.type, ".png"))
    qq(caqtls$pvalue, main = paste0("Q-Q plot:", cell.type, "pvalue"))
    dev.off()
    
    # Q-Q FDR
    png(paste0(RESULTS.DIR, cell.type, "qq_plot_fdr_", cell.type, ".png"))
    qq(caqtls$FDR, main = paste0("Q-Q plot:", cell.type, "fdr"))
    dev.off()
    
    ## Save results
    
    # caQTLs
    write.table(caqtls, file = OUT.FILE, sep = "\t", quote = FALSE)
    
    # Stats

    stats <- c(job.id, cell.type, me$cis$ntests, me$cis$neqtls, sum(caqtls$FDR < ALPHA))
    
    row.out.df <- as.data.frame(t(stats), stringsAsFactors = FALSE)
    colnames(row.out.df) <- c('JOB_ID', 'CELL_TYPE', 'N_TESTS', 'N_CAQTLS_NOM', 'N_CAQTLS_FDR')

    # Check if the file exists
    if (file.exists(STATS.FILE)) {
        
        # If the file exists, append the new row
        write.table(row.out.df, file = STATS.FILE, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
        
    } else {
        
        # If the file does not exist, create it with the new row
        write.table(row.out.df, file = STATS.FILE, sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE, quote = FALSE)
        
    }
    
    
    
    sessionInfo()
}, error = handle_error)

cat("Script completed successfully.\n")
