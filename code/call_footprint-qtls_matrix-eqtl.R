# Perform QTL calling using matrix-eQTL
#
# Arguments:
#	- GT.FILE (str): genotype file (tsv), snps x samples
#	- SNP.LOC.FILE (str): snp locations file, snp x (snp_id, chr, pos)
# 	- PHE.FILE (str): phenotype file (tsv), peaks x samples
#	- PEAK.LOC.FILE (str): peak locations file (tsv), peaks x (peak_id, chr, start, end)
#	- COV.FILE (str): covariates file (tsv), covariate x samples
#	- OUT.DIR (str)
#	- JOB.ID
#
# Output:
#
# Assumptions:
#
# Notes:
#    - When running manually e.g. within Rstudio, set the PROJ.ROOT, args and some variables manually (lines 9, 41-47, 79-80)


### Setup ###


# Environment

source('code/helpers/R/utils.R')
source('code/helpers/R/matrix-eqtl.R')

tryCatch({



    # Set working dir

    if ('./code' %in% list.dirs(recursive = FALSE)) {
        
        PROJ.ROOT <- '.'
        
    } else if (basename(getwd()) == 'code') {
        
        PROJ.ROOT <- '..'
        
    } else {
        
        cat('Must be run either under project root (/) or /code')
        quit(status = 1)
        
    }


    # renv setup

    library(renv)

    if (file.exists('.Rprofile')) {

	    source('.Rprofile')

    }
    renv::activate()


    ## Libraries
    
    library(tools)
    library(MatrixEQTL)
    library(data.table)
    library(dplyr)
    library(tibble)
    library(stringr)
    library(qqman)

}, error=handle_error)



### Script ###

tryCatch({

    source('code/glob_vars.R') # MATRIX_EQTL_OUTPUT_DIR
    

    ## Args
    
    args <- commandArgs(trailingOnly = TRUE)
    
    # Print args

    for (i in seq_along(args)) {

	    cat('Arg ', i, ': ', args[i], '\n', sep='')

    }
    
    # Register args

    date <- format(Sys.Date(), '%Y-%m-%d')

    if (length(args) == 9) {
        
        cat("Running as script.\n")

        PHE.FILE <- as.character(args[1])
        PEAK.LOC.FILE <- as.character(args[2])
   	
   	GT.FILE <- as.character(args[3])
   	SNP.LOC.FILE <- as.character(args[4])
   	  
	COV.FILE <- as.character(args[5])

   	OUT.DIR <- as.character(args[6])
	
        JOB.ID <- as.character(args[7])

	SKIP.FORMAT.CHECK <- as.logical(args[8])

	SLICE.SIZE <- as.numeric(args[9])

	cat('Arguments supplied:\n')
	cat('Skip input format check:', SKIP.FORMAT.CHECK, '\n')
	cat('Phenotype file:', PHE.FILE, '\n')
	cat('Peak locations file:', PEAK.LOC.FILE, '\n')
	cat('Genotype file:', GT.FILE, '\n')
	cat('Snp locations file:', SNP.LOC.FILE, '\n')
	cat('Covariates file:', COV.FILE, '\n')
	cat('Output dir:', OUT.DIR, '\n')
	cat('Job id:', JOB.ID, '\n')
	cat('Slice size:', SLICE.SIZE, '\n')

        
    } else if (length(args) == 0) {
        
        cat('Running manually.\n')

        cell.type <- 'Neural-progenitors'
	mode <- 'bulk-tests' # {bulk-tests, single-tests}

        PHE.FILE <-  paste0(MATRIX_EQTL_INPUT_DIR, "/", mode, "/", cell.type, "/footprints.tsv")
        PEAK.LOC.FILE <- paste0(MATRIX_EQTL_INPUT_DIR, "/", mode, "/", cell.type, "/peak_locations.tsv")

        GT.FILE <-  paste0(MATRIX_EQTL_INPUT_DIR, "/", mode,  "/", cell.type, "/genotype_NA.tsv")
        SNP.LOC.FILE <-  paste0(MATRIX_EQTL_INPUT_DIR, "/", mode,  "/", cell.type, "/snp_locations.tsv")

        COV.FILE <-  paste0(MATRIX_EQTL_INPUT_DIR, "/", mode, "/", cell.type, "/footprint_pcs.tsv")
    
	OUT.DIR <- paste0(MATRIX_EQTL_INPUT_DIR, "/", mode, "/", cell.type)
        
        JOB.ID <- paste0("manual-run_01234_", date)

	SLICE.SIZE <- 2000 # matrix eQTL default

	cat('Variables:\n')
	cat('Skip input format check:', SKIP.FORMAT.CHECK, '\n')
	cat('Phenotype file:', PHE.FILE, '\n')
	cat('Peak locations file:', PEAK.LOC.FILE, '\n')
	cat('Genotype file:', GT.FILE, '\n')
	cat('Snp locations file:', SNP.LOC.FILE, '\n')
	cat('Covariates file:', COV.FILE, '\n')
	cat('Output dir:', OUT.DIR, '\n')
	cat('Job id:', JOB.ID, '\n')
	cat('Slice size:', SLICE.SIZE, '\n')

        
    } else {
        
        cat('Wrong number of arguments.\n')
        quit(status = 1)
        
    }
    
    STATS.FILE <- paste0(OUT.DIR, '/qtl_test_stats.tsv')
    OUT.FILE.SIMPLE <- paste0(OUT.DIR, '/qtls_simple.tsv')
    OUT.FILE <- paste0(OUT.DIR, '/qtls.tsv')
    

    input_format_check(gt.tsv = GT.FILE,
		       snp.loc.tsv = SNP.LOC.FILE,
		       phe.tsv = PHE.FILE,
		       peak.loc.tsv = PEAK.LOC.FILE,
		       cov.tsv = COV.FILE)


    ## Make IO directories
    
    files <- c(PHE.FILE,
              PEAK.LOC.FILE,
	      #
              GT.FILE,
              SNP.LOC.FILE,
	      COV.FILE,
	      #
    	      STATS.FILE,
              OUT.FILE.SIMPLE,
              OUT.FILE)
    
    for (d in files) {
        
        create_parent_dir(d)
        
    }

    
    ##  User variables
    
    MODEL <-  modelLINEAR
    ALPHA <-  1e-2 # Aygun-2023 caQTLS: 0.05
    ERROR.COV <-  numeric()
    CIS.DIST <-  0 # Aygun-2023 caQTLS: within peak
    

    ## Prep matrix-eQTL input

    gt = SlicedData$new()
    gt$fileDelimiter = "\t"      # the TAB character
    gt$fileOmitCharacters = "NaN" # denote missing values
    gt$fileSkipRows = 1          # one row of column labels
    gt$fileSkipColumns = 1       # one column of row labels
    gt$fileSliceSize = SLICE.SIZE # read file in pieces of 2,000 rows
    gt$LoadFile( GT.FILE )
    
    snp.pos <- read.table(SNP.LOC.FILE, sep = '\t', header = TRUE)
    
    phe = SlicedData$new()
    phe$fileDelimiter = "\t"      # the TAB character
    phe$fileOmitCharacters = "NaN" # denote missing values
    phe$fileSkipRows = 1          # one row of column labels
    phe$fileSkipColumns = 1       # one column of row labels
    phe$fileSliceSize = SLICE.SIZE # read file in pieces of 2,000 rows
    phe$LoadFile( PHE.FILE )
    
    peak.pos <- read.table(PEAK.LOC.FILE, sep = '\t', header = TRUE)
    
    cov = SlicedData$new()
    cov$fileDelimiter = "\t"      # the TAB character
    cov$fileOmitCharacters = "NaN" # denote missing values
    cov$fileSkipRows = 1          # one row of column labels
    cov$fileSkipColumns = 1       # one column of row labels
    cov$fileSliceSize = SLICE.SIZE # read file in pieces of 2,000 rows
    cov$LoadFile( COV.FILE )
    
    
    # Run matrix-eQTL
    
    me = Matrix_eQTL_main(
        snps = gt,
        snpspos = snp.pos,
        gene = phe,
        genepos = peak.pos,
        cvrt = cov,
        output_file_name.cis = OUT.FILE.SIMPLE,
        pvOutputThreshold = 0, # Skip trans-eQTL mapping
        pvOutputThreshold.cis = ALPHA,
        cisDist = CIS.DIST,
        useModel = MODEL,
        errorCovariance = ERROR.COV,
        verbose = TRUE,
        pvalue.hist = TRUE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE)
    
    unlink(OUT.FILE.SIMPLE)
    

    ## Results
    
    # stdout

    cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');

    cat('Nr. of nominal QTLs:\n');
    show(me$cis$neqtls)

    cat('\nNr. of fdr-corrected QTLs:\n');
    show(me$cis$eqtls)
    

    # P-val histogram

    png(paste0(OUT.DIR, '/p-value_histogram.png'))
    
    plot(me)
    
    dev.off()


    # Make QTL df
    cis.qtls <- me$cis$eqtls
    cis.qtls <-  cis.qtls %>%
        mutate(snp_chr = sapply(str_split(snps,'_'), function(z) as.numeric(str_sub(z[1], 4)))) %>% 
        mutate(snp_pos = sapply(str_split(snps,'_'), function(z) as.numeric(z[2]) - 1)) %>% 
        mutate(snp_ref = sapply(str_split(snps,'_'), function(z) z[3])) %>% 
        mutate(snp_alt = sapply(str_split(snps,'_'), function(z) z[4])) %>% 
        mutate(peak_chr = sapply(str_split(gene,':'), function(z) as.numeric(str_sub(z[1], 4)))) %>%
        mutate(peak_start = sapply(str_split(gene,':'), function(z) as.numeric(z[2]) - 1)) %>% 
        mutate(peak_end = sapply(str_split(gene,':'), function(z) as.numeric(z[3]))) %>% 
        mutate(peak_len = sapply(str_split(gene,':'), function(z) as.numeric(z[4]))) %>%
        arrange(desc(FDR)) 
        

    # Manhattan plot nominal
    png(paste0(OUT.DIR, '/manhattan-plot_nominal.png'))
    manhattan(cis.qtls, chr = 'snp_chr', bp = 'snp_pos', p = 'pvalue', snp = 'snps')
    dev.off()
    
    # Manhattan plot FDR
    png(paste0(OUT.DIR, '/manhattan-plot_fdr.png'))
    manhattan(cis.qtls, chr = 'snp_chr', bp = 'snp_pos', p = 'FDR', snp = 'snps')
    dev.off()
    
    # Q-Q plot nominal
    png(paste0(OUT.DIR, '/qq-plot_nominal.png'))
    qq(cis.qtls$pvalue, main = paste0('p-value nominal'))
    dev.off()
    
    # Q-Q plot FDR
    png(paste0(OUT.DIR, '/qq-plot_fdr.png'))
    qq(cis.qtls$FDR, main = paste0('p-value FDR'))
    dev.off()
    

    ## Save results
    
    write.table(cis.qtls, file = OUT.FILE, sep = '\t', quote = FALSE)
    
    # Stats

    stats <- c(JOB.ID, me$cis$ntests, me$cis$neqtls, sum(cis.qtls$FDR < ALPHA))
    
    row.out.df <- as.data.frame(t(stats), stringsAsFactors = FALSE)
    colnames(row.out.df) <- c('JOB_ID', 'N_TESTS', 'N_QTLS_NOM', 'N_QTLS_FDR')

    # Check if the file exists
    if (file.exists(STATS.FILE)) {
        
        # If the file exists, append the new row
        write.table(row.out.df, file = STATS.FILE, sep = '\t', row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
        
    } else {
        
        # If the file does not exist, create it with the new row
        write.table(row.out.df, file = STATS.FILE, sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE, quote = FALSE)
        
    }
    
}, error = handle_error)


sink(paste0(OUT.DIR, 'session_info.txt'))
sessionInfo()
sink()

cat("Script completed successfully.\n")
