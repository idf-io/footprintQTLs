# Perform QTL calling using matrix-eQTL


### Setup ###


## Environment

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


    if (file.exists('.Rprofile')) {

		# renv
	    source('.Rprofile')

    }


    library(tools)
    library(data.table)
    library(dplyr)
    library(tibble)
    library(stringr)

}, error = function(e) {handle_error(e, quit=TRUE)})



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

    if (length(args) == 10) {
        
        cat("Running as script.\n")

		IN.DIR <- as.character(args[1])
   	    OUT.DIR <- as.character(args[2])
	    OUT.FILE <- as.character(args[3])

		MODE <- as.character(args[4])

		ALPHA <- as.numeric(args[5])
		CIS.DIST <- as.numeric(args[6])
	    SLICE.SIZE <- as.numeric(args[7])

		RANGE.LOWER <- as.numeric(args[8])
		RANGE.UPPER <- as.numeric(args[9])
	
	    FORMAT.CHECK <- as.logical(args[10])


    } else if (length(args) == 0) {
        
        cat('Running manually.\n')

        cell.type <- 'DL-EN'

		IN.DIR <- file.path(MATRIX_EQTL_INPUT_DIR, cell.type, mode)
	    OUT.DIR <- file.path(MATRIX_EQTL_OUTPUT_DIR, cell.type, mode)
	    OUT.FILE <- 'chr9_99819926_G_A%chr9:99819749:99820249:501:*:17810'

	    MODE <- 'peak-tests' # {bulk-tests, single-tests, peak-tests}

		ALPHA <- 0.01 # Aygun-2023 caQTLS: 0.05
		CIS.DIST <- 0 # Aygun-2023 caQTLS: within peak
	    SLICE.SIZE <- 2 # matrix eQTL default

		RANGE.LOWER <- as.numeric(args[7])
		RANGE.UPPER <- as.numeric(args[8])
        
	    FORMAT.CHECK <- TRUE


    } else {
        
        stop('Wrong number of arguments.\n')
        
    }

    OUT.FILE <- file.path(OUT.DIR, 'tests', OUT.FILE)
	STATS.FILE <- file.path(OUT.DIR, 'qtl-testing_stats.tsv')
    

    ## Results
	cat('results\n')

	input.list <- input_loader(input.dir = IN.DIR, mode = MODE)


	cis.qtls.list <- list() # To populate

	for (input in input.list[RANGE.LOWER:RANGE.UPPER]) {

		cis.qtls <- call_qtls_meqtl(gt.file = input$gt.file,
									snp.loc.file = input$snp.loc.file,
									phe.file = input$phe.file,
									peak.loc.file = input$peak.loc.file,
									cov.file = input$cov.file,
									cis.dist = CIS.DIST,
									alpha = ALPHA,
									format.check = FORMAT.CHECK,
									slice.size = SLICE.SIZE)[['cis.qtls']]

		cis.qtls.list <- append(cis.qtls.list, list(cis.qtls))

	}

	cis.qtls.all <- bind_rows(cis.qtls.list)

        
	## Save results
	
	create_parent_dir(OUT.FILE)

	if (MODE == 'bulk-tests') {

		write.table(cis.qtls.all, file = OUT.FILE, sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)


		# Stats

		stats <- c(me$cis$ntests, me$cis$neqtls, sum(cis.qtls$FDR < ALPHA))

		row.out.df <- as.data.frame(t(stats), stringsAsFactors = FALSE)
		colnames(row.out.df) <- c('N_TESTS', 'N_QTLS_NOM', 'N_QTLS_FDR')

		write.table(row.out.df, file = STATS.FILE, sep = "\t", row.names = FALSE, quote = FALSE)
	

	} else {

		write.table(cis.qtls.all, file = OUT.FILE, sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

	}
    

}, error = function(e) {handle_error(e, sys.calls(), quit=TRUE)})
