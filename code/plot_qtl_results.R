# Plot results of QTL calling


### Setup ###


## Environment

source('code/helpers/R/utils.R')
source('code/helpers/R/matrix-eqtl.R') # plot_results
source('code/glob_vars.R') # MATRIX_EQTL_OUTPUT_DIR

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


	library(dplyr)
    library(data.table)
    library(qqman)

}, error = function(e) {handle_error(e, sys.calls(), quit=TRUE)}
)



### Script ###

tryCatch({

    

    ## Args
    
    args <- commandArgs(trailingOnly = TRUE)
    
    # Print args

    for (i in seq_along(args)) {

	    cat('Arg ', i, ': ', args[i], '\n', sep='')

    }
    
    # Register args

    if (length(args) == 4) {
        
        cat("Running as script.\n")

		IN.FILE <- as.character(args[1])
   	    OUT.DIR <- as.character(args[2])

		MODE <- as.character(args[3])
		ALPHA <- as.numeric(args[4])


    } else if (length(args) == 0) {
        
        cat('Running manually.\n')

        cell.type <- 'DL-EN'

		IN.FILE <- file.path(MATRIX_EQTL_OUTPUT_DIR, cell.type, mode, 'qtls_all.tsv')
		OUT.DIR <- file.path(MATRIX_EQTL_OUTPUT_DIR, cell.type, mode, 'plots')

	    MODE <- 'peak-tests' # {bulk-tests, single-tests, peak-tests}
		ALPHA <- 0.01 # Aygun-2023 caQTLS: 0.05


    } else {
        
        stop('Wrong number of arguments.\n')
        
    }


	if (MODE == 'bulk-tests') {

		## Load qtls tsv

		out.cols <- c(snp = 'character',
					  gene = 'character',
					  fdr = 'numeric',
					  pvalue = 'numeric',
					  beta = 'numeric',
					  snp_chr = 'numeric',
					  snp_pos = 'numeric',
					  snp_ref = 'character',
					  snp_alt = 'character',
					  peak_chr = 'numeric',
					  peak_start = 'numeric',
					  peak_end = 'numeric',
					  peak_len = 'numeric')

		qtls.df <- fread(IN.FILE, sep='\t', header = FALSE, col.names = names(out.cols))


		## Metadata

		n.qtls.nom <- qtls.df %>% filter(pvalue < ALPHA) %>% nrow
		n.qtls.fdr <- qtls.df %>% filter(fdr < ALPHA) %>% nrow

		stats.file <- file.path(dirname(IN.FILE), 'qtls_stats.tsv')
		lines <- c('\t\tn_qtls',
					paste0('nominal_pval\t', n.qtls.nom),
					paste0('fdr_pval\t', n.qtls.fdr))

		print(lines)
		writeLines(lines, stats.file)


		## Plot

		plot_results(qtls.df = qtls.df,
					 out.dir = OUT.DIR,
					 alpha = ALPHA)

	} else if (MODE %in% c('single-tests', 'peak-tests')) {

		## Load qtls tsv

		out.cols <- c(snp = 'character',
					  gene = 'character',
					  fdr_meqtl = 'numeric',
					  pvalue = 'numeric',
					  beta = 'numeric',
					  snp_chr = 'numeric',
					  snp_pos = 'numeric',
					  snp_ref = 'character',
					  snp_alt = 'character',
					  peak_chr = 'numeric',
					  peak_start = 'numeric',
					  peak_end = 'numeric',
					  peak_len = 'numeric')

		qtls.df <- fread(IN.FILE, sep='\t', header = FALSE, col.names = names(out.cols))


		## FDR

		qtls.df$fdr <- p.adjust(qtls.df[['pvalue']], 'BH')
		qtls.df <- qtls.df %>% filter(snp != 'placeholder_snp') # remove placeholder loci

		new.df.file <- file.path(dirname(IN.FILE), 'qtls_all_fdr.tsv')
		write.table(qtls.df, file = new.df.file, sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
		

		## Metadata

		n.qtls.nom <- qtls.df %>% filter(pvalue < ALPHA) %>% nrow
		n.qtls.fdr <- qtls.df %>% filter(fdr < ALPHA) %>% nrow

		stats.file <- file.path(dirname(IN.FILE), 'qtls_stats.tsv')
		lines <- c('\t\tn_qtls',
					paste0('nominal_pval\t', n.qtls.nom),
					paste0('fdr_pval\t', n.qtls.fdr))

		print(lines)
		writeLines(lines, stats.file)


		## Plot

		plot_results(qtls.df = qtls.df,
					 out.dir = OUT.DIR,
					 alpha = ALPHA)

	}
    
}, error = function(e) {handle_error(e, sys.calls(), quit=TRUE)}
)
