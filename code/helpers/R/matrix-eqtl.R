.args <- commandArgs(trailingOnly = FALSE)
.script_path <- normalizePath(sub('--file=', '', .args[grep('--file=', .args)]))
.script_dir <- dirname(.script_path)

source(file.path(.script_dir, '/helpers/R/utils.R'))

library(ggplot2)
library(MatrixEQTL)
#library(dplyr)
#library(data.table)



input_format_check <- function(gt.tsv, snp.loc.tsv, phe.tsv, peak.loc.tsv, cov.tsv) {
    # Check the tabular format of the input for correctness.
    #
    # Notes:
    #    - columns = sample labels
    #    - rows = snp & peak labels
    #
    # Todos:
    #    - Check for missing values and if they are the same as the user specified and expected

    # Get column and row labels

	tryCatch({

		gt.header <- colnames(fread(gt.tsv, nrows = 1, header=TRUE))[-1]
		phe.header <- colnames(fread(phe.tsv, nrows=1, header=TRUE))[-1]
		cov.header <- colnames(fread(cov.tsv, nrows=1, header=TRUE))[-1]

		gt.idxs <- fread(gt.tsv, select = 1, header = TRUE)[[1]]
		snp.loc.idxs <- fread(snp.loc.tsv, select = 1, header = TRUE)[[1]]

		phe.idxs <- fread(phe.tsv, select = 1, header = TRUE)[[1]]
		peak.loc.idxs <- fread(peak.loc.tsv, select = 1, header = TRUE)[[1]]
		

		# Columns are unique

		if (length(gt.header) != length(unique(gt.header))) {
			stop('FORMAT ERROR: genotype file columns are note unique.')
		}

		if (length(phe.header) != length(unique(phe.header))) {
			stop('FORMAT ERROR: phenotype file columns are note unique.')
		}

		if (length(cov.header) != length(unique(cov.header))) {
			stop('FORMAT ERROR: covariate file columns are note unique.')
		}


		# Columns are identical and have the same order
		
		if (all(gt.header == phe.header)) {
			
			print('CORRECT: Same columns and their order.')
			
		} else if (all(sort(gt.header) == sort(phe.header))) {
			
			stop('FORMAT ERROR: Genotype and phenotype files have same headers but different order!')
			
		} else {
			
		cat(gt.header, '\n')
		cat(phe.header, '\n')
			stop('FORMAT ERROR: Genotype and phenotype file headers don\'t match!')
			
		}
		

		if (all(gt.header == cov.header)) {
			
			print('CORRECT: Same columns and their order.')
			
		} else if (all(sort(gt.header) == sort(cov.header))) {
			
			stop('FORMAT ERROR: Genotype and covariate files have same headers but different order!')
			
		} else {
			
			stop('FORMAT ERROR: Genotype and covariate file headers don\'t match!')
			
		}
		

		if (all(phe.header == cov.header)) {
			
			print('CORRECT: Same columns and their order.')
			
		} else if (all(sort(phe.header) == sort(cov.header))) {
			
			stop('FORMAT ERROR: Phenotype and covariate files have same headers but different order!')
			
		} else {
			
			stop('FORMAT ERROR: Phenotype and covariate file headers don\'t match!')
			
		}
		

		# Rows are unique

		if (length(gt.idxs) != length(unique(gt.idxs))) {

			stop('FORMAT ERROR: genotype file rows are note unique.')

		}

		if (length(snp.loc.idxs) != length(unique(snp.loc.idxs))) {

			stop('FORMAT ERROR: snp location file rows are note unique.')

		}

		if (length(phe.idxs) != length(unique(phe.idxs))) {

			stop('FORMAT ERROR: phenotype file rows are note unique.')

		}

		if (length(peak.loc.idxs) != length(unique(peak.loc.idxs))) {

			stop('FORMAT ERROR: peak location file rows are note unique.')

		}


		# Rows are identical 
		
		if (!all(gt.idxs %in% snp.loc.idxs)) {
			
			stop('FORMAT ERROR: Genotype and snp location file indexes (snp labels) don\'t match!')
			
		}
		
		if (! all(phe.idxs %in% peak.loc.idxs)) {
			
			stop('FORMAT ERROR: Phenotype and peak file indexes (peak labels) don\'t match!')
			
		}

	}, error = function(e) {handle_error(e, sys.calls(), quit=TRUE)})

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
        
        print('Colnames match.')
        
    } else if (all(sort(colnames(gt.pcs)) == sort(colnames(phe.pcs)))) {
        
        stop('Same headers/columns, different order!')
        
    } else {
        
        stop('Header don\'t match! Check order and equivalent sets.')
        
    }
    
    cvrts <- rbind(gt.pcs[phe.header], phe.pcs[phe.header]) %>%
        rownames_to_column(var = 'id')
    
    CVRT.FILE <- tsv.out

    return(CVRT.FILE)
    #fwrite(cvrts, CVRT.FILE, sep='\t')

}


call_qtls_meqtl <- function(gt.file,
							snp.loc.file,
							phe.file,
							peak.loc.file,
							cov.file,
							cis.dist,
							alpha,
							format.check = FALSE,
							slice.size = 2000) {
	# Call QTLs using matrix-eQTL

	tryCatch({

		## Checks

		if (format.check == TRUE) {

			input_format_check(gt.tsv = gt.file,
							   snp.loc.tsv = snp.loc.file,
							   phe.tsv = phe.file,
							   peak.loc.tsv = peak.loc.file,
							   cov.tsv = cov.file)

		}
		

		##  User variables
		
		model <-  modelLINEAR
		error.cov <-  numeric()
		

		## Prep matrix-eQTL input

		gt = SlicedData$new()
		gt$fileDelimiter = '\t'      # the TAB character
		gt$fileOmitCharacters = 'NaN' # denote missing values
		gt$fileSkipRows = 1          # one row of column labels
		gt$fileSkipColumns = 1       # one column of row labels
		gt$fileSliceSize = slice.size # read file in pieces of 2,000 rows
		gt$LoadFile( gt.file )
		
		snp.pos <- read.table(snp.loc.file, sep = '\t', header = TRUE)
		
		phe = SlicedData$new()
		phe$fileDelimiter = '\t'      # the TAB character
		phe$fileOmitCharacters = 'NaN' # denote missing values
		phe$fileSkipRows = 1          # one row of column labels
		phe$fileSkipColumns = 1       # one column of row labels
		phe$fileSliceSize = slice.size # read file in pieces of 2,000 rows
		phe$LoadFile( phe.file )
		
		peak.pos <- read.table(peak.loc.file, sep = '\t', header = TRUE)
		
		cov = SlicedData$new()
		cov$fileDelimiter = '\t'      # the TAB character
		cov$fileOmitCharacters = 'NaN' # denote missing values
		cov$fileSkipRows = 1          # one row of column labels
		cov$fileSkipColumns = 1       # one column of row labels
		cov$fileSliceSize = slice.size # read file in pieces of 2,000 rows
		cov$LoadFile( cov.file )
		
		
		## Run matrix-eQTL
		
		me = Matrix_eQTL_main(
			snps = gt,
			snpspos = snp.pos,
			gene = phe,
			genepos = peak.pos,
			cvrt = cov,
			pvOutputThreshold = 0, # Skip trans-eQTL mapping
			pvOutputThreshold.cis = alpha,
			cisDist = cis.dist,
			useModel = model,
			errorCovariance = error.cov,
			verbose = TRUE,
			pvalue.hist = TRUE,
			min.pv.by.genesnp = FALSE,
			noFDRsaveMemory = FALSE)

		cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n------\n');


		## Make QTL dataframe

		cis.qtls <- me$cis$eqtls

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

		if (nrow(cis.qtls) > 0) {

			cis.qtls <-  cis.qtls %>%
				mutate(snp_chr = sapply(str_split(snps,'_'), function(z) as.numeric(str_sub(z[1], 4)))) %>% 
				mutate(snp_pos = sapply(str_split(snps,'_'), function(z) as.numeric(z[2]) - 1)) %>% 
				mutate(snp_ref = sapply(str_split(snps,'_'), function(z) z[3])) %>% 
				mutate(snp_alt = sapply(str_split(snps,'_'), function(z) z[4])) %>% 
				mutate(peak_chr = sapply(str_split(gene,':'), function(z) as.numeric(str_sub(z[1], 4)))) %>%
				mutate(peak_start = sapply(str_split(gene,':'), function(z) as.numeric(z[2]) - 1)) %>% 
				mutate(peak_end = sapply(str_split(gene,':'), function(z) as.numeric(z[3]))) %>% 
				mutate(peak_len = sapply(str_split(gene,':'), function(z) as.numeric(z[4]))) %>%
				rename(snp = snps) %>%
				rename(fdr = FDR) %>%
				arrange(desc(fdr)) %>%
				select(all_of(names(out.cols)))

			print(str(cis.qtls))
			print(cis.qtls)

		} else {

			cis.qtls <- data.frame(lapply(out.cols, function(type) vector(type, length=0)))
			colnames(cis.qtls) <- names(out.cols)

		}
	
		return(list(me = me, cis.qtls = cis.qtls))

	}, error = function(e) {handle_error(e, sys.calls(), quit=TRUE)})

}


input_loader <- function(input.dir, mode) {
	# Load multiple instances of input files as a list of lists.

	tryCatch({

		## Checks

		stopifnot(dir.exists(input.dir))
		stopifnot(mode %in% c('bulk-tests', 'single-tests', 'peak-tests'))


		## Populate input list

		input.list = list() # [[input_file_X, ], ] labeled both inner and outer list

		if (mode == 'bulk-tests') {

			input.list[['bulk']] <- list(gt.file = file.path(input.dir, 'genotype_NA.tsv'),
										 snp.loc.file = file.path(input.dir, 'snp_location.tsv'),
										 phe.file = file.path(input.dir, 'phenotype.tsv'),
										 peak.loc.file = file.path(input.dir, 'peak_location.tsv'),
										 cov.file = file.path(input.dir, 'covariates.tsv'))


		} else if (mode == 'single-tests') {

			snp.peak.pairs <- fread(file.path(input.dir, 'tests_snp_peak_pairs.bed'), sep = '\t', header = FALSE)

			for (i in 1:nrow(snp.peak.pairs)) {

				snp <- as.character(snp.peak.pairs[[i, 4]])
				peak <- as.character(snp.peak.pairs[[i, 8]])

				input.files <- list(gt.file = paste0(input.dir, '/genotypes/genotype_NA%', snp, '.tsv'),
							  		snp.loc.file = paste0(input.dir, '/snp_locations/snp_location%', snp, '.tsv'),
							  		phe.file = paste0(input.dir, '/phenotypes/phenotype%', peak, '.tsv'),
							  		peak.loc.file = paste0(input.dir, '/peak_locations/peak_location%', peak , '.tsv'),
							  		cov.file = paste0(input.dir, '/covariates/covariates%', snp, '%', peak, '.tsv'))

				input.list[[paste0(snp, '_>_', peak)]] <- input.files

			}


		} else if (mode == 'peak-tests') {

			snp.peak.pairs <- fread(file.path(input.dir, 'tests_snp_peak_pairs.bed'), sep = '\t', header = FALSE)
			peaks <- unique(as.character(snp.peak.pairs[[8]]))

			for (peak in peaks) {

				input.files <- list(gt.file = paste0(input.dir, '/genotypes/genotype_NA%', peak, '.tsv'),
							  		snp.loc.file = paste0(input.dir, '/snp_locations/snp_location%', peak, '.tsv'),
							  		phe.file = paste0(input.dir, '/phenotypes/phenotype%', peak, '.tsv'),
							  		peak.loc.file = paste0(input.dir, '/peak_locations/peak_location%', peak , '.tsv'),
							  		cov.file = paste0(input.dir, '/covariates/covariates%', peak, '.tsv'))

				input.list[[peak]] <- input.files

			}

		}
		

		return(input.list)


	}, error = function(e) {handle_error(e, sys.calls(), quit=TRUE)})

}


plot_results <- function(qtls.df, out.dir, alpha) {
    # Plot results of qtl calling data

	tryCatch({

		## Make directories
		
		create_parent_dir(out.dir)


		## Nominal p-val histogram

		p <- ggplot(qtls.df, aes(x = pvalue)) +
				geom_histogram(aes(y = after_stat(density)), bins = 100) +
				geom_hline(yintercept = 1, color = 'blue', linetype = 'dashed', linewidth = 1) +
				theme_classic() +
				labs(y = 'P-values', x = 'Density') +
				theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))

		ggsave(file.path(out.dir, 'histogram_pval_nominal.png'), plot = p, width = 10, height = 6, dpi = 300)


		## FDR p-val histogram

		p <- ggplot(qtls.df, aes(x = fdr)) +
				geom_histogram(aes(y = after_stat(density)), bins = 100) +
				geom_hline(yintercept = 1, color = 'blue', linetype = 'dashed', linewidth = 1) +
				theme_classic() +
				labs(y = 'P-values', x = 'Density') +
				theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))

		ggsave(file.path(out.dir, 'histogram_pval_fdr.png'), plot = p, width = 10, height = 6, dpi = 300)


		## Manhattan plot nominal

		png(file.path(out.dir, 'manhattan-plot_nominal.png'))
		manhattan(qtls.df,
				  chr = 'snp_chr',
				  bp = 'snp_pos',
				  p = 'pvalue',
				  snp = 'snp',
				  suggestiveline = FALSE,
				  genomewideline = -log10(alpha))
		dev.off()

	
		## Manhattan plot FDR
		
		png(file.path(out.dir, 'manhattan-plot_fdr.png'))
		manhattan(qtls.df,
				  chr = 'snp_chr',
				  bp = 'snp_pos',
				  p = 'fdr',
				  snp = 'snp',
				  suggestiveline = FALSE,
				  genomewideline = -log10(alpha))
		dev.off()
	

		## Q-Q plot nominal

		png(file.path(out.dir, 'qq-plot_nominal.png'))
		qq(qtls.df$pvalue, main = 'p-value nominal')
		dev.off()
	

		## Q-Q plot FDR

		png(file.path(out.dir, 'qq-plot_fdr.png'))
		qq(qtls.df$fdr, main = 'p-value FDR')
		dev.off()

	}, error = function(e) {handle_error(e, sys.calls(), quit=TRUE)})

}
