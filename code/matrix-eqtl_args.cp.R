## Setup

PROJ.ROOT <- "/home/fichtner/projects/footprintQTL"

# Libraries

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

source("code/helpers/helpers.R")


#  User parameters

# cell.type <- "Midbrain EN"
cell.type <- "Differentiating RG"
local.run <- "sep8"

MODEL <-  modelLINEAR
ALPHA <-  1e-2
ERROR.COV <-  numeric()
CIS.DIST <-  1000 # As Aygun?


# Atgs

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 1) {
    
    cell.type.r <- as.character(args[1])
    run <- "EMPTY_ARGS_RUN"
    
} else if (length(args) == 0) {
    
    
    cell.type.r <- cell.type
    run <- local.run
    
} else if (length(args) == 2) {
    
    cell.type.r <- as.character(args[1])
    run <- as.character(args[2])
    
} else {
    
    print("Wrong number of arguments!")
    stop()
    
}

cell.type = ct.format(cell.type.r)
cell.type.alt = ct.format.alt(cell.type.r)

# Variables

DATA.DIR <-  "data/datasets/hca_brain-organoids_processed/"
RESULTS.DIR <-  "results/"

STATS.FILE <- "results/caqtls_stats.tsv"
OUT.FILE.SIMPLE <- paste0(RESULTS.DIR, "peak_caqtls_simple_", cell.type, ".tsv")
OUT.FILE <- paste0(RESULTS.DIR, "peak_caqtls_", cell.type, ".tsv")

SNP.GT.FILE <-  paste0(DATA.DIR, "covariates/genotype_NA_", cell.type , ".tsv")
SNP.LOC.FILE <-  paste0(DATA.DIR, "covariates/snp_locations_", cell.type, ".tsv")
PEAKS.CA.FILE <- paste0(DATA.DIR, "chromatin_accessibility/peak-matrix_rna-qc-cells_norm-reads-in-tss_", cell.type, ".tsv")
PEAKS.LOC.FILE <- paste0(DATA.DIR, "covariates/peak_locations_", cell.type, ".tsv")

GT.PCS.FILE <- paste0(DATA.DIR, "covariates/genotype_pcs_", cell.type, ".tsv")
CA.PCS.FILE <- paste0(DATA.DIR, "covariates/ca_pcs_", cell.type, ".tsv")







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

fwrite(cvrts, paste0(DATA.DIR, "covariates/covariates_", cell.type, ".tsv"), sep="\t")

CVRT.FILE <- paste0(DATA.DIR, "covariates/covariates_", cell.type, ".tsv")

# Sanity check: Ensure correct rows(indexes) and columns
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

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected QTLs:', '\n');
show(me$cis$eqtls)

# Times
# Midbrain EN - DRD4 - 11 016 cis-tests = 33s
# Added CA PCs as covariate: cis-tests = 44s


## Plot the histogram of all p-values
png(paste0(RESULTS.DIR, "caQTLs_hist_", cell.type, ".png"))

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
as.data.frame(table(gwasResults$CHR))

## Plots

# Manhattan FDR
png(paste0(RESULTS.DIR, "caQTLs_manhattan-plot_", cell.type, ".png"))
manhattan(caqtls, chr = 'snp_chr', bp = 'snp_pos', p = 'FDR', snp = 'snps')
dev.off()

# Manhattan nominal
png(paste0(RESULTS.DIR, "caQTLs_manhattan-plot_nominal_", cell.type, ".png"))
manhattan(caqtls, chr = 'snp_chr', bp = 'snp_pos', p = 'pvalue', snp = 'snps')
dev.off()

# Q-Q nominal
png(paste0(RESULTS.DIR, "caQTLs_qq_plot_nominal_", cell.type, ".png"))
qq(caqtls$pvalue, main = paste0("Q-Q plot:", cell.type, "pvalue"))
dev.off()

# Q-Q FDR
png(paste0(RESULTS.DIR, "caQTLs_qq_plot_fdr_", cell.type, ".png"))
qq(caqtls$FDR, main = paste0("Q-Q plot:", cell.type, "pvalue"))
dev.off()


## Save results

# caQTLs
write.table(caqtls, file = OUT.FILE, sep = "\t", quote = FALSE)

# Stats

stats <- c(run, cell.type, me$cis$ntests, me$cis$neqtls, sum(caqtls$FDR < ALPHA))

row.out.df <- as.data.frame(t(stats), stringsAsFactors = FALSE)
colnames(row.out.df) <- c('RUN', 'CELL_TYPE', 'N_TESTS', 'N_CAQTLS_NOM', 'N_CAQTLS_FDR')

# Check if the file exists
if (file.exists(STATS.FILE)) {
    # If the file exists, append the new row
    write.table(row.out.df, file = STATS.FILE, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
} else {
    # If the file does not exist, create it with the new row
    write.table(row.out.df, file = STATS.FILE, sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE, quote = FALSE)
}



sessionInfo()
