## Setup

setwd('/home/fichtner/projects/footprintQTL/')

library(renv)
renv::activate()

library(MatrixEQTL)
library(data.table)
library(dplyr)
library(tibble)

source("code/helpers/helpers.R")

## Variables

cell_type = "Midbrain EN"
cell_type = ct.format(cell_type)

## User parameters

MODEL <-  modelLINEAR
ALPHA <-  1e-2
ERROR.COV <-  numeric()
CIS.DIST <-  100000 # As in Cuomo-2021

## Variables

DATA.DIR <-  "data/datasets/hca_brain-organoids_processed/"
RESULTS.DIR <-  "results/"

SNP.GT.FILE <-  paste0(DATA.DIR, "covariates/genotype_NA.tsv")
SNP.LOC.FILE <-  paste0(DATA.DIR, "covariates/snp_locations.tsv")
PEAKS.CA.FILE <- paste0(DATA.DIR, "chromatin_accessibility/peak-matrix_rna-qc-cells_norm-reads-in-tss.tsv")
PEAKS.LOC.FILE <- paste0(DATA.DIR, "covariates/peak_locations.tsv")

OUT.FILE <- paste0(RESULTS.DIR, "peak_caQTLs.txt")

GT.PCS.FILE <- paste0(DATA.DIR, "covariates/genotype_pcs.tsv")
CA.PCS.FILE <- paste0(DATA.DIR, "covariates/ca_pcs.tsv")

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
    print("Same headers, different order!")
} else {
    print("Headers don't match! Check order and equivalent sets.")
}

cvrts <- rbind(gt.pcs[snp.gt.header], ca.pcs[snp.gt.header]) %>%
    rownames_to_column(var = 'id')

fwrite(cvrts, paste0(DATA.DIR, "covariates/covariates.tsv"), sep='\t')

CVRT.FILE <- paste0(DATA.DIR, "covariates/covariates.tsv")

# Sanity check: Ensure correct rows(indexes) and columns
cvrt.header <- colnames(fread(CVRT.FILE, nrows=1, header=TRUE))[-1]

if (all(snp.gt.header == peaks.ca.header)) {
    print("Good.")
} else if (all(sort(snp.gt.header) == sort(peaks.ca.header))) {
    print("Same headers, different order!")
} else {
    print("Headers don't match! Check order and equivalent sets.")
}

if (all(snp.gt.header == cvrt.header)) {
    print("Good.")
} else if (all(sort(snp.gt.header) == sort(cvrt.header))) {
    print("Same headers, different order!")
} else {
    print("Headers don't match! Check order and equivalent sets.")
}

if (all(peaks.ca.header == cvrt.header)) {
    print("Good.")
} else if (all(sort(peaks.ca.header) == sort(cvrt.header))) {
    print("Same headers, different order!")
} else {
    print("Headers don't match! Check order and equivalent sets.")
}


snps.gt.idxs <- fread(SNP.GT.FILE, select = 1, header = TRUE)[[1]]
snps.loc.idxs <- fread(SNP.LOC.FILE, select = 1, header = TRUE)[[1]]

if (!all(snps.gt.idxs %in% snps.loc.idxs)) {
    print("Indexes don't match! Check order and equivalent sets.")
}


peaks.ca.idxs <- fread(PEAKS.CA.FILE, select = 1, header = TRUE)[[1]]
peaks.loc.idxs <- fread(PEAKS.LOC.FILE, select = 1, header = TRUE)[[1]]

if (! all(peaks.ca.idxs %in% peaks.loc.idxs)) {
    print("Indexes don't match! Check order and equivalent sets.")
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
    output_file_name.cis = OUT.FILE,
    pvOutputThreshold.cis = ALPHA,
    cisDist = CIS.DIST,
    useModel = MODEL,
    errorCovariance = ERROR.COV,
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected QTLs:', '\n');
show(me$all$eqtls)

# Times
# Midbrain EN - DRD4 - 11 016 cis-tests = 33s
# Added CA PCs as covariate: cis-tests = 44s


## Plot the histogram of all p-values

plot(me)
