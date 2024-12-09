caqtls <- read.csv("data/datasets/hca_brain-organoids_processed/chromatin_accessibility/peaks_stats.tsv", sep = '\t', header = TRUE, row.names = NULL)
caqtls

cell.types <- unique(caqtls$CELL_TYPE)

PROJ.ROOT = "/home/fichtner/projects/footprintQTL"
HCA.BORGS.PROJ = "data/datasets/hca_brain-organoids/"

eqtls"eQTL_mapping/eSNPs_significant_all_celltypes_HVGs.tsv"


# 
# caqtl$ePeak_fraction <- caqtls$N_PEAKS_EQTLS / caqtls$MANUAL_N_CA_QTLS





for (ct in cell.types) {
    
    
    print(ct)
    
}
