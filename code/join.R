input.stats <- read.table("data/datasets/hca_brain-organoids_processed/chromatin_accessibility/peaks_stats.tsv", sep = '\t', header = TRUE)
input.stats$key <- paste0(caqtls.stats$RUN, "_", caqtls.stats$CELL_TYPE)

caqtls.stats <- read.table("results/caqtls_stats.tsv", sep = '\t', header = TRUE)
caqtls.stats$key <- paste0(caqtls.stats$RUN, "_", caqtls.stats$CELL_TYPE)

library(dplyr)

joined_df <- input.stats %>%
    inner_join(caqtls.stats, by = "key")

write.table(joined_df, file = "results/stats_joined.tsv", sep = "\t", append = FALSE, quote = FALSE)
