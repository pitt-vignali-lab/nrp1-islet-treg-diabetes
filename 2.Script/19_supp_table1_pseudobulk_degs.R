# ----------------------------------------------------------------------------
# Supp Table 1 | Collect significant pseudobulk DEGs into the supplementary table
#
# Paper    : Treg insufficiency in autoimmune diabetes is driven by selective
#            loss of Neuropilin-1 on intra-islet regulatory T cells
# Original : 3.Output/3.ScVI/pseudoBuld/csv/pdfGener.R
# Input    : 3.Output/revision/3.ScVI/pseudoBuld/csv/
# Output   : 3.Output/revision/Supp Table 1. Treg pesudo bulk significant DEGs replicate wised.xlsx
# Upstream : 09_umap_statistics_pseudobulk.R
# Note     : Keeps genes with adjusted P < 0.05, one sheet per comparison
# ----------------------------------------------------------------------------

library(tidyverse)
library(openxlsx)

csvs <- list.files('3.Output/revision/3.ScVI/pseudoBuld/csv/', pattern = 'csv')
t <- read.csv(csvs[1]) %>% 
  arrange(padj)

csvLst <- list()
for (csv in csvs) {
  csvLst[[gsub('.csv', '', csv)]] <- read.csv(csv) %>% 
    filter(padj < 0.05) %>% 
    arrange(padj)
}

write.xlsx(csvLst, 
           file = '3.Output/revision/Supp Table 1. Treg pesudo bulk significant DEGs replicate wised.xlsx')
