# ----------------------------------------------------------------------------
# Setup | Combine per-sample experimental metadata with the CellRanger matrices
#
# Paper    : Treg insufficiency in autoimmune diabetes is driven by selective
#            loss of Neuropilin-1 on intra-islet regulatory T cells
# Original : 2.Script/1.1.inforCombn.R
# Input    : 1.Data/revision/filtered_feature_bc_matrix/
#            1.Data/revision/experimantalSettings.csv
# Output   : 1.Data/revision/datProcessed.rds
# Note     : CellRanger 6.0.1, mm10
#            samples resolved from cell hashes with CITE-seq-Count 1.4.3
# ----------------------------------------------------------------------------

if(!dir.exists('3.Output/revision/1.DatInt')){
  dir.create('3.Output/revision/1.DatInt')
}

library(Seurat)
library(tidyverse)
library(clustree)


# Gene expression matrix
rawDat <- Read10X("1.Data/revision/filtered_feature_bc_matrix/")
colnames(rawDat) <- gsub("-1", "", colnames(rawDat))

# hashtag matix
rawCite <- Read10X("1.Data/revision/umi_count/", gene.column =1)
rawCite <- rawCite[grep("m_", rownames(rawCite))[1:10],]

# Select cell barcodes detected by both RNA and HTO In the example datasets we have already
# filtered the cells for you, but perform this step for clarity.
cellsUse <- intersect(colnames(rawDat), colnames(rawCite))

# Subset RNA and HTO counts by joint cell barcodes
rawDat <- rawDat[, cellsUse]
rawCite <- as.matrix(rawCite[, cellsUse])

# Confirm that the HTO have the correct names
rownames(rawCite) <- paste0("HTO", "_", 1:10) # look back check code

rownames(rawCite)

# Setup Seurat object
datTag <- CreateSeuratObject(counts = rawDat)

# Normalize RNA data with log normalization
datTag <- NormalizeData(datTag)

# Find and scale variable features
datTag <- FindVariableFeatures(datTag, selection.method = "mean.var.plot")
datTag <- ScaleData(datTag, features = VariableFeatures(datTag))

# Add HTO data as a new assay independent from RNA
datTag[["HTO"]] <- CreateAssayObject(counts = rawCite)

# Demultiplexing 
datTag <- NormalizeData(datTag, assay = "HTO", normalization.method="CLR")

datTag <- HTODemux(datTag, assay = "HTO", positive.quantile = 0.1)

table(datTag$HTO_classification.global)

dat <- subset(datTag, idents = "Doublet")

# tcr data
tcr <- read.csv(paste0("1.Data/revision/TCR/filtered_contig_annotations.csv")) %>% 
  mutate(barcode = gsub("-1", "", barcode),
         barcode = gsub(' ', '', barcode),
         cdr3 = gsub(' ', '', cdr3)) %>% 
  filter(chain == 'TRB',
         high_confidence == 'true',
         is_cell == 'true',
         full_length == 'true',
         productive == 'true') %>% 
  select(barcode, cdr3) %>% 
  unique() 

redunTcr <- names(table(tcr$barcode)[table(tcr$barcode) != 1])

tcrSub <- tcr %>% filter(! barcode %in% redunTcr)

dat <- AddMetaData(object = dat, metadata = tcrSub)

# Nrp1 protein citeseq data 

proCite <- Read10X(paste0("1.Data/revision/Nrp1umi_count/"), 
                   gene.column = 1) 
rownames(proCite)[1:2] <- c("antiNrp1", "antiCtrl")
proCite <- data.frame(t(as.matrix(proCite)))

dat <- AddMetaData(object = dat, metadata = proCite)

saveRDS(dat, "1.Data/revision/datProcessed.rds")

