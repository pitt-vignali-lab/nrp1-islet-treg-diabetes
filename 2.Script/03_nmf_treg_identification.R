# ----------------------------------------------------------------------------
# Setup | Identify Tregs by non-negative matrix factorisation
#
# Paper    : Treg insufficiency in autoimmune diabetes is driven by selective
#            loss of Neuropilin-1 on intra-islet regulatory T cells
# Original : 2.Script/2.1.clusterLabNMF(NMF define Treg).R
# Input    : 1.Data/revision/datClean.RDS
# Output   : 1.Data/revision/datNmf.rds
#            3.Output/revision/2.NMFCluster/
# Upstream : 02_qc_and_cleaning.R
# Note     : NMF 0.23.0
#            separates Tregs from the conventional T cell mixture
# ----------------------------------------------------------------------------

if(!dir.exists("3.Output/revision/2.NMFCluster")){
  dir.create("3.Output/revision/2.NMFCluster")
}

library(Seurat)
library(msigdbr)
library(fgsea)
library(tidyverse)
library(ggpattern)
library(NMF)

dat <- readRDS("1.Data/revision/datClean.RDS")

# ndLN ----

datLn <- subset(dat, tissue == "ndLN")

datLn=CreateSeuratObject(
  counts = datLn@assays$RNA@counts,
  meta.data = datLn@meta.data)

datLn <- NormalizeData(datLn) %>% 
  FindVariableFeatures() %>% 
  ScaleData(do.center = F)

vmLn <- datLn@assays$RNA@scale.data

n = 4

resLn <- nmf(vmLn, rank = n) 

fsLn <- extractFeatures(resLn, 10L)

fsLn <- lapply(fsLn,function(x)rownames(resLn)[x])

fsLn

datLn$cluster <- apply(NMF::coefficients(resLn)[1:n,],2,which.max)

Idents(datLn) <- "cluster"

write.csv(datLn@meta.data["cluster"], "3.Output/revision/2.NMFCluster/datLnNMFCluster.csv")

saveRDS(fsLn, "3.Output/revision/2.NMFCluster/datLnNMFClusterMarkers.RDS")

# Islet ----

datIs <- subset(dat, tissue == "islet")

datIs=CreateSeuratObject(
  counts = datIs@assays$RNA@counts,
  meta.data = datIs@meta.data)

datIs <- NormalizeData(datIs) %>% 
  FindVariableFeatures() %>% 
  ScaleData(do.center = F)

vmIs <- datIs@assays$RNA@scale.data

n = 5

resIs <- nmf(vmIs, rank = n) 

fsIs <- extractFeatures(resIs, 50L)

fsIs <- lapply(fsIs,function(x)rownames(resIs)[x])

fsIs[grep("Foxp3", fsIs)]

datIs$cluster <- apply(NMF::coefficients(resIs)[1:n,],2,which.max)

write.csv(datIs@meta.data["cluster"], "3.Output/revision/2.NMFCluster/datIsNMFCluster.csv")

saveRDS(fsIs, "3.Output/revision/2.NMFCluster/datIsNMFClusterMarkers.RDS")
