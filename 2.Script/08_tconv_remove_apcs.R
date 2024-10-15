# ----------------------------------------------------------------------------
# Setup | Remove antigen-presenting cell contamination from the Tconv object
#
# Paper    : Treg insufficiency in autoimmune diabetes is driven by selective
#            loss of Neuropilin-1 on intra-islet regulatory T cells
# Original : 2.Script/3.3.1.refineTcovneRemoveAPCs.R
# Input    : 1.Data/revision/scVITconv.rds
# Output   : 1.Data/revision/scVITconvNoAPCs.rds
# Upstream : 07_cluster_labeling.R
# ----------------------------------------------------------------------------

library(Seurat)
library(COSG)
library(presto)

dat <- readRDS('1.Data/revision/scVITconv.rds')

dat <- RunPCA(dat)

dat <- FindNeighbors(dat, dims = 1:10)

dat <- RunUMAP(dat, dims = 1:10)

dat <- FindClusters(dat, resolution = 0.2)

DimPlot(dat, label = T)

cosGs <- cosg(dat)

datSub <- subset(dat, idents = c(3, 5:9), invert = T)

datSub <- RunPCA(datSub)

datSub <- RunUMAP(datSub, dims = 1:10)

datSub <- FindClusters(datSub, resolution = 0.15)

DimPlot(datSub, label = T) + NoLegend()

cosGs <- cosg(datSub)

datSub@meta.data <- datSub@meta.data %>% 
  mutate(predicted.id = ifelse(seurat_clusters == 0, 'Memory',
                               ifelse(seurat_clusters == 1, 'Activated',
                                      ifelse(seurat_clusters == 2, 'Proliferating', 'Effector'
                                             ))))

saveRDS(datSub, '1.Data/revision/scVITconvNoAPCs.rds')

DimPlot(datSub, group.by = 'predicted.id')
DotPlot(datSub, features = 'Ifit1') +
  scale_color_gradient2(high = 'red', low = 'blue', mid = 'white', midpoint = 0)

deG <- wilcoxauc(datSub)

deG %>% 
  filter(group == 1,
         logFC > 0) %>% 
  top_n(50, -padj) %>% 
  arrange(padj)

