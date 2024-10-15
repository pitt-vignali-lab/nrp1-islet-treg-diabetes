# ----------------------------------------------------------------------------
# Setup | Build integrated Treg and Tconv objects
#
# Paper    : Treg insufficiency in autoimmune diabetes is driven by selective
#            loss of Neuropilin-1 on intra-islet regulatory T cells
# Original : 2.Script/3.1.datScVIInit.R
# Input    : 1.Data/revision/datNmf.rds
# Output   : 1.Data/revision/scVITreg.rds
#            1.Data/revision/scVITconv.rds
# Upstream : 04_nmf_cluster_with_foxp3_cells.R
# Note     : Calls scvi-tools through reticulate and sceasy
#            needs a conda env with scvi-tools
# ----------------------------------------------------------------------------

if(!dir.exists("3.Output/revision/3.ScVI")){
  dir.create("3.Output/revision/3.ScVI")
}

library(Seurat)
library(tidyverse)
library(reticulate)
library(sceasy)

sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)

datAll <- readRDS("1.Data/revision/datNmf.rds")
table(datAll$cellType)

dat <- subset(datAll, cellType == "Tconv")

top2000 <- head(VariableFeatures(dat), 2000)

datSub <- dat[top2000]

adata <- convertFormat(datSub, 
                       from="seurat", 
                       to="anndata", 
                       main_layer="counts", 
                       drop_single_values=FALSE)
print(adata) # Note generally in Python, dataset conventions are obs x var

# run setup_anndata
scvi$model$SCVI$setup_anndata(adata) # , batch_key = 'datset'

# create the model
model = scvi$model$SCVI(adata)

# train the model
model$train()

# get the latent represenation
latent = model$get_latent_representation()

# put it back in our original Seurat object
latent <- as.matrix(latent)
rownames(latent) = colnames(datSub)
dat[["scvi"]] <- CreateDimReducObject(embeddings = latent, 
                                      key = "scvi_", 
                                      assay = DefaultAssay(dat))

# Find clusters, then run UMAP, and visualize
dat <- FindNeighbors(dat, dims = 1:10, reduction = "scvi")

dat <- RunUMAP(dat, dims = 1:10, reduction = "scvi", n.components = 2, umap.method = "umap-learn")


dat <- AddMetaData(dat, predictions['predicted.id'])

# dat <- FindClusters(dat, resolution = 0.3, algorithm = 4)

DimPlot(dat, reduction = "umap", label = T, pt.size = 5, group.by = 'predicted.id') +
  coord_fixed() +
  NoLegend()

DimPlot(dat, reduction = "umap", label = T, group.by = 'datset') +
  coord_fixed() +
  NoLegend()


saveRDS(dat, '1.Data/revision/scVITreg.rds')

# cell remove for tconve
umap <- Embeddings(dat, 'umap')

DimPlot(dat, reduction = "umap", label = T) +
  coord_fixed() +
  NoLegend() +
  geom_vline(xintercept = 6.5)

dat <- subset(dat, 
              cells = rownames(umap %>% as.data.frame %>% filter(UMAP_1 < 6.5)))
dat <- FindNeighbors(dat, dims = 1:10, reduction = "scvi")

dat <- RunUMAP(dat, dims = 1:10, reduction = "scvi", n.components = 2)

DimPlot(dat, reduction = "umap", label = T, pt.size = 5, group.by = 'predicted.id') +
  coord_fixed() +
  NoLegend()

saveRDS(dat, '1.Data/revision/scVITconv.rds')