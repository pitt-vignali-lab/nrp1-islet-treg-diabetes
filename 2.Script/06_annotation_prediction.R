# ----------------------------------------------------------------------------
# Setup | Transfer cell-type annotations onto the integrated objects
#
# Paper    : Treg insufficiency in autoimmune diabetes is driven by selective
#            loss of Neuropilin-1 on intra-islet regulatory T cells
# Original : 2.Script/3.2.annotationPrediction.R
# Input    : 1.Data/revision/scVITconv.rds
#            1.Data/initial/scVITconv.rds
# Output   : 1.Data/revision/
# Upstream : 05_scvi_integration.R
# ----------------------------------------------------------------------------

library(simspec)
library(Seurat)

datAll <- readRDS("1.Data/revision/datNmf.rds")
table(datAll$cellType)

datC <- subset(datAll, cellType == "Tconv")
datC <- NormalizeData(datC)
datC <- FindVariableFeatures(datC, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(datC)
datC <- ScaleData(datC, features = all.genes)
datC <- RunPCA(datC, features = VariableFeatures(object = datC))
datC <- FindNeighbors(datC, dims = 1:10)
datC <- RunUMAP(datC, dims = 1:10)

datC$sample <- 'current'

datR <- readRDS('1.Data/initial/scVITconv.rds')

datR$sample <- 'ref'

panc8 <- merge(datC, datR)

pancreas.list <- SplitObject(panc8, split.by = "sample")

for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000,
                                             verbose = FALSE)
}

reference.list <-  SplitObject(datR, split.by = "datset")

pancreas.anchors <- FindIntegrationAnchors(reference.list, dims = 1:30)

pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)

DefaultAssay(pancreas.integrated) <- "integrated"

pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)

DimPlot(pancreas.integrated, reduction = "umap", group.by = "annot")

pancreas.query <- pancreas.list[["current"]]

pancreas.anchors <- FindTransferAnchors(reference = pancreas.integrated, query = pancreas.query,
                                        dims = 1:30, reference.reduction = "pca")

predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.integrated$annot,
                            dims = 1:30)

pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)

DimPlot(pancreas.query, group.by = 'predicted.id')

pancreas.query <- RunPCA(pancreas.query, features = VariableFeatures(object = pancreas.query))
pancreas.query <- RunUMAP(pancreas.query, reduction = "scvi", dims = 1:10)

DimPlot(datC)

datC <- AddMetaData(datC, predictions['predicted.id'])
DimPlot(datC, group.by = 'predicted.id')

# naive
DotPlot(datC, 
        features = c('Ccr7', 'Tcf7', 'Sell'),
        group.by = 'predicted.id') + 
  # Effector
  DotPlot(datC, 
          features = c('Ifng'),
          group.by = 'predicted.id') +
  # Memory
  DotPlot(datC, 
          features = c('Tnfrsf4', 'Cxcr4'),
          group.by = 'predicted.id') 

library(COSG)

mkG <- cosg(datC)
