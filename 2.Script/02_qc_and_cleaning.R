# ----------------------------------------------------------------------------
# Setup | Quality control, filtering, normalisation and dimensionality reduction
#
# Paper    : Treg insufficiency in autoimmune diabetes is driven by selective
#            loss of Neuropilin-1 on intra-islet regulatory T cells
# Original : 2.Script/1.2.datClean.R
# Input    : 1.Data/revision/datProcessed.rds
# Output   : 1.Data/revision/datClean.RDS
#            3.Output/revision/1.DatInt/
# Upstream : 01_sample_info_combine.R
# Note     : Seurat 4.1.0 on R 4.1.3
# ----------------------------------------------------------------------------

library(Seurat)
library(harmony)
library(tidyverse)
library(clustree)

# replicate1 <- readRDS("1.Data/revision/set1/set1KmDemultiTcr.rds")
dat <- readRDS("1.Data/revision/datProcessed.rds")

exSet <- read.csv("1.Data/revision/experimantalSettings.csv")
exSet$hashTag <- gsub("-", "_", exSet$hashTag)

# replicate1$replicate = "set1"
# replicate2$replicate = "set2"
# merge(replicate1, replicate2, add.cell.ids = c("set1", "set2"))
dat$HTO_classification <- gsub("-", "_", dat$HTO_classification)
dat <- subset(dat, HTO_classification %in% exSet$hashTag)

dat$hashTag <- dat$HTO_classification

exSet <- dat@meta.data["hashTag"] %>%
  rownames_to_column("cells") %>%
  left_join(exSet, by = "hashTag") %>%
  column_to_rownames("cells") %>%
  select(-1)

dat <- AddMetaData(dat, exSet)

nrow(dat@meta.data)

table(dat$type, dat$tissue)

rm(list = ls()[which(ls() != "dat")])

DefaultAssay(dat) <- "RNA"

dat$log10GenesPerUMI <- log10(dat$nFeature_RNA) / log10(dat$nCount_RNA)
dat$mitoRatio <- PercentageFeatureSet(object = dat, pattern = "^mt-")/100
dat$ribonRatio <- PercentageFeatureSet(dat, pattern = "^Rp[sl]")/100
dat$hbRatio <- PercentageFeatureSet(dat, pattern = "^Hb[^(p)]")/100

VlnPlot(dat, features = c("nFeature_RNA")) +
  scale_y_continuous(trans='log10') +
  geom_hline(yintercept = c(700))  
  
VlnPlot(dat, features = c("nCount_RNA")) +
  scale_y_continuous(trans='log10') +
  geom_hline(yintercept = c(1200)) 



dat@meta.data %>% 
  ggplot(aes(x=replicate, fill=replicate)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

ggsave("3.Output/revision/1.DatInt/Cell_number_distr_samples.png", 
       width = 8, height = 5, units = "in")


dat@meta.data %>% 
  ggplot(aes(x=nCount_RNA, fill = replicate, col = replicate)) + 
  geom_density(alpha = 0.2) + 
  # facet_grid(~type) +
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  expand_limits(x = 0) +
  geom_vline(xintercept = c(300)) + 
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nCount_RNA")
ggsave("3.Output/revision/1.DatInt/Cell_UMIs_distr_reps.png", 
       width = 8, height = 5, units = "in")


dat@meta.data %>% 
  ggplot(aes(x=nFeature_RNA, fill = replicate, col = replicate)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = c(200)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nFeature_RNA")
ggsave("3.Output/revision/1.DatInt/Genes_per_cell_distr_reps.png", 
       width = 8, height = 5, units = "in")

dat@meta.data %>% 
  ggplot(aes(y=log10(nFeature_RNA), x= replicate, fill= replicate)) + 
  geom_boxplot(alpha = 0.2) + 
  # scale_x_log10() + 
  theme_classic() +
  ylab("lgo10(nGene)") +
  geom_vline(xintercept = 300) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nFeature_RNA")

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs ----
dat@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = c(550, 10000)) +
  geom_hline(yintercept = c(400, 3000)) +
  #facet_wrap(~replicate) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nCount_RNA vs. nFeature_RNA")

dat@meta.data %>% 
  ggplot(aes(color=replicate, x=ribonRatio, fill=replicate)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.15) +
  theme(plot.title = element_text(hjust=0.2, face="bold")) +
  ggtitle("Ribosome counts ratio")

dat@meta.data %>% 
  ggplot(aes(color=replicate, x=hbRatio, fill=replicate)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 9e-4) +
  theme(plot.title = element_text(hjust=0.2, face="bold")) +
  ggtitle("Hemoglobin counts ratio")

# Visualize the distribution of mitochondrial gene expression detected per cell ----

dat@meta.data %>% 
  ggplot(aes(color=replicate, x=mitoRatio, fill=replicate)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Mitochondrial counts ratio")
ggsave("3.Output/revision/1.DatInt/Mitochondrial_gene_per_cell_reps.png", 
       width = 8, height = 5, units = "in")

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI ----
dat@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = replicate, fill=replicate)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.80) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Complexity")


# This part is to evaluate the low nCount_RNA cells populations, and it seems like these cells are none-immune cells, but just some nomal pancrece cells ----
# dat$test <- ifelse(dat$nCount_RNA < 1200, "low", "high")
# Idents(dat) <- "test"
# mks <- FindMarkers(dat, ident.1 = "low", ident.2 = "high")
# mksSig <- mks[mks$p_val_adj < 0.05,]
# rownames(mksSig[mksSig$avg_log2FC > 0,])
# write.csv(mksSig, "3.Output/revision/1.DatInt/nCountLowVSHigh.csv")


# Filter out low quality reads using selected thresholds - these will change with experiment ----
datFiltered <- subset(x = dat, 
                      subset= nCount_RNA >= 300 & #nCount_RNA <= 18000 &
                        nFeature_RNA >= 200 & #nFeature_RNA <= 4500 &
                        log10GenesPerUMI > 0.80 &
                        mitoRatio < 0.2)

# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = datFiltered, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keepGenes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filteredCts <- counts[keepGenes, ]

# Reassign to filtered Seurat object
datFiltered <- datFiltered[, colnames(datFiltered) %in% colnames(filteredCts)]

table(datFiltered$replicate)

rm(counts, nonzero, keepGenes, filteredCts)

datFiltered <- NormalizeData(datFiltered)  

# Select the top 1000 most variable features
datFiltered <- FindVariableFeatures(datFiltered, selection.method = "mean.var.plot") 

datFiltered@assays$RNA@data@x[is.na(datFiltered@assays$RNA@data@x)] <- 0

# Scaling RNA data, we only scale the variable features here for efficiency
datFiltered <- ScaleData(datFiltered, features = VariableFeatures(datFiltered)) 

# Run PCA
datFiltered <- RunPCA(datFiltered, features = VariableFeatures(datFiltered))

# datFiltered <- RunHarmony(datFiltered, group.by.vars = "replicate")

DimPlot(datFiltered)

ElbowPlot(datFiltered, ndims = 50) 

# Determine percent of variation associated with each PC
pct <- datFiltered[["pca"]]@stdev / sum(datFiltered[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# Minimum of the two calculation
pcs <- min(co1, co2)
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))
# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()


PCAPlot(datFiltered)

PCAPlot(datFiltered , split.by = "replicate", ncol = 2)

datFiltered <- RunUMAP(datFiltered, 
                       dims = 1:pcs,
                       reduction = "pca")

# sctransform 
datFilteredTest <- PercentageFeatureSet(datFiltered, features = "Nrp1", col.name = "Nrp1")
datFilteredTest <- SCTransform(datFilteredTest, vars.to.regress = "Nrp1", method = "glmGamPoi")


FeaturePlot(datFiltered, c("Nrp1", "Foxp3"), split.by = "type")

datFiltered <- FindNeighbors(datFiltered, 
                             dims = 1:pcs,
                             annoy.metric = "manhattan")

datFiltered <- FindClusters(datFiltered, 
                            resolution = 0.35)

# clustree(datFiltered, prefix = "RNA_snn_res.")
# 
# Idents(datFiltered) <- "RNA_snn_res.0.35"

DimPlot(datFiltered,  label = T, reduction = "umap", pt.size = 1) + NoLegend()

DimPlot(datFiltered,  label = T, reduction = "umap", split.by = "replicate", pt.size = 1) + NoLegend()

DimPlot(datFiltered,  label = T, reduction = "umap", split.by = "tissue", pt.size = 1) + NoLegend()

DimPlot(datFiltered,  label = T, reduction = "umap", split.by = "type", pt.size = 1) + NoLegend()

DimPlot(datFiltered,  label = T, reduction = "umap", pt.size = 1, split.by = "RNA_snn_res.0.3")  + theme_light() + NoLegend()

saveRDS(datFiltered, "1.Data/revision/datClean.RDS")
