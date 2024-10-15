# ----------------------------------------------------------------------------
# Setup | NMF clustering including Foxp3-positive cells
#
# Paper    : Treg insufficiency in autoimmune diabetes is driven by selective
#            loss of Neuropilin-1 on intra-islet regulatory T cells
# Original : 2.Script/2.2.NMFclusterIncludeFoxp3PostiveCells.R
# Input    : 1.Data/revision/datClean.RDS
# Output   : 1.Data/revision/datNmf.rds
#            3.Output/revision/2.NMFCluster/
# Upstream : 03_nmf_treg_identification.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)

dat <- readRDS("1.Data/revision/datClean.RDS")

nmfIs <- read.csv("3.Output/revision/2.NMFCluster/datIsNMFCluster.csv")
nmfLn <- read.csv("3.Output/revision/2.NMFCluster/datLnNMFCluster.csv")

nmfIs[nmfIs$cluster == 1,]$cluster <- "Treg"
nmfLn[nmfLn$cluster == 4,]$cluster <- "Treg"

nmfClst <- rbind(nmfIs, nmfLn) %>% 
  column_to_rownames(var = "X") %>% 
  rename(nmfClst = cluster)

dat <- AddMetaData(dat, nmfClst)

DimPlot(dat, group.by = "nmfClst")

DotPlot(dat, features = "Foxp3", group.by = "nmfClst")

dat$foxp3 <- PercentageFeatureSet(dat, features = "Foxp3") 

table(dat$foxp3 == 0)

dat$foxp3 <- dat$foxp3 != 0

dat$cellType <- dat$foxp3 != 0 | dat$nmfClst == "Treg"

dat$cellType <- ifelse(dat$cellType == TRUE, "Treg", "Tconv")

DotPlot(dat, features = "Foxp3", group.by = "cellType")

saveRDS(dat, "1.Data/revision/datNmf.rds")

DimPlot(dat, group.by = "nmfClst", split.by = "cellType") +
  theme_bw() +
  coord_fixed()

# plot for manual 
merge(dat@meta.data[c("nmfClst", "cellType")],
      as.data.frame(Embeddings(dat, reduction = "umap")),
      by = 0) %>% 
  mutate(cellType = replace(cellType, cellType == "Treg", "T[reg]~cell")) %>% 
  mutate(cellType = replace(cellType, cellType == "Tconv", "T[conv]~cell")) %>% 
  mutate(cellType = factor(cellType, 
                           levels = c("T[reg]~cell", "T[conv]~cell"))) %>% 
  rename(`NMF cluster` = nmfClst) %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, col = `NMF cluster`))+
  geom_point(alpha = 0.5, size = 2) +
  facet_grid(~cellType, labeller = label_parsed) +
  coord_fixed() +
  theme_bw() +
  theme(aspect.ratio=1,
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 15, face = 2),
        strip.background = element_blank())

ggsave("3.Output/revision/2.NMFCluster/NMFCluster.png", 
       width = 9, height = 4, units = "in", dpi = 300)

# marker heatmap ----

library(patchwork)

dat@meta.data$cellType <- factor(dat@meta.data$cellType,
                                 levels = c("Treg", "Tconv"))

VlnPlot(dat, 
        features = c("Foxp3", "Ctla4", "Tnfrsf4", "Tigit"), 
        group.by = "cellType", 
        ncol = 2) &
  theme_bw() +
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_text(size = 12, face = 2),
        plot.title = element_text(hjust = 0.5, size = 12, face = 2))

ggsave("3.Output/revision/2.NMFCluster/NMFMakersVln.png", 
       width = 4, height = 4, units = "in", dpi = 300)


DotPlot(dat, 
        features = c("Foxp3", "Ctla4", "Tnfrsf4", "Tigit"), 
        group.by = "cellType") &
  theme_bw() +
  theme(aspect.ratio=1,
        strip.background = element_blank(),
        axis.text = element_text(size = 10, face = 2),
        plot.title = element_text(hjust = 0.5, size = 10, face = 2))

ggsave("3.Output/revision/2.NMFCluster/NMFMakersDot.png", 
       width = 6, height = 4, units = "in", dpi = 300)
