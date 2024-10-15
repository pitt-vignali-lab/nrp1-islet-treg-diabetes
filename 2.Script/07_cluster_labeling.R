# ----------------------------------------------------------------------------
# Setup | Assign final Treg and Tconv cluster identities
#
# Paper    : Treg insufficiency in autoimmune diabetes is driven by selective
#            loss of Neuropilin-1 on intra-islet regulatory T cells
# Original : 2.Script/3.3.scRNALabel.R
# Input    : 1.Data/revision/scVITconv.rds
# Output   : 1.Data/revision/
# Upstream : 06_annotation_prediction.R
# Note     : Markers from presto wilcoxauc
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(COSG)
library(presto)

# Tconv
dat <- readRDS('1.Data/revision/scVITconv.rds')
datSub <- subset(dat, predicted.id == 'Antigen \npresentation')
DimPlot(datSub, group = 'predicted.id', label = T) +
  coord_fixed() + 
  geom_vline(xintercept = -2)

dimDf <- data.frame(Embeddings(datSub, 'umap'),
                    id =  datSub$predicted.id)%>% 
  mutate(id = ifelse(UMAP_1 < -2, 'Antigen \npresentation', 'Proliferating'), sep = '_')

newMetaDf <- rbind(data.frame(predicted.id = dat$predicted.id) %>% 
                     filter(!rownames(.) %in% rownames(dimDf)),
                   dimDf %>% select(id) %>% rename(predicted.id = 1))

dat <- AddMetaData(dat, newMetaDf)

DimPlot(dat, group = 'predicted.id', label = T) +
  coord_fixed() + 
  geom_vline(xintercept = -2)

Idents(dat) <- 'predicted.id'
mkG <- cosg(dat)
deG <- wilcoxauc(dat) %>% 
  filter(padj < 0.05,
         logFC > 0)

DimPlot(dat, label = T) +
  coord_fixed()

DimPlot(dat, group.by = 'predicted.id', label = T) +
  coord_fixed()


# naive
DotPlot(dat, 
        features = c('Ccr7', 'Tcf7', 'Sell')) + 
  # Effector
  DotPlot(dat, 
          features = c('Ifng', 'Cxcr3')) +
  # Memory
  DotPlot(dat, 
          features = c('Tnfrsf4', 'Cxcr4')) +
  # Antigen
  DotPlot(dat, 
          features = c('Cd74', 'H2-Eb1', 'H2-Ab1')) +
  # Proliferating
  DotPlot(dat, 
          features = c('Top2a', 'Mki67'))

saveRDS(dat, '1.Data/revision/scVITconv.rds')


# Treg 
dat <- readRDS('1.Data/revision/scVITreg.rds')
Idents(dat) <- 'predicted.id'
mkG <- cosg(dat)

DimPlot(dat, label = T) +
  coord_fixed()

# naive
DotPlot(dat, 
        features = c('Ccr7', 'Tcf7', 'Sell')) + 
  # Effector
  DotPlot(dat, 
          features = c('Ifng', 'Cxcr3')) +
  # Memory
  DotPlot(dat, 
          features = c('Tnfrsf4', 'Cxcr4')) 
