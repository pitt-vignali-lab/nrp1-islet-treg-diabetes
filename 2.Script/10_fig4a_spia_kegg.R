# ----------------------------------------------------------------------------
# Fig 4a | Unbiased signaling pathway impact analysis on KEGG
#
# Paper    : Treg insufficiency in autoimmune diabetes is driven by selective
#            loss of Neuropilin-1 on intra-islet regulatory T cells
# Original : 2.Script/3.5.SPIA.R
# Input    : 3.Output/revision/3.ScVI/pseudoBuld/csv/
#            1.Data/revision/RefData/
# Output   : 3.Output/revision/3.ScVI/spiaBar.png
#            3.Output/revision/3.ScVI/spiaTc_islet_Effector.csv
# Upstream : 09_umap_statistics_pseudobulk.R
# Note     : SPIA 2.44.0
#            pathways with Bonferroni-adjusted pGFWER < 0.05 are significant. Needs a local KEGG data directory at 1.Data/revision/RefData/
# ----------------------------------------------------------------------------

library(Seurat)
library(scRepertoire)
library(ggpubr)
library(presto)
library(SPIA)
library(tidyverse)
library(pathview)

dat <- readRDS('1.Data/revision/scVITreg.rds')

DimPlot(dat, group.by = 'predicted.id') +
  coord_fixed()

datSub <- subset(dat, tissue == 'islet' & predicted.id == 'Effector')

refG <- read.csv("1.Data/revision/RefData/mmuRefGeneGeneIDSymbol.csv")

Idents(datSub) <- 'type'

deG <- wilcoxauc(datSub) %>% 
  filter(group == "L/+") %>% 
  # arrange(desc(logFC)) %>% 
  select(feature, logFC, pval, padj, auc) %>% 
  merge(refG, by.x = "feature", by.y = "Symbol")

ranks <- deG[deG$logFC != 0, ]$logFC
names(ranks) <- deG[deG$logFC != 0, ]$GeneID
ranks <- sort(ranks, decreasing = T)

allG <- deG$GeneID

spiaRes <- spia(de = ranks, all = allG, 
                organism = "mmu", 
                data.dir = "1.Data/revision/RefData/", 
                nB = 2000)

write.csv(spiaRes, '3.Output/revision/3.ScVI/spiaTreg_islet_Effector.csv')

spiaRes <- read.csv('3.Output/revision/3.ScVI/spiaTreg_islet_Effector.csv')
dfSub <- spiaRes[c(14, 20, 40, 42, 57, 62, 86, 101, 160, 16, 19, 23, 31, 47, 48, 122, 123, 155),]

ggplot(dfSub, aes(y = reorder(Name, tA), x = tA, fill = Status)) +
  geom_bar(stat = "identity") +
  xlim(-2, 5) +
  geom_text(aes(label = Name, 
                hjust = ifelse(tA > 0, 1, 0), 
                x = ifelse(tA > 0, -0.1, 0.1)),
            size = 3, fontface = 2) +
  scale_fill_manual(labels = c(expression(bold(Activated~'in'~Rosa^Nrp1)), 
                               expression(bold(Inhibited~'in'~Rosa^Nrp1))), 
                    values = alpha(c("red", 'blue'), 0.8)) +
  coord_fixed(1/6) +
  labs(fill = '') +
  theme_minimal() +
  labs( y = "Pathways", 
        x = "Pathway shifting accumulation", 
        fill = "Pathway status") +
  theme(axis.text.y = element_blank(),
        text = element_text(size = 10, face = 2),
        legend.position = "top")

ggsave("3.Output/revision/3.ScVI/spiaBar.png", width = 8, height = 4, units = "in", dpi = 300, bg = "white")

# single cell inform readin ----
datTc <- readRDS('1.Data/revision/scVITconv.rds')

DimPlot(datTc, group.by = 'predicted.id') +
  coord_fixed()

datSub <- subset(datTc, tissue == 'islet' & predicted.id == 'Effector')

table(datSub$type)

refG <- read.csv("1.Data/revision/RefData/mmuRefGeneGeneIDSymbol.csv")

Idents(datSub) <- 'type'

deG <- wilcoxauc(datSub) %>% 
  filter(group == "L/+") %>% 
  # arrange(desc(logFC)) %>% 
  select(feature, logFC, padj) %>% 
  merge(refG, by.x = "feature", by.y = "Symbol")

ranks <- deG[deG$logFC != 0, ]$logFC
names(ranks) <- deG[deG$logFC != 0, ]$GeneID
ranks <- sort(ranks, decreasing = T)

allG <- deG$GeneID

spiaRes <- spia(de = ranks, all = allG, 
                organism = "mmu", 
                data.dir = "1.Data/revision/RefData/", 
                nB = 2000)

write.csv(spiaRes, '3.Output/revision/3.ScVI/spiaTc_islet_Effector.csv')

dfSub <- spiaRes[spiaRes$Name %in% strsplit("T cell receptor signaling pathway, Cell cycle, Th17 cell differentiation, mTOR signaling pathway, FoxO signaling pathway, Natural killer cell mediated cytotoxicity, Th1 and Th2 cell differentiation, VEGF signaling pathway, PI3K-Akt signaling pathway, Ras signaling pathway, IL-17 signaling pathway, HIF-1 signaling pathway, Apelin signaling pathway, NF-kappa B signaling pathway, Human immunodeficiency virus 1 infection, Apoptosis - multiple species", split = ", ")[[1]],]

ggplot(dfSub, aes(y = reorder(Name, tA), x = tA, fill = Status)) +
  geom_bar(stat = "identity") +
  # xlim(-5, 5) +
  geom_text(aes(label = Name, 
                hjust = ifelse(tA > 0, 1, 0), 
                x = ifelse(tA > 0, -0.1, 0.1)),
            size = 4, fontface = 2) +
  theme_minimal() +
  xlim(-15, 15) +
  coord_fixed(2/4) +
  labs( y = "Pathways", 
        x = "Pathway shifting accumulation", 
        fill = "Pathway status") +
  theme(axis.text.y = element_blank(),
        text = element_text(size = 12, face = 2),
        legend.position = "bottom")


pathview(gene.data = rankR, pathway.id = '04930',
         species = "mmu", out.suffix = "gse16873", kegg.native = T)
