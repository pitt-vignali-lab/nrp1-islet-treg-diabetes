# ----------------------------------------------------------------------------
# Analysis | Treg TCR repertoire diversity
#
# Paper    : Treg insufficiency in autoimmune diabetes is driven by selective
#            loss of Neuropilin-1 on intra-islet regulatory T cells
# Original : 2.Script/6.2.TCRTreg.R
# Input    : 1.Data/revision/TCR/filtered_contig_annotations.csv
#            1.Data/revision/scVITreg.rds
# Output   : 3.Output/revision/6.TCR/TregTcr_DiversityAnalysis.png
# Upstream : 05_scvi_integration.R
# Note     : scRepertoire 1.3.5 on full-length productive TCRs from CellRanger vdj
# ----------------------------------------------------------------------------

library(Seurat)
library(scRepertoire)
library(ggpubr)
library(tidyverse)

# single cell inform readin ----
datTc <- readRDS('1.Data/revision/scVITreg.rds')
miloMeta <- readRDS('3.Output/revision/7.MiloR/miloMetaTreg.rds') %>% 
  column_to_rownames('CellID')

datTc <- AddMetaData(datTc, miloMeta)

cells <- data.frame(genotype = datTc$type,
                    tissue = datTc$tissue,
                    celltype = datTc$predicted.id,
                    miloType = datTc$Nhood.Group,
                    sampletype = 'Treg')

# tcr data readin ----
tcrDf <- read.csv('1.Data/revision/TCR/filtered_contig_annotations.csv') %>% 
  mutate(barcode = gsub('-1', '', barcode)) %>% 
  merge(cells, by.x = 'barcode', by.y = 0)

table(table(tcrDf$barcode))

tcrDf$GenoCell <- paste(tcrDf$genotype, tcrDf$celltype, sep = "_")

contig_list <- split(tcrDf, f = tcrDf$miloType)

combined <- combineTCR(contig_list, 
                       samples = names(contig_list), 
                       ID = rep(c("L/+", 'WT'), 2))




daDf <- clonalDiversity(combined, 
                cloneCall = "gene", 
                group.by = "sample", 
                x.axis = "ID", 
                n.boots = 100,
                exportTable = T)

daDf %>% as.data.frame() %>% 
  select(c(shannon, sample, ID)) %>% 
  pivot_longer(!c(sample, ID), names_to = "index", values_to = "count") %>% 
  mutate(ID = factor(ID, levels = c('WT', 'L/+')),
         sample = gsub('WT|L/\\+', '', sample)) %>% 
  # separate(sample, into = c('Geno', 'Cell'), sep = '_') %>% 
  # mutate(Cell = factor(Cell, levels = c("Naïve",
  #                                       'Memory',
  #                                       'Effector',
  #                                       'Proliferating',
  #                                       "Antigen \npresentation")),
  #        Geno = factor(Geno, levels = c('WT', 'L/+'))) %>% 
  ggplot(aes(ID, count)) +
  geom_boxplot() +
  geom_point(aes(shape = sample, col = ID), size = 5,  alpha = 0.7) +
  scale_color_manual(values = c('blue', 'red')) +
  guides(color = 'none' ) +
  scale_x_discrete(labels=c("WT" = expression(bold("Cre control")), 
                            "L/+" = expression(bold(Rosa^{Nrp1})))) +
  # facet_wrap(.~index, nrow = 1, scale = 'free_y') +
  theme_bw() +
  labs(x = '', y = 'Index Score') +
  theme(aspect.ratio=2.2,
        legend.position = 'bottom', 
        strip.background = element_blank(),
        text = element_text(size = 15, face = 2))

ggsave('3.Output/revision/6.TCR/TregTcr_DiversityAnalysis_shannon.png', width = 3, height = 6)


# t <- daDf %>% as.data.frame() %>% 
#   select(c(shannon, sample, ID)) %>% 
#   pivot_longer(!c(sample, ID), names_to = "index", values_to = "count")
# 
# 
# t.test(t[t$ID == 'L/+', ]$count, t[t$ID == 'WT', ]$count)
# 
# daDf %>% as.data.frame() %>% 
#   select(-ID) %>% 
#   pivot_longer(!sample, names_to = "index", values_to = "count") %>% 
#   separate(sample, into = c('Geno', 'Cell'), sep = '_') %>% 
#   mutate(Cell = factor(Cell, levels = c("Naïve",
#                                         'Memory',
#                                         'Effector',
#                                         'Proliferating',
#                                         "Antigen \npresentation")),
#          Geno = factor(Geno, levels = c('WT', 'L/+'))) %>% 
#   ggplot(aes(Geno, count)) +
#   geom_boxplot() +
#   geom_point(aes(col = Cell), size = 5) +
#   scale_color_brewer(palette="Dark2") +
#   facet_wrap(.~index, nrow = 1, scale = 'free_y') +
#   theme_bw() +
#   labs(x = '', y = 'Index Score') +
#   theme(aspect.ratio=1,
#         legend.position = 'bottom', 
#         strip.background = element_blank(),
#         text = element_text(size = 12, face = 2))
# 
# ggsave('3.Output/revision/6.TCR/TregTcr_DiversityAnalysis.png', width = 7, height = 5)
# abundanceContig(combined, 
#                 cloneCall = "gene", 
#                 group.by = "sample", 
#                 scale = FALSE, exportTable = T) %>% 
#   separate(values, into = c('Geno', 'Cell', 'Geno2'), sep = '_') %>% 
#   group_by(Geno, Cell) %>% 
#   arrange(Abundance) %>% 
#   mutate(order = 1:n()) %>% 
#   mutate(Cell = factor(Cell, levels = c("Naïve",
#                                         'Memory',
#                                         'Effector',
#                                         'Proliferating',
#                                         "Antigen \npresentation"))) %>% 
#   ggplot(aes(Abundance, order, color = Cell)) +
#   geom_line() +
#   geom_point() +
#   facet_wrap(Geno ~ Cell, scale = 'free_x', nrow = 2) +
#   scale_color_brewer(palette="Dark2") +
#   # coord_fixed(ratio = 1/25) +
#   theme_bw() +
#   scale_y_continuous(trans = 'log10', breaks = seq(0, 600, 100)) +
#   theme(legend.position = 'none',
#         strip.background = element_blank(),
#         text = element_text(size = 12, face = 2))
# 
# clonalHomeostasis(combined, cloneCall = "aa") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))