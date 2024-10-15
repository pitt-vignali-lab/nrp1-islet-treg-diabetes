# ----------------------------------------------------------------------------
# Analysis | Tconv TCR repertoire diversity, clonotypes and clonal homeostasis
#
# Paper    : Treg insufficiency in autoimmune diabetes is driven by selective
#            loss of Neuropilin-1 on intra-islet regulatory T cells
# Original : 2.Script/6.2.TCRTconv.R
# Input    : 1.Data/revision/TCR/filtered_contig_annotations.csv
#            1.Data/revision/scVITconvNoAPCs.rds
# Output   : 3.Output/revision/6.TCR/TconvTcr_*.png
# Upstream : 08_tconv_remove_apcs.R
# ----------------------------------------------------------------------------

library(Seurat)
library(scRepertoire)
library(ggpubr)
library(tidyverse)

# single cell inform readin ----
datTc <- readRDS('1.Data/revision/scVITconvNoAPCs.rds')

cells <- data.frame(genotype = datTc$type,
                    tissue = datTc$tissue,
                    celltype = datTc$predicted.id,
                    sampletype = 'Tconv')

# tcr data readin ----
tcrDf <- read.csv('1.Data/revision/TCR/filtered_contig_annotations.csv') %>% 
  mutate(barcode = gsub('-1', '', barcode)) %>% 
  merge(cells, by.x = 'barcode', by.y = 0)

table(table(tcrDf$barcode))

tcrDf$GenoCell <- paste(tcrDf$genotype, tcrDf$tissue, sep = "_")

contig_list <- split(tcrDf, f = tcrDf$GenoCell)

combined <- combineTCR(contig_list, 
                       samples = names(contig_list), 
                       ID = rep(c("L/+", 'WT'), each = 2))


clonalProportion(combined, 
                 cloneCall = "gene",
                 clonalSplit = c(5, 10, 50, 100, 500, 1000), exportTable = T) %>% 
  as.data.frame() %>% 
  rownames_to_column('sample') %>% 
  pivot_longer(names_to = 'ci',
               values_to = 'value',
               !sample) %>% 
  group_by(sample) %>% 
  separate(sample, into = c('genotype', 'tissue', 'no'), sep = '_', remove = F) %>% 
  mutate(ratio = value/sum(value),
         ci = factor(ci, levels = c('[1:5]', '[6:10]', '[11:50]', '[51:100]', '[101:500]', '[501:1000]')),
         variable=as.numeric(levels(ci)), 
         tissue = factor(tissue, levels = c('ndLN', 'islet')),
         genotype = factor(genotype, levels = c('WT', 'L/+'))) %>%
  ggplot(aes(x = genotype, y = ratio, fill = ci)) +
  geom_bar(stat = 'identity', col = 'gray20', linewidth = 1)+
  viridis::scale_fill_viridis(discrete = TRUE) +
  scale_x_discrete(labels=c("WT" = expression(bold("Cre control")), 
                            "L/+" = expression(bold(Rosa^{Nrp1})))) +
  facet_grid(~tissue) +
  labs(x = '',
       y = 'Occupied Repertoire Space',
       fill = 'Clonal Indices') +
  theme_bw() +
  theme(text = element_text(size = 15, face = 2, color = 'black'),
        axis.text = element_text(color = 'black', size = 15),
        strip.background = element_blank(),
        strip.text = element_text(size = 15),
        legend.position = 'bottom')


ggsave('3.Output/revision/6.TCR/TconvTcr_clonalProportion.png', width = 6, height = 6)

clonalHomeostasis(combined, cloneCall = "gene",
                  cloneSize = c(Rare = 0.005, 
                                Small = 0.01, 
                                Medium = 0.03, 
                                Large = 0.05, 
                                Hyperexpanded = 0.1), 
                  exportTable = T) %>% 
  as.data.frame() %>% 
  rownames_to_column('sample') %>% 
  pivot_longer(names_to = 'ci',
               values_to = 'value',
               !sample) %>% 
  group_by(sample) %>% 
  separate(sample, into = c('genotype', 'tissue', 'no'), sep = '_', remove = F) %>% 
  mutate(ratio = value/sum(value),
         ci = factor(ci, levels = c('Rare (0 < X <= 0.005)', 'Small (0.005 < X <= 0.01)', 'Medium (0.01 < X <= 0.03)', 'Large (0.03 < X <= 0.05)', 'Hyperexpanded (0.05 < X <= 0.1)')),
         genotype = factor(genotype, levels = c('WT', 'L/+'))) %>%
  ggplot(aes(x = genotype, y = ratio, fill = ci)) +
  geom_bar(stat = 'identity', col = 'gray20', linewidth = 1)+
  viridis::scale_fill_viridis(discrete = TRUE) +
  scale_x_discrete(labels=c("WT" = expression(bold("Cre control")), 
                            "L/+" = expression(bold(Rosa^{Nrp1})))) +
  facet_grid(~tissue) +
  labs(x = '',
       y = 'Occupied Repertoire Space',
       fill = 'Clonal Indices') +
  guides(fill = guide_legend(nrow = 3)) +
  theme_bw() +
  theme(text = element_text(size = 12, face = 2, color = 'black'),
        strip.background = element_blank(),
        legend.position = 'bottom')

ggsave('3.Output/revision/6.TCR/TconvTcr_clonalHomeostasis.png', width = 6.5, height = 6)

clonalQuant(combined, cloneCall="gene+nt", scale = T)

clonalQuant(combined, 
            cloneCall="gene", 
            group.by = "ID", 
            scale = TRUE)



clonalHomeostasis(combined, 
                  cloneCall = "gene",
                  cloneSize = c(Rare = 0.01, Small = 0.05, Medium = 0.1, Large = 0.3, Hyperexpanded =
                                  1))
clonalRarefaction(combined,
                  plot.type = 2,
                  hill.numbers = 1,
                  n.boots = 2)

clonalQuant(combined, cloneCall="gene+nt", scale = T, exportTable = T) %>% 
  separate(values, into = c('Geno', 'Cell', 'Geno2'), sep = '_') %>% 
  mutate(Cell = factor(Cell, levels = c("Naïve",
                                        'Memory',
                                        'Effector',
                                        'Proliferating',
                                        "Antigen \npresentation")),
         Geno = factor(Geno, levels = c('WT', 'L/+'))) %>% 
  ggplot(aes(Geno, scaled, fill = Cell)) +
  geom_bar(stat = 'identity') +
  labs(y = 'Percent of Unique Clonotype',
       x = 'Genotype') +
  facet_grid(~ Cell) +
  scale_fill_brewer(palette="Dark2") +
  coord_fixed(ratio = 1/25) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        text = element_text(size = 12, face = 2))
ggsave('3.Output/revision/6.TCR/TconvTcr_QuantifyClonotypes.png', width = 7, height = 5)


daDf <- clonalDiversity(combined, 
                cloneCall = "gene", 
                group.by = "sample", 
                x.axis = "ID", 
                n.boots = 100,
                exportTable = T)

daDf %>% as.data.frame() %>% 
  select(c(shannon, sample)) %>% 
  pivot_longer(!sample, names_to = "index", values_to = "count") %>% 
  separate(sample, into = c('Geno', 'Cell'), sep = '_') %>% 
  mutate(Cell = factor(Cell, levels = c("Naïve",
                                        'Memory',
                                        'Effector',
                                        'Proliferating',
                                        "Antigen \npresentation")),
         Geno = factor(Geno, levels = c('WT', 'L/+'))) %>% 
  ggplot(aes(Geno, count)) +
  geom_boxplot() +
  geom_point(aes(col = Cell), size = 5) +
  scale_color_brewer(palette="Dark2") +
  # facet_wrap(.~index, nrow = 1, scale = 'free_y') +
  theme_bw() +
  labs(x = '', y = 'Index Score') +
  theme(aspect.ratio=1,
        legend.position = 'bottom', 
        strip.background = element_blank(),
        text = element_text(size = 12, face = 2))

ggsave('3.Output/revision/6.TCR/TconvTcr_DiversityAnalysis_shannon.png', width = 6, height = 5)

daDf %>% as.data.frame() %>% 
  select(-ID) %>% 
  pivot_longer(!sample, names_to = "index", values_to = "count") %>% 
  separate(sample, into = c('Geno', 'Cell'), sep = '_') %>% 
  mutate(Cell = factor(Cell, levels = c("Naïve",
                                        'Memory',
                                        'Effector',
                                        'Proliferating',
                                        "Antigen \npresentation")),
         Geno = factor(Geno, levels = c('WT', 'L/+'))) %>% 
  ggplot(aes(Geno, count)) +
  geom_boxplot() +
  geom_point(aes(col = Cell), size = 5) +
  scale_color_brewer(palette="Dark2") +
  facet_wrap(.~index, nrow = 1, scale = 'free_y') +
  theme_bw() +
  labs(x = '', y = 'Index Score') +
  theme(aspect.ratio=1,
        legend.position = 'bottom', 
        strip.background = element_blank(),
        text = element_text(size = 12, face = 2))

ggsave('3.Output/revision/6.TCR/TconvTcr_DiversityAnalysis.png', width = 7, height = 5)
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