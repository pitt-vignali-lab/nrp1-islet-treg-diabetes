# ----------------------------------------------------------------------------
# Analysis | Treg clonotype quantification
#
# Paper    : Treg insufficiency in autoimmune diabetes is driven by selective
#            loss of Neuropilin-1 on intra-islet regulatory T cells
# Original : 2.Script/6.2.1.TCRTreg.R
# Input    : 1.Data/revision/TCR/filtered_contig_annotations.csv
# Output   : 3.Output/revision/6.TCR/TregTcr_QuantifyClonotypes.png
# Upstream : 16_tcr_treg_diversity.R
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
                    rep = datTc$replicate,
                    sampletype = 'Treg')

# tcr data readin ----
tcrDf <- read.csv('1.Data/revision/TCR/filtered_contig_annotations.csv') %>% 
  mutate(barcode = gsub('-1', '', barcode)) %>% 
  merge(cells, by.x = 'barcode', by.y = 0)

table(table(tcrDf$barcode))

tcrDf$repMilo <- paste(tcrDf$rep, tcrDf$miloType, sep = "_")

contig_list <- split(tcrDf, f = tcrDf$repMilo)

contig_list <- contig_list[!grepl('NA', names(contig_list))]

sampleG <- gsub('NO.[0-9]+_', '', names(contig_list))
sampleG <- gsub('let', 'let_', sampleG)
sampleG <- gsub('LN', 'LN_', sampleG)

combined <- combineTCR(contig_list, 
                       samples = names(contig_list), 
                       ID = sampleG)

# clonalQuant(combined, cloneCall="gene+nt", scale = T)

my_comparisons <- list( c('WT', 'L/+'))

clonalQuant(combined, cloneCall="gene", exportTable = T) %>% 
  separate(values, into = c('rep', 'condi', 'Cell', 'Geno'), sep = '_') %>% 
  mutate(condi = factor(condi, levels = c('ndLNWT', 'ndLNL/+', 'isletWT', 'isletL/+')),
         sample = gsub('WT|L/\\+', '', Geno),
         Cell = factor(Cell, levels = c('ndLN', 'islet')),
         Geno = factor(Geno, levels = c('WT', 'L/+'))) %>%
  ggplot(aes(Geno, total, col = Geno)) +
  geom_boxplot() +
  geom_jitter() +
  facet_grid(~Cell) +
  stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test', 
                     paired = T, label = "p.signif") +
  # geom_bar(stat = 'identity',  alpha = 0.7) +
  labs(y = 'Total number of Unique Clonotype',
       x = 'Genotype') +
  # facet_grid(~ Cell) +
  scale_color_manual(breaks = c('L/+', 'WT'),
                     values = alpha(c('red', 'blue'), 0.8)) +
  coord_fixed(ratio = 1/25) +
  scale_x_discrete(labels=c("WT" = expression(bold("Cre control")), 
                            "L/+" = expression(bold(Rosa^{Nrp1})))) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        text = element_text(size = 15, face = 2, color = 'black'),
        axis.text = element_text(size = 12, face = 2, color = 'black'))

ggsave('3.Output/revision/6.TCR/TregTcr_QuantifyClonotypes.png', width = 5, height = 5)


