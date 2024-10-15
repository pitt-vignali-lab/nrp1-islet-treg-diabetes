# ----------------------------------------------------------------------------
# Fig 4 | UMAPs, cluster composition statistics and pseudobulk aggregation
#
# Paper    : Treg insufficiency in autoimmune diabetes is driven by selective
#            loss of Neuropilin-1 on intra-islet regulatory T cells
# Original : 2.Script/3.4.scRNAstatistical.R
# Input    : 1.Data/revision/scVITreg.rds
#            1.Data/revision/scVITconvNoAPCs.rds
# Output   : 3.Output/revision/3.ScVI/TconvUMAPnoAPCs.png
#            3.Output/revision/3.ScVI/pseudoBuld/csv/
# Upstream : 08_tconv_remove_apcs.R
# Note     : Writes the per-genotype per-tissue pseudobulk tables used by Fig 4 and Supp Table 1
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(ggpubr)
library(COSG)
library(export)

# T conv -----
datTc <- readRDS('1.Data/revision/scVITconvNoAPCs.rds')

data.frame(Embeddings(datTc, 'umap'), id =  datTc$predicted.id) %>% 
  as_tibble() %>% 
  mutate(id = factor(id, levels = c('Memory',
                                    'Activated',
                                    'Effector',
                                    'Proliferating'))) %>% 
  ggplot(aes(umap_1, umap_2, col = id)) +
  geom_point(size = 2, alpha = 0.8) +
  coord_fixed() +
  scale_color_brewer(palette="Dark2") +
  theme_bw() +
  guides(
    color = guide_legend(override.aes = list(size=5),
                         nrow = 2))+
  labs(col = 'Cell type') +
  theme(legend.position = 'bottom',
        text = element_text(size = 15, face = 2),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

ggsave('3.Output/revision/3.ScVI/TconvUMAPnoAPCs.png', width = 5, height = 5, bg = 'white')

table(datTc$predicted.id, datTc$tissue, datTc$type, datTc$replicate) %>% 
  as.data.frame() %>% 
  mutate(Var1 = factor(Var1, levels = c('Memory',
                                        'Activated',
                                        'Effector',
                                        'Proliferating')),
         Var2 = factor(Var2, levels = c('ndLN', 'islet')),
         Var3 = factor(Var3, levels = c('WT', 'L/+'))) %>% 
  group_by(Var4) %>% 
  mutate(ratio = Freq/sum(Freq)) %>% 
  filter(Freq != 0) %>% 
  ggplot(aes(x = Var3, y = ratio, col = Var1)) +
  geom_boxplot() +
  geom_point() +
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                         comparisons = list(c('L/+', 'WT')), 
                     test =  'wilcox.test',
                     face = 2) +
  
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.52)) +
  facet_grid(Var2 ~ Var1) +
  scale_color_brewer(palette="Dark2") +
  scale_x_discrete(labels=c("WT" = expression(bold("Cre control")), 
                            "L/+" = expression(bold(Rosa^{Nrp1})))) +
  # coord_fixed(ratio = 1/15) +
  labs(x = '',
       y = 'Percentage') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 12, face = 2, color = 'black'),
        axis.text = element_text(size = 10, face = 2, color = 'black'),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = 2, color = 'black'))

ggsave('3.Output/revision/3.ScVI/TconvBarplotnoAPCswithStat.png', width = 7.5, height = 5, bg = 'white')

# table(datTc$predicted.id, datTc$tissue, datTc$type, datTc$replicate) %>% 
#   as.data.frame() %>% 
#   mutate(Var1 = factor(Var1, levels = c('Memory',
#                                         'Activated',
#                                         'Effector',
#                                         'Proliferating'))) %>% 
#   group_by(Var2, Var3) %>% 
#   mutate(ratio = Freq/sum(Freq)) %>% 
#   ggplot(aes(ratio, Var3, fill = Var1)) +
#   geom_bar(stat = 'identity',position=position_fill(reverse=T)) +
#   facet_grid(Var2~.) +
#   coord_fixed(1/8)+
#   scale_fill_brewer(palette="Dark2") +
#   theme_bw() +
#   labs(fill = 'Cell type',
#        y = ' ') +
#   guides(
#     fill = guide_legend(nrow = 2,
#                         keyheight = 1.2,
#                         keywidth = 1.2))+
#   scale_y_discrete(labels=c("WT" = expression(bold("Cre control")), 
#                             "L/+" = expression(bold(Rosa^{Nrp1})))) +
#   labs(col = 'Cell type') +
#   theme(legend.position = 'bottom',
#         text = element_text(size = 12, face = 2),
#         strip.background = element_blank())
# 
# ggsave('3.Output/revision/3.ScVI/TconvBarplotnoAPCs.png', width = 6, height = 4, bg = 'white')

Idents(datTc) <- 'predicted.id'
# naive Ccr7 ++
rbind(DotPlot(datTc, 
        features = c('Ccr7', 'Tcf7', 'Sell'))$data %>% 
          mutate(celltype = 'Memory'), 
  # Effector
  DotPlot(datTc, 
          features = c('Ifng', 'Cxcr3', 'Pdcd1'))$data %>% 
    mutate(celltype = 'Effector'),
  # Proliferating
  DotPlot(datTc, 
          features = c('Top2a', 'Mki67'))$data %>% 
    mutate(celltype = 'Proliferating'), 
  # Activated
  DotPlot(datTc, 
          features = c('S100a9', 'S100a8'))$data %>% 
          mutate(celltype = 'Activated',
                 pct.exp = pct.exp * 10)) %>% 
  mutate(id = factor(id, levels = c('Memory',
                                    'Activated',
                                    'Effector',
                                    'Proliferating')),
         celltype = factor(celltype, levels = c('Memory',
                                    'Activated',
                                    'Effector',
                                    'Proliferating'))) %>% 
  ggplot(aes(features.plot, id)) +
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21) +
  scale_fill_gradient2(high = 'red', low = 'blue', mid = 'white', midpoint = 0) +
  facet_grid(~celltype, scales = "free", space = 'free') +
  theme_bw() +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, 
                                theme = theme(
                                  legend.key.width  = unit(10, "lines")
                                )),
         size = guide_legend(title.position="top", title.hjust = 0.5)) +
  labs(x = 'Marker gene',
       y = '',
       size = 'Expression percentage',
       fill = 'Scaled average expresion') +
  theme(legend.position = 'bottom',
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 12, face = 2, color = 'black'),
        text = element_text(size = 12, face = 2),
        strip.background = element_blank(),
        strip.text = element_text(color = 'black', size = 15))

ggsave('3.Output/revision/3.ScVI/TconvMkG.png', width = 7.5, height = 4, bg = 'white')

# Treg ----
datTr <- readRDS('1.Data/revision/scVITreg.rds')

datTr@meta.data <- datTr@meta.data %>% 
  mutate(predicted.id = ifelse(predicted.id == "Naïve", 'Memory',
                               ifelse(predicted.id == 'Memory', 'Activated', 'Effector')))

Idents(datTr) <- 'predicted.id'
data.frame(Embeddings(datTr, 'umap'), id =  datTr$predicted.id) %>% 
  as_tibble() %>% 
  mutate(id = factor(id, levels = c('Memory',
                                    'Activated',
                                    'Effector'))) %>% 
  ggplot(aes(UMAP_1, UMAP_2, col = id)) +
  geom_point(size = 2, alpha = 0.8) +
  coord_fixed() +
  scale_color_brewer(palette="Dark2") +
  theme_bw() +
  guides(
    color = guide_legend(override.aes = list(size=5),
                         nrow = 1))+
  labs(col = 'Cell type') +
  theme(legend.position = 'bottom',
        text = element_text(size = 20, face = 2),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

ggsave('3.Output/revision/3.ScVI/TregUMAP.png', width = 7, height = 5, bg = 'white')


table(datTr$predicted.id, datTr$tissue, datTr$type) %>% 
  as.data.frame() %>% 
  mutate(Var1 = factor(Var1, levels = c('Memory',
                                        'Activated',
                                        'Effector')),
         Var2 = factor(Var2, levels = c('ndLN', 'islet'))) %>% 
  group_by(Var2, Var3) %>% 
  mutate(ratio = (100*Freq)/sum(Freq)) %>% 
  ggplot(aes(ratio, Var3, fill = Var1)) +
  scale_x_continuous(labels = scales::percent) +
  geom_bar(stat = 'identity',position=position_fill(reverse=T)) +
  scale_y_discrete(labels=c("WT" = expression(bold("Cre control")), 
                            "L/+" = expression(bold(Rosa^{Nrp1})))) +
  facet_grid(Var2~.) +
  coord_fixed(1/8)+
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  labs(fill = 'Cell type',
       y = '',
       x = 'Percentage') +
  guides(
    fill = guide_legend(nrow = 1,
                        keyheight = 1.2,
                        keywidth = 1.2))+
  labs(col = 'Cell type') +
  theme(legend.position = 'bottom',
        text = element_text(size = 15, face = 2),
        axis.text = element_text(color = 'black', size = 13),
        strip.text = element_text(size = 15),
        strip.background = element_blank())

ggsave('3.Output/revision/3.ScVI/TregBarplot.png', width = 6, height = 4, bg = 'white')

# naive Ccr7 ++
rbind(DotPlot(datTr, 
              features = c('Ccr7', 'Tcf7', 'Sell'),
              group.by = )$data %>% 
        mutate(celltype = 'Memory'), 
      # Effector
      DotPlot(datTr, 
              features = c('Ifng', 'Cxcr3'))$data %>% 
        mutate(celltype = 'Effector'),
      # Memory Ccr7 +
      DotPlot(datTr,
              features = c('S100a8', 'S100a9'))$data %>%
        mutate(celltype = 'Activated',
               pct.exp = pct.exp * 150)) %>% 
  mutate(id = factor(id, levels = c('Memory',
                                    'Activated',
                                    'Effector',
                                    'Proliferating')),
         celltype = factor(celltype, levels = c('Memory',
                                                'Activated',
                                                'Effector',
                                                'Proliferating'))) %>% 
  ggplot(aes(features.plot, id)) +
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21) +
  scale_fill_gradient2(high = 'red', low = 'blue', mid = 'white', midpoint = 0) +
  facet_grid(~celltype, scales = "free", space = 'free') +
  theme_bw() +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, 
                                theme = theme(
                                  legend.key.width  = unit(10, "lines")
                                )),
         size = guide_legend(title.position="top", title.hjust = 0.5)) +
  labs(x = 'Marker gene',
       y = '',
       size = 'Expression percentage',
       fill = 'Scaled average expresion') +
  theme(legend.position = 'bottom',
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 13, color = 'black'),
        text = element_text(size = 12, face = 2),
        strip.background = element_blank(),
        strip.text = element_text(size = 15))

ggsave('3.Output/revision/3.ScVI/TregMkG.png', width = 6, height = 4, bg = 'white')


if(!dir.exists('3.Output/revision/3.ScVI/pseudoBuld')){
  dir.create('3.Output/revision/3.ScVI/pseudoBuld')
}

datTr$condi <- paste(datTr$type, datTr$tissue, datTr$replicate, sep = '_')

unique(datTr$condi)

Idents(datTr) <- 'condi'

deG <- presto::wilcoxauc(datTr) %>% 
  mutate(group = gsub('L\\/\\+', 'RosaNrp1', group),
         group = gsub('WT', 'CreControl', group),
         group = gsub('NO.', 'Replicate.', group))

for (i in unique(deG$group)) {
  write.csv(deG %>% 
              filter(group == i) %>% 
              select(-group), 
            paste0('3.Output/revision/3.ScVI/pseudoBuld/csv/', i ,'.csv'),
            row.names = F, quote = F)
}


deGTrans <- deG %>% 
  pivot_wider(id_cols = feature, names_from = group, values_from = logFC) %>% 
  column_to_rownames('feature')

pca_res <- prcomp(deGTrans)

pca_res[["rotation"]] %>% 
  as.data.frame() %>% 
  rownames_to_column('sample') %>% 
  mutate(sample1 = factor(gsub('_NO.[0-9]+', '', sample))) %>% 
  ggplot(aes(PC2, PC1, fill = sample1)) +
  geom_point(shape = 21, size = 5) +
  coord_fixed() +
  theme_bw() +
  theme(# aspect.ratio = 1,
        legend.position = 'none',
        text = element_text(size = 15, face = 2, color = 'black'),
        axis.text = element_text(size = 15, face = 2, color = 'black'))
ggsave('3.Output/revision/3.ScVI/TregSamplePCA.png', width = 6, height = 5)
