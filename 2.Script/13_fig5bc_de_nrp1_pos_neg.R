# ----------------------------------------------------------------------------
# Fig 5b-c | Differential expression between Nrp1-positive and Nrp1-negative Tregs
#
# Paper    : Treg insufficiency in autoimmune diabetes is driven by selective
#            loss of Neuropilin-1 on intra-islet regulatory T cells
# Original : 2.Script/5.2.DeGForproteinExp.R
# Input    : 1.Data/revision/datwithNrp1ProteinExpStat.RDS
# Output   : 3.Output/revision/5.DEanalysis/
# Upstream : 12_fig5a_citeseq_nrp1_split.R
# Note     : presto wilcoxauc
#            feeds the Venn diagram and the four-quadrant volcano plots
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(ggridges)
library(scales)
library(presto)
library(ggh4x)

dat <- readRDS("1.Data/revision/datNmf.rds")

protDf <- na.omit(data.frame(nrp1 = dat$antiNrp1, 
                             ctrl = dat$antiCtrl)) 

protDf$dis <- protDf$nrp1 - protDf$ctrl

protDf$clst <- ifelse(protDf$dis > summary(protDf$dis)["Median"], 
                      "Positive", "Negative")

subP <- protDf %>% filter(clst == "Positive") %>% 
  arrange(desc(dis))

subP <- subP[1:(nrow(subP) * 0.95),]

subN <- protDf %>% filter(clst == "Negative") %>% 
  arrange(-desc(dis)) 

subN <- subN[1: (nrow(subN)*0.95),]

protDfnew <- rbind(subP, subN)

dat <- AddMetaData(dat, 
                   protDfnew %>% 
                     select(clst) %>% 
                     rename(nrp1Prot = clst))

dat <- subset(dat, cells = rownames(protDfnew))

for (tp in c('L/+', 'WT')) {
  for (ts in c('islet', 'ndLN')){
    datSub <- subset(dat, 
                     tissue == ts & 
                       type == tp & 
                       cellType == "Treg")
    
    cat(dim(datSub), ' ')
    
    Idents(datSub) <- 'nrp1Prot'
    
    cat(ts, tp, '\n')
    
    wilcoxauc(datSub) %>% 
      filter(group == 'Positive') %>% 
      write.csv(paste0('3.Output/revision/5.DEanalysis/',
                       paste('Treg', 
                             ifelse(tp == 'L/+', 'Rosa', 'WT'),
                             ts,
                             'nrp1_pos_vs_nrp1_neg.csv', sep = '_')))
  }
}

cntMat <- GetAssayData(dat, assay = 'RNA', layer = 'counts')

for (tp in c('L/+', 'WT')) {
  for (ts in c('islet', 'ndLN')){
    for (nrp in c('Negative', 'Positive')) {
      datSub <- subset(dat, 
                       tissue == ts & 
                         type == tp & 
                         nrp1Prot == nrp &
                         cellType == "Treg")
      cat(tp, ts, nrp, length(Cells(datSub)),
      sum(cntMat[,Cells(datSub)] != 0)/length(Cells(datSub)), '\n')
    }
  }
}

sum(cntMat[, 1:2] != 0)/2

df <- rbind(
  read.csv("3.Output/revision/5.DEanalysis/Treg_Rosa_islet_nrp1_pos_vs_nrp1_neg.csv") %>% 
    filter(logFC != 0, pval < 0.05) %>% 
    add_column(sample = expression("bold(Rosa^Nrp1~islet)")), 
  read.csv("3.Output/revision/5.DEanalysis/Treg_WT_islet_nrp1_pos_vs_nrp1_neg.csv") %>% 
    filter(logFC != 0, pval < 0.05) %>% 
    add_column(sample = "bold(Cre~Control~islet)"),
  read.csv("3.Output/revision/5.DEanalysis/Treg_Rosa_ndLN_nrp1_pos_vs_nrp1_neg.csv") %>% 
    filter(logFC != 0, pval < 0.05) %>% 
    add_column(sample = "bold(Rosa^Nrp1~ndLN)"),
  read.csv("3.Output/revision/5.DEanalysis/Treg_WT_ndLN_nrp1_pos_vs_nrp1_neg.csv") %>% 
    filter(logFC != 0, pval < 0.05) %>% 
    add_column(sample = "bold(Cre~Control~ndLN)")) %>% 
  mutate(sample = factor(sample, 
                         levels = c(
                           "bold(Cre~Control~ndLN)",
                           "bold(Rosa^Nrp1~ndLN)",
                           "bold(Cre~Control~islet)",  
                           "bold(Rosa^Nrp1~islet)")))

featureUp <- strsplit("Ikzf2, Inpp4b, Tnfrsf4, Il2rb, Il2ra, Foxp3, Cd5, Zfp36l1, Ccr7, Junb, Jund, Jun, Fos, Ptprc, Nrp1, Ebi3", split = ", ")[[1]]


df <- df %>% 
  mutate(genotype = gsub('~ndLN|~islet', '', sample),
         tissue = factor(gsub('Rosa\\^Nrp1~|Cre~Control~', '', sample),
                         levels = c("bold(ndLN)", "bold(islet)" ))) 

df %>% 
  filter(! feature %in% c(featureUp, 
                          df %>% 
                            filter(logFC > 0) %>% 
                            arrange(-desc(pval)) %>% 
                            group_by(sample) %>%
                            group_modify(~ head(.x, 5)) %>% 
                            filter(! feature %in% featureUp) %>% 
                            pull(feature),
                          df %>% 
                            filter(logFC < 0) %>% 
                            arrange(-desc(pval)) %>% 
                            group_by(sample) %>%
                            group_modify(~ head(.x, 5)) %>% 
                            filter(! feature %in% featureUp) %>% 
                            pull(feature))) %>% 
  ggplot(aes(x = logFC, 
             y = -log10(pval))) +
  geom_point(col = "gray", alpha = 0.5) +
  # selected genes 
  geom_point(data = df %>% filter(feature %in% featureUp),
             aes(logFC, -log10(pval), col = "up", size = pct_in), alpha = 0.6) +
  ggrepel::geom_text_repel(data = df %>% filter(feature %in% featureUp),
                           aes(logFC, -log10(pval), label = feature, col = "up"),
                           fontface = 2, show.legend = F) +
  # up genes
  geom_point(data = df %>% 
               filter(logFC > 0) %>% 
               arrange(-desc(pval)) %>% 
               group_by(sample) %>%
               group_modify(~ head(.x, 5)) %>% 
               filter(! feature %in% featureUp),
             aes(logFC, -log10(pval), col = "up")) +
  ggrepel::geom_text_repel(data = df %>% 
                             filter(logFC > 0) %>% 
                             arrange(-desc(pval)) %>% 
                             group_by(sample) %>%
                             group_modify(~ head(.x, 5))%>% 
                             filter(! feature %in% featureUp),
                           
                           aes(logFC, -log10(pval), label = feature, col = "up"), 
                           show.legend = F, max.overlaps = Inf) +
  # down genes
  geom_point(data = df %>% 
               filter(logFC < 0) %>% 
               arrange(-desc(pval)) %>% 
               group_by(sample) %>%
               group_modify(~ head(.x, 5)) %>% 
               filter(! feature %in% featureUp),
             aes(logFC, -log10(pval), col = "dn", size = pct_in), alpha = 0.6) +
  ggrepel::geom_text_repel(data = df %>% 
                             filter(logFC < 0) %>% 
                             arrange(-desc(pval)) %>% 
                             group_by(sample) %>%
                             group_modify(~ head(.x, 5))%>% 
                             filter(! feature %in% featureUp),
                           aes(logFC, -log10(pval), label = feature, col = "dn"), 
                           show.legend = F) +
  scale_color_manual(name = "Status",
                     breaks = c("up", "dn"),
                     labels = c("Up-regualted genes", "Donw-regulated genes"),
                     values = c("red", "blue")) + 
  
  facet_grid2(tissue ~ genotype,  
             labeller=label_parsed, independent = "all", 
             scales = "free") +
  labs(y = bquote(bold(Significance~(-log^{italic(p)~value}))), 
       x = "Genes logFC value",
       size = 'Expressin \nratio (%)') +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_bw() +
  theme(aspect.ratio=2/3,
        legend.position = "bottom", 
        legend.box="vertical",
        legend.text = element_text(size = 12, face = 2, color = 'black'),
        text = element_text(size = 12, face = 2, color = 'black'),
        axis.text = element_text(size = 12, face = 2, color = 'black'), 
        strip.text = element_text(size = 15),
        strip.background = element_blank(),
        strip.placement = "inside")

ggsave("3.Output/revision/5.DEanalysis/Rosa_WT_islet_ndLN_genes_reNew1.png",
       width = 8, height = 7, units = "in", dpi = 300, bg = "white")
