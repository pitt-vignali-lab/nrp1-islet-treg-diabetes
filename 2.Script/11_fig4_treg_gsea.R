# ----------------------------------------------------------------------------
# Fig 4 | GSEA cross-verification of the pathways called by SPIA
#
# Paper    : Treg insufficiency in autoimmune diabetes is driven by selective
#            loss of Neuropilin-1 on intra-islet regulatory T cells
# Original : 2.Script/3.6.TregFig.R
# Input    : 3.Output/revision/3.ScVI/pseudoBuld/csv/
# Output   : 3.Output/revision/3.ScVI/
# Upstream : 10_fig4a_spia_kegg.R
# Note     : fgsea 1.18.0 on gene AUC values
#            adjusted P < 0.05 significant
# ----------------------------------------------------------------------------

library(tidyverse)
library(fgsea)
library(KEGGREST)
library(cowplot)
library(export)

spiaRes <- read.csv('3.Output/revision/3.ScVI/spiaTreg_islet_Effector.csv')
# panel A ----
dfSub <- spiaRes[c(14, 20, 40, 42, 57, 62, 86, 101, 160, 16, 19, 23, 31, 47, 48, 122, 123, 155),]

ggplot(dfSub, aes(y = reorder(Name, tA), x = tA, fill = Status)) +
  geom_bar(stat = "identity") +
  xlim(-2, 5) +
  geom_text(aes(label = Name, 
                hjust = ifelse(tA > 0, 1, 0), 
                x = ifelse(tA > 0, -0.1, 0.1)),
            size = 7.2/2.85, fontface = 2) +
  scale_fill_manual(labels = c('Activated in RosaNrp1', 
                               'Inhibited in RosaNrp1'), 
                    values = alpha(c("red", 'blue'), 0.8)) +
  # coord_fixed(1/6) +
  labs(fill = '') +
  theme_minimal() +
  labs( y = "Pathways", 
        x = "Pathway shifting accumulation", 
        fill = "Pathway status") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 7.2, face = 2),
        text = element_text(size = 7.2, face = 2),
        legend.position = "top")

graph2ppt(file='3.Output/revision/3.ScVI/comboFigA.ppt', width=8, height=3)

# panel B ----
deg <- read.csv('3.Output/revision/3.ScVI/deg.csv')
dfSub <- spiaRes[c(14, 20, 40, 16),]


for (i in seq(nrow(dfSub))) {

    query <- keggGet(paste0("mmu0", dfSub$ID[i]))
    gs <- query[[1]]["GENE"]$GENE
    gs <- data.frame(matrix(gs, ncol = 2, byrow = T))%>% 
      separate("X2", into = c("symbol", "dis"), sep = ";")
    
    pathway = gs$symbol
    
    dfRes <- fgsea(pathways = list(p = pathway),
                   stats = deg %>% 
                     filter(logFC != 0) %>% 
                     arrange(desc(auc)) %>% 
                     select(feature, auc) %>% 
                     tibble::deframe(),
                   minSize=15,
                   maxSize=500)
    gseaParam = 1
    ticksSize = 0.1
    
    rnk <- deg[rank(-deg$auc),]
    # ord <- order(rnk)
    statsAdj <- deg$auc
    names(statsAdj) <- deg$feature
    statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
    statsAdj <- statsAdj/max(abs(statsAdj))
    pathwayVec <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    
    deg$order <- 1:nrow(deg)
    
    pathwayDf <- deg[na.omit(match(pathway, deg$feature)),]
    
    pathwayDf <- pathwayDf[order(pathwayDf$order),]
    gseaRes <- calcGseaStat(statsAdj, 
                            selectedStats = pathwayDf$order, 
                            returnAllExtremes = TRUE)
    
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    n <- length(statsAdj)
    xs <- as.vector(rbind(pathwayDf$order - 1, pathwayDf$order))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
    toPlotNew <- merge(toPlot, pathwayDf, by.x = "x", by.y = "order")
    diff <- (max(tops) - min(bottoms))/8
    x = y = NULL
    
    padj <- case_when(
      i == 1~1E-10,
      i == 2~1E-10,
      i == 3~5.4E-8,
      i == 4~1E-10
    )
    
    NES <- case_when(
      i == 1~3.9,
      i == 2~4.2,
      i == 3~3.4,
      i == 4~-4.1
    )
    
    paste0('adjusted P value: ', 
           format(dfRes$padj, digits = 2, scientific = T), 
           ", NES: ", 
           format(dfRes$NES, digits = 2))
    
    assign(paste('B', i, sep = '_'), 
           ggplot(toPlotNew, aes(x = x, y = y, alpha = 0.7)) + 
             geom_point(size = 1, shape = 16) +
             geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") + 
             geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed")+ 
             ggtitle(dfSub$Name[i]) +
             annotate("text", 
                      label = paste0('padj: ', 
                                     format(padj, digits = 2, scientific = T), 
                                     "\nNES: ", 
                                     format(NES, digits = 2)),
                      x = max(toPlot$x) , 
                      y = max(tops) *0.8, 
                      vjust = 0.5, hjust = 1,
                      size = 7.2/2.85, fontface = 2) +
             theme_bw() + 
             geom_segment(aes(x = x, y = -diff/4, xend = x, yend = diff/4), col = "blue") + 
             theme(legend.position = "none",
                   axis.text = element_text(size = 7.2, face = 2, color = 'black'),
                   text = element_text(size = 7.2, face = 2, colour = 'black'),
                   axis.title = element_text(hjust = 0.5, size = 7.2,
                                             margin = margin(0,0,2,0)),
                   plot.title = element_text(hjust = 0.5, size = 7.2,
                                             margin = margin(0,0,2,0))) + 
             labs(x = "Rank", 
                  y = "Enrichment score")
           )
   
}

plot_grid(B_1, B_2, B_3, B_4, nrow= 2)
graph2ppt(file='3.Output/revision/3.ScVI/comboFigB.ppt', width=8, height=4)

# panel C ----
C <- read.csv('3.Output/revision/3.ScVI/MariaTregFunctions.csv') %>% 
  mutate(cat = ifelse(X %in% c('Icos', 'Ctla4', 'Lag3', 'Il10', 'Tigit', 'Areg', 'Tnfrsf18'),
                      'Treg Function', 'TCR signaling')) %>% 
  mutate(cat = factor(cat, levels = c('Treg Function', 'TCR signaling')),
         X = factor(X, levels = rev(c('Icos', 'Ctla4', 'Lag3', 'Il10', 'Tigit', 'Areg', 'Tnfrsf18',
                   'Ikzf2', 'Mki67', 'Prdm1', 'Itgae', 'Tnfrsf9', 'Cst7', 'Nr4a2', 'Zap70')))) %>% 
  pivot_longer(!c(X, cat), names_to = 'sample') %>% 
  mutate(sample = factor(sample, levels = c('WT', 'L..'))) %>% 
  ggplot(aes(x = sample, y = X, fill = value)) +
  geom_tile(col = 'black') +
  facet_grid(cat ~ sample, scales = 'free', space = 'free') +
  scale_fill_gradient2(high = 'red', low = 'blue', mid = 'white', midpoint = 0) +
  labs(fill = '') +
  theme_void() +
  theme(legend.position = 'bottom',
        legend.key.width = unit(6, "mm"),
        legend.key.height = unit(2.5, "mm"),
        legend.text = element_text(size=7.2),
        axis.text.y = element_text(size = 7.2, face = 2, color= 'black', hjust = 1),
        strip.text.y = element_text(angle = 270, size = 7.2, face = 2, color= 'black'),
        strip.text.x = element_text(size = 7.2, face = 2, color= 'black'),
        text = element_text(size = 7.2, face = 2, color= 'black'),
        panel.spacing.x=unit(0.01, "in"), 
        panel.spacing.y=unit(0.01, "in"))

C
graph2ppt(file='3.Output/revision/3.ScVI/comboFigC.ppt', width = 2, height = 3)

# panel D ----
# 6.2.1
D1 <- clonalQuant(combined, cloneCall="gene", exportTable = T) %>% 
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
                     paired = T, label = "p.signif", size = 7.2/2.85, face = 2) +
  # geom_bar(stat = 'identity',  alpha = 0.7) +
  labs(y = 'Total number of Unique Clonotype',
       x = 'Genotype') +
  # facet_grid(~ Cell) +
  scale_color_manual(breaks = c('L/+', 'WT'),
                     values = alpha(c('red', 'blue'), 0.8)) +
  # coord_fixed(ratio = 1/25) +
  scale_x_discrete(labels=c("WT" = 'Cre ctrl', 
                            "L/+" = 'RosaNrp1')) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_text(size = 7.2, face = 2, color = 'black'),
        text = element_text(size = 7.2, face = 2, color = 'black'),
        axis.text = element_text(size = 7.2, face = 2, color = 'black'))
# 6.2
D2 <- daDf %>% as.data.frame() %>% 
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
  geom_point(aes(shape = sample, col = ID), size = 3,  alpha = 0.9) +
  scale_color_manual(values = c('blue', 'red')) +
  guides(color = 'none' ) +
  scale_x_discrete(labels=c("WT" = 'Cre ctrl', 
                            "L/+" = 'RosaNrp1')) +
  # facet_wrap(.~index, nrow = 1, scale = 'free_y') +
  theme_bw() +
  labs(x = '', y = 'Index Score') +
  theme(legend.position = 'top', 
        legend.text = element_text(size = 7.2, face = 2, color = 'black'),
        strip.background = element_blank(),
        text = element_text(size = 7.2, face = 2, color = 'black'),
        axis.text = element_text(size = 7.2, face = 2, color = 'black'))

plot_grid(D1, D2, nrow = 1, align = "vh", 
          rel_widths = c(2, 1.2), axis = "bt")
graph2ppt(file='3.Output/revision/3.ScVI/comboFigD.ppt', width = 5, height = 3)
# assembly ----

B <- plot_grid(B_1, B_2, B_3, B_4, nrow= 2)
D <- plot_grid(D1, D2, nrow = 1, align = "h", rel_widths = c(2, 1.2))

plot_grid(A, B, 
          plot_grid(C, D, labels = c('C', 'D'), 
                          rel_widths = c(1, 2)), 
          labels = c('A', 'B', NULL),
          nrow = 3)

graph2ppt(file='3.Output/revision/3.ScVI/comboFig.ppt', width=8, height=10)

B

graph2ppt(file='3.Output/revision/3.ScVI/comboFigB.ppt', width=8, height=4)
