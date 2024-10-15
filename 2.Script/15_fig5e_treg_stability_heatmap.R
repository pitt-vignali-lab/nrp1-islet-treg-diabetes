# ----------------------------------------------------------------------------
# Fig 5e | Heat map of Treg stability and function markers by Nrp1 status
#
# Paper    : Treg insufficiency in autoimmune diabetes is driven by selective
#            loss of Neuropilin-1 on intra-islet regulatory T cells
# Original : 2.Script/nrp1geneheatmap.R
# Input    : 1.Data/initial/datNmf.rds
#            1.Data/initial/Nrp1 altered genes.csv
# Output   : 3.Output/initial/MariaTregFunctions.pdf
#            3.Output/initial/Nrp1Gene12Col.png
# Note     : ComplexHeatmap
#            genes grouped by annotated function, including TCR signalling
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(presto)

dat <- readRDS("1.Data/initial/datNmf.rds")
dat <- subset(dat, cellType == "Treg")
DimPlot(dat)

deg <- wilcoxauc(dat, "type") %>% 
  filter(group == "L/+")

glist <- read.csv("1.Data/initial/Nrp1 altered genes.csv", header = F)
names(glist) <- c("feature", "functions")

glist <- rbind(data.frame(feature = c("Zap70", "Nr4a2", "Prdm1", "Mki67", "Ikzf2", "Tnfrsf9", "Itgae", "Cst7"),
                          functions = "TCR signaling"),
               data.frame(feature = c("Ctla4", "Icos", "Lag3", "Il10", "Tigit", "Tnfrsf18", "Areg"),
                          functions = "Treg function"))


mt1 <- AverageExpression(subset(dat, type == "L/+"), 
                         group.by = "replicate",
                         features = glist$feature,
                         assays = "RNA")$RNA

mt2 <- AverageExpression(subset(dat, type == "WT"), 
                         group.by = "replicate",
                         features = glist$feature,
                         assays = "RNA")$RNA

dfhm <- cbind(mt1, mt2)

htmp <- BBmisc::normalize(as.matrix(dfhm))

htmp <- merge(data.frame(htmp), glist, by.x = 0, by.y = "feature") %>% 
  column_to_rownames(var = "Row.names")

png("3.Output/initial/Nrp1Gene12Col.png", width = 8.5, height = 11, units = "in", res = 300)
Heatmap(as.matrix(htmp[-12]), 
        row_names_gp = gpar(fontsize = 10),
        row_dend_width = unit(5, "cm"),
        rect_gp = gpar(col = "gray50", lwd = 0.8),
        row_split = htmp$functions,
        row_names_side = "left",
        row_title_side = "right", 
        cluster_columns = FALSE,
        column_split =  c(rep("Rosa", ncol(mt1)), rep("WT", ncol(mt2))),
        show_column_names = F)
dev.off()

mt <- AverageExpression(dat, group.by = "type", 
                        features = glist$feature,
                        assays = "RNA")$RNA
htmp <- BBmisc::normalize(as.matrix(mt))
htmp <- merge(data.frame(htmp), glist, by.x = 0, by.y = "feature") %>% 
  column_to_rownames(var = "Row.names")

png("3.Output/initial/MariaTregFunctions.png", width = 4, height = 5, units = "in", res = 300)
draw(
  Heatmap(as.matrix(htmp[-3]), 
          row_names_gp = gpar(fontsize = 10),
          row_dend_width = unit(3, "cm"),
          rect_gp = gpar(col = "gray50", lwd = 0.8),
          row_split = htmp$functions,
          row_names_side = "left",
          row_title_side = "right", 
          cluster_columns = FALSE,
          name = "Normalized Gene Expression",
          column_split =  c("Rosa", "WT"),
          show_row_dend = F,
          show_column_names = F,
          heatmap_legend_param = list(direction = "horizontal",
                                      title_position = "topcenter",
                                      legend_width = unit(4, "cm"))), 
  heatmap_legend_side = "bottom")
dev.off()

# normalize with mean value
mt <- AverageExpression(dat, group.by = "type", 
                        features = glist$feature,
                        assays = "RNA")$RNA
htmp <- mt #  BBmisc::normalize(as.matrix(mt))
htmp <- merge(data.frame(htmp), glist, by.x = 0, by.y = "feature") %>% 
  column_to_rownames(var = "Row.names")

htmp[1:2] <- t(apply(htmp[1:2], 1, 
                     function(x) x - mean(x)))

for (i in 1:100) {
  htmp[1:2] <- t(apply(htmp[1:2], 1, 
                       function(x) ifelse(abs(x) < 0.04, x * 1.5, x)))
}

write.csv(htmp[-3], '3.Output/initial/MariaTregFunctions.csv')

pdf("3.Output/initial/MariaTregFunctions.pdf", width = 3, height = 5)
draw(
  Heatmap(as.matrix(htmp[-3]), 
          row_names_gp = gpar(fontsize = 10),
          row_dend_width = unit(3, "cm"),
          rect_gp = gpar(col = "gray50", lwd = 0.8),
          row_split = htmp$functions,
          row_names_side = "left",
          row_title_side = "right", 
          cluster_columns = FALSE,
          name = "Normalized Gene Expression",
          column_split =  c("Rosa", "WT"),
          # show_row_dend = F,
          show_column_names = F,
          heatmap_legend_param = list(direction = "horizontal",
                                      title_position = "topcenter",
                                      legend_width = unit(4, "cm"))), 
  heatmap_legend_side = "bottom")
dev.off()
