# ----------------------------------------------------------------------------
# Fig 5a | Split Tregs into Nrp1-positive and Nrp1-negative using CITE-seq signal
#
# Paper    : Treg insufficiency in autoimmune diabetes is driven by selective
#            loss of Neuropilin-1 on intra-islet regulatory T cells
# Original : 2.Script/5.1.Nrp1CiteSeqTregSplit.R
# Input    : 1.Data/revision/scVITreg.rds
#            1.Data/revision/Nrp1umi_count/
# Output   : 1.Data/revision/datwithNrp1ProteinExpStat.RDS
#            3.Output/revision/5.Nrp1ProteinExpressionTest/
# Upstream : 05_scvi_integration.R
# Note     : TotalSeq-C0552 anti-Nrp1 with matched IgG2a isotype control
# ----------------------------------------------------------------------------

if(! dir.exists('3.Output/revision/5.DEanalysis')) {
  dir.create('3.Output/revision/5.DEanalysis')
}



library(Seurat)
library(tidyverse)
library(ggridges)
library(scales)
library(presto)

dat <- readRDS("1.Data/revision/datNmf.rds")
# Treg Tconv comparison ----
protDf <- na.omit(rbind(data.frame(ab = "AntiNrp1",
                                   abVal = dat$antiNrp1,
                                   cellType = dat$cellType),
                        data.frame(ab = "Control",
                                   abVal = dat$antiCtrl,
                                   cellType = dat$cellType)))

median(protDf[protDf$ab == "AntiNrp1" & protDf$cellType == "Treg",]$abVal) - 
  median(protDf[protDf$ab == "Control" & protDf$cellType == "Treg",]$abVal)

median(protDf[protDf$ab == "AntiNrp1" & protDf$cellType == "Tconv",]$abVal) - 
  median(protDf[protDf$ab == "Control" & protDf$cellType == "Tconv",]$abVal)

medianDf <- protDf %>% group_by(ab, cellType) %>% summarise(medians = median(abVal))

ggplot(protDf) +
  geom_density(alpha = 0.2, aes(x = abVal, fill = ab, col = ab))+
  geom_vline(data = medianDf, aes(xintercept = medians, col = ab),
             linetype="dashed") +
  geom_text(data = medianDf, aes(x = medians, label = medians, y = 1.5)) +
  scale_x_log10("x",
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand = c(0, 0)) +
  ylab("Density") +
  facet_wrap(~ cellType, ncol = 1) +
  theme_bw() +
  coord_fixed()+
  theme(axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # legend.position = c(0.8, 0.95),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        strip.background = element_blank(),
        text = element_text(size = 12, face = 2))

ggsave("3.Output/revision/5.Nrp1ProteinExpressionTest/Nrp1_Protein_expression_comparison_between_Treg_Tconv.png",
       width = 6, height = 6, units = "in", dpi = 300, bg = "white")

# method 1 median ----
# dat <- subset(dat, cellType == "Treg")

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

ggplot(protDfnew, aes(x = dis, y = -0.2)) +
  
  # horizontal boxplots & density plots
  geom_boxplot(aes(fill = clst, col = clst), 
               alpha = 0.3, width = 0.3, show.legend = F) +
  geom_density(aes(x = dis, fill = clst, col = clst), 
               inherit.aes = FALSE, alpha = 0.3) +
  scale_x_log10("x",
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand = c(0, 0)) +
  geom_vline(lty = 4, color = "black",
             aes(xintercept = median(protDf$dis), size = "Threshold \n(median value)")) +
  scale_fill_manual(name = "Categry",
                    breaks = c("Positive", "Negative"),
                    values = c("red", "blue")) +
  scale_color_manual(name = "Categry",
                     breaks = c("Positive", "Negative"),
                     values = c("red", "blue"))+
  scale_size_manual(name = "", values = 1) +
  labs(x = "", y = "Density") +
  theme_bw() +
  coord_fixed()+
  theme(axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # legend.position = c(0.8, 0.95),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        strip.background = element_blank(),
        text = element_text(size = 10, face = 2))
ggsave("3.Output/revision/5.Nrp1ProteinExpressionTest/Treg_cells_Protein_expression_threshold_settings.png",
       width = 5, height = 3, units = "in", dpi = 300, bg = "white")

dat <- AddMetaData(dat, protDfnew %>% select(clst) %>% rename(nrp1Prot = clst))

# method 2 Kmean ----
# t <- kmeans(protDf$dis, 2)
# 
# protDf$clst <- t[['cluster']]
# 
# if (mean(protDf[protDf$clst == 1,]$dis) > mean(protDf[protDf$clst == 2,]$dis)){
#   protDf[protDf$clst == 1,]$clst = "Positive"
#   protDf[protDf$clst == 2,]$clst = "Negative"
# } else{
#   protDf[protDf$clst == 2,]$clst = "Positive"
#   protDf[protDf$clst == 1,]$clst = "Negative"
# }


# ----
table(protDfnew$clst)

protDfnew %>% filter(clst == "Positive") %>% summarise(median(nrp1) - median(ctrl))
protDfnew %>% filter(clst == "Negative") %>% summarise(median(nrp1) - median(ctrl))

rbind(protDfnew %>% select(c(nrp1, clst)) %>% add_column(Antibody = "Anti-Nrp1") %>% rename(prot = nrp1),
      protDfnew %>% select(c(ctrl, clst)) %>% add_column(Antibody = "Control") %>% rename(prot = ctrl))  %>% 
  ggplot(aes(x = prot, fill = Antibody, col = Antibody)) +
  geom_density(alpha = 0.2) +
  scale_x_log10("x",
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand = c(0, 0)) +
  facet_grid(~clst) +
  theme_bw() +
  coord_fixed()+
  labs(x = "", y = "Density") +
  theme(axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # legend.position = c(0.8, 0.95),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        strip.background = element_blank(),
        text = element_text(size = 12, face = 2))

ggsave("3.Output/revision/5.Nrp1ProteinExpressionTest/Treg_cells_Protein_expression_comparison_between_pos_neg.png",
       width = 8, height = 3, units = "in", dpi = 300, bg = "white")

# add meta data ----

# dat <- AddMetaData(dat, protDf %>% select(clst) %>% rename(nrp1Prot = clst) %>% 
#   merge(data.frame(barcode = Cells(dat)), by.x = 0, by.y = "barcode", all = T) %>% 
#   column_to_rownames(var = "Row.names"))

# split figure ----
# dat <- subset(dat, cellType == "Treg")
saveRDS(dat, "1.Data/revision/datwithNrp1ProteinExpStat.RDS")


saveRDS(subset(dat, cellType == "Treg"), 
        "1.Data/revision/datTregwithNrp1ProteinExpStat.RDS")
