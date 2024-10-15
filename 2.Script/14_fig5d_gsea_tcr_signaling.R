# ----------------------------------------------------------------------------
# Fig 5d | GSEA of Nrp1-positive intra-islet Tregs against TCR signalling gene sets
#
# Paper    : Treg insufficiency in autoimmune diabetes is driven by selective
#            loss of Neuropilin-1 on intra-islet regulatory T cells
# Original : 2.Script/gseaplot.r
# Input    : 1.Data/initial/4.newDatTreg.RDS
#            1.Data/initial/Nrp1 altered genes.csv
# Output   : 3.Output/initial/GseaFigures/
#            3.Output/initial/gseaBar.png
# Note     : Uses the MSigDB TREG_VS_TCONV gene sets
#            also reads two published supplementary tables (see README)
# ----------------------------------------------------------------------------

library(Seurat)
library(msigdbr)
library(fgsea)
library(tidyverse)
library(future)
library(presto)
library(SPIA)
library(tidyverse)
library(KEGGREST)
library(readxl)

if(!dir.exists("3.Output/initial/GseaFigures")){
  dir.create("3.Output/initial/GseaFigures")
}

dat <- readRDS("1.Data/initial/4.newDatTreg.RDS")
Idents(dat) <- "type"
deg <- wilcoxauc(dat, "type") %>%
  # filter(pval < 0.05) %>% 
  filter(group == "L/+") %>% 
  arrange(desc(auc)) %>% 
  select(feature, auc, logFC) 

write.csv(deg, '3.Output/initial/GseaFigures/deg.csv')
# genes from other source ----

## Greg up genes
df <- read.csv("1.Data/initial/Nrp1 altered genes.csv", header = F)
pathway <- df$V1

## Abby's

df <- read_excel("1.Data/initial/Nrp1 Fl vs WT Til Tregs abbys paper.xlsx",
                 sheet = 1)

names(df)[1] <- "gene"

df <- df[df$`log2 fold change` < 0,]


## Maria 
deg <- deg %>% filter(logFC != 0)
excel_sheets("1.Data/initial/ITRhi_vs_ITRlo.xlsx")
df <- read_excel("1.Data/initial/ITRhi_vs_ITRlo.xlsx", 
                 sheet = 4)



df <- df[df$`log2(fold_change)` > 0,]

# GSEA
library(msigdbr)

m_df<- msigdbr(species = "Mus musculus", category = "C7")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name) 

fgsea_sets <- fgsea_sets[grep("TREG_VS_TCONV", names(fgsea_sets))]

deg <- deg %>% filter(logFC != 0)

fgseaRes <- fgsea(pathways = fgsea_sets,
                  stats = deg %>% 
                    arrange(desc(auc)) %>% 
                    select(feature, auc) %>% 
                    tibble::deframe() ,
                  minSize=15,
                  maxSize=500)


fgseaResSub <- fgseaRes[grep("_DN", fgseaRes$pathway),]
fgseaResSub <- fgseaResSub[fgseaResSub$NES < 0,]

data.table::fwrite(fgseaRes, "tregVsTconv.csv")



pathway <- df$gene[!is.na(df$gene)]

df$gseaAdjp <- NA

df$gseaNES <- NA

dfRes <- fgsea(pathways = list(p = pathway),
               stats = deg %>% 
                 # filter(logFC != 0) %>% 
                 arrange(desc(auc)) %>% 
                 select(feature, auc) %>% 
                 tibble::deframe(),
               minSize=15,
               maxSize=500)

df$gseaAdjp <- dfRes$padj
df$gseaNES <- dfRes$NES

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

ggplot(toPlotNew, aes(x = x, y = y, alpha = 0.7)) + 
  geom_point(size = 1, shape = 16) +
  geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") + 
  geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed")+ 
  ggtitle("DN in CD5hi Tregs", 
          subtitle = paste0('adjusted P value: ', 
                            format(dfRes$padj, digits = 2, scientific = T), 
                            ", NES: ", 
                            format(dfRes$NES, digits = 2))) +
  theme_bw() + 
  geom_segment(aes(x = x, y = -diff/4, xend = x, yend = diff/4), col = "blue") + 
  theme(legend.position = "none",
        text = element_text(size = 15, face = 2)) + 
  labs(x = "Rank", 
       y = "Enrichment score")
ggsave(paste0( "3.Output/initial/GseaFigures/", "MariaDN", ".png"), width = 6, height = 4,
       units = "in", dpi = 300)

# spia ----
df <- read.csv("3.Output/initial/3.PseudoBulk/islet_Treg_SPIA_KEGG_by_LogFC(logFCnot0).csv")
df <- df[df$X %in% c(14, 16, 20, 21, 31, 40, 42, 48, 57, 62, 82, 86, 87, 101, 122, 123),]
df <- df[df$X %in% c(14, 16, 20, 40),]

df$gseaAdjp <- NA
df$gseaNES <- NA

for (i in seq(nrow(df))) {
  tryCatch({
    query <- keggGet(paste0("mmu0", df$ID[i]))
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
    
    df$gseaAdjp[i] <- dfRes$padj
    df$gseaNES[i] <- dfRes$NES
    
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
    
    
    ggplot(toPlotNew, aes(x = x, y = y, alpha = 0.7)) + 
      geom_point(size = 1, shape = 16) +
      geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") + 
      geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed")+ 
      ggtitle(df$Name[i], 
              subtitle = paste0('adjusted P value: ', 
                                format(dfRes$padj, digits = 2, scientific = T), 
                                ", NES: ", 
                                format(dfRes$NES, digits = 2))) +
      theme_bw() + 
      geom_segment(aes(x = x, y = -diff/4, xend = x, yend = diff/4), col = "blue") + 
      theme(legend.position = "none",
            text = element_text(size = 15, face = 2)) + 
      labs(x = "Rank", 
           y = "Enrichment score")
    # header <- ifelse(dfRes$ES * df$tA[i] < 0, "Diff", "")
    
    write.csv(toPlotNew, paste0("3.Output/initial/GseaFigures/", df$Name[i], ".csv"))
    
    # ggsave(paste0( "3.Output/initial/GseaFigures/", df$Name[i], ".png"), width = 6, height = 4,
    #        units = "in", dpi = 300)
  }, error=function(e){})
}

# bar ----
df$gseaStatus <- ifelse(df$gseaNES > 0, "Activated", "Inhibited")
ggplot(df, aes(x = gseaNES, y = reorder(Name, gseaNES), fill = gseaStatus)) +
  geom_bar(stat="identity") +
  xlab("Normalized enrichment score") +
  ylab("Pathway") +
  theme_minimal() +
  theme(text = element_text(size = 15, face = 2),
        legend.position = "none")

ggsave("3.Output/initial/gseaBar.png", width = 12, height = 4, units = "in", dpi = 300, bg = "white")



# pie ----
data <- data.frame(table(toPlotNew$logFC > 0))
data$category <- ifelse(data$Var1 == T, "UP", "DN")
data$fraction = data$Freq / sum(data$Freq)
data$ymax = cumsum(data$fraction)
data$ymin = c(0, head(data$ymax, n=-1))
data$labelPosition <- (data$ymax + data$ymin) / 2
data$label <- paste0(data$category, ": \n", data$Freq)

p2 <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category, alpha = 0.8)) +
  geom_bar() +
  geom_text( x=2, aes(y=labelPosition, label=label, color=category), size=6) + # x here controls label position (inner / outer)
  # coord_polar(theta="y") +
  scale_fill_manual(breaks = c("UP", "DN"), 
                    values=c("#D50A0A", "#013369")) +
  scale_color_manual(breaks = c("UP", "DN"), 
                    values=c("#D50A0A", "#013369")) +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none") 



p2 <- data %>%
  ggplot(aes(x = category,
             y = Freq,
             fill=category,
             label=Freq))+
  geom_col(width=0.5)+
  scale_fill_manual(breaks = c("UP", "DN"), 
                    values=c("#D50A0A", "#013369"))+
  theme(legend.position="none")+
  labs(x="",
       y="Mean Load Debt")+
  geom_text(nudge_y= -50,
            color="white",
            size = 5,
            fontface="bold")+
  theme_void()+
  theme(legend.position = "none") 

p1 + p2

require(gtable)
library(grid)

p1 <- ggplotGrob(p1)
p2 <- ggplotGrob(p2)

plot.new()
gt <- gtable(widths = unit(c(1, 1), "null"), heights = unit(c(1, 1, 1), "null"))
gt <- gtable_add_grob(gt, p1, t = 1, b = 3, l = 1, r = 2)
gt <- gtable_add_grob(gt, p2, t = 1, l = 2)
grid.draw(gt)

