## I. load packages
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(readxl)
library(rlist)
library(data.table)
library(pheatmap)
library(VennDiagram)
library(ggplotify)
library(cowplot)
library(pheatmap)
library(clusterProfiler)
library(RColorBrewer)
library(paletteer)

## II. load data
tissues.expr <- readRDS("./tissues.expr.rds")
heart.expr <- tissues.expr$heart
liver.expr <- tissues.expr$liver
serum.expr <- tissues.expr$serum

## III. colors
colAL<- "#fe8b29"
colCR<- "#4fa8f8"
group_color <- c(colAL, colCR)
names(group_color) <- c("AL", "CR")

tissue_color <- c("#BEBADA", "#8DD3C7", "#80B1D3")
names(tissue_color) <- c("Liver", "Heart", "Serum")

########################### Figure 3A #############################

# pathway enrichment results from MetaboAnalyst5.0 (https://www.metaboanalyst.ca/)
h.res <- read_xlsx("./heart_pathway_results.xlsx")
l.res <- read_xlsx("./liver_pathway_results.xlsx")
s.res <- read_xlsx("./serum_pathway_results.xlsx")

kegg.res <- rbind(h.res, l.res, s.res)
kegg.res$Tissue <- c(rep("Heart", 9), rep("Liver", 9), rep("Serum", 4))
kegg.res <- kegg.res[order(kegg.res$logp, decreasing = T), ]
kegg.res$Path <- factor(kegg.res$Path,
                        levels=rev(kegg.res$Path[!duplicated(kegg.res$Path)]))
kegg.res$Tissue <- factor(kegg.res$Tissue, levels = c("Heart","Liver","Serum"))
kegg.res$Ratio <- kegg.res$Hits/kegg.res$Total

kegg.res %>% filter()
ggplot(kegg.res,aes(x=Tissue,y=Path))+
  geom_point(aes(size=Impact,color=logp))+
  theme_bw()+
  labs(x="",y="")+
  scale_color_gradient("logp",low = "#e5c9c8", high = "#e70022")+
  scale_size_continuous(name="RichFactor",
                        breaks = c(0.1, 0.2,0.3,0.4,0.5))

########################### Figure 3C #############################
# 1. Heart
h.lipid <- c("PC(18:0/18:2)","PC(18:1/16:1)","PC(18:1/18:2)","PC(18:2/22:6)",
             "O-Phosphorylethanolamine","Glycerophosphocholine")
h.expr <- heart.expr[str_to_upper(rownames(heart.expr)) %in% str_to_upper(h.lipid), ]

anno_col <- data.frame(Treat=c(rep("AL",4),rep("CR",5)))
anno_col$Treat <- factor(anno_col$Treat, levels = c("AL","CR"))
rownames(anno_col) <- colnames(heart.expr)

pheatmap(h.expr, scale="row", cluster_rows = T, cluster_cols = F,
         show_rownames = T, show_colnames = T,
         color = colorRampPalette(c("#002f5f","#1565ab","#3894bf","#8cc5e2",
                                    "#d1e4ef", "#f9f9f9","#ffdbcd","#fba486","#dc5e51","#b91234","#6e001c"))(100),
         cutree_rows = 2,
         gaps_col = 4, 
         annotation_col = anno_col,
         annotation_colors = list(Treat=group_color))

# 2. Liver
l.lipid <- c("PE(18:1/20:3)","PE(18:1/20:5)","PE(18:1/22:5)","PE(18:1/22:6)",
             "PE(P-16:0/22:6)","PC(16:0/18:1)","PC(16:0/20:4)","PC(16:0/22:6)",
             "PC(18:0/14:0)","PC(18:0/18:1)","PC(18:0/18:2)","PC(18:0/22:6)",
             "PC(18:1/16:1)","PC(18:1/18:2)","PC(18:1/22:6)","PC(O-16:0/22:6)",
             "Acetylcholine","Choline","Phosphorylcholine")
l.expr <- liver.expr[str_to_upper(rownames(liver.expr)) %in% str_to_upper(l.lipid), ]

anno_col <- data.frame(Treat=c(rep("AL",5),rep("CR",6)))
anno_col$Treat <- factor(anno_col$Treat, levels = c("AL","CR"))
rownames(anno_col) <- colnames(liver.expr)

pheatmap(l.expr, scale="row", cluster_rows = T, cluster_cols = F,
         show_rownames = T, show_colnames = T,
         color = colorRampPalette(c("#002f5f","#1565ab","#3894bf","#8cc5e2",
                                    "#d1e4ef","#f9f9f9","#ffdbcd","#fba486","#dc5e51","#b91234","#6e001c"))(100),
         gaps_col = 5, 
         cutree_rows = 2,
         annotation_col = anno_col,
         annotation_colors = list(Treat=group_color))

# 3. Serum
s.lipid <- c("PE(18:1/20:3)","PE(P-16:0/22:5)","PE(P-16:0/22:6)",
             "PE(P-18:0/22:6)","PC(18:0/20:4)","PC(18:0/22:6)","PC(18:2/22:6)",
             "PC(O-16:0/22:6)", "Acetylcholine","O-Phosphorylethanolamine",
             "Glycerol 3-phosphate","Glycerophosphocholine")
s.expr <- serum.expr[str_to_upper(rownames(serum.expr)) %in% str_to_upper(s.lipid), ]

anno_col <- data.frame(Treat=c(rep("AL",5),rep("CR",5)))
anno_col$Treat <- factor(anno_col$Treat, levels = c("AL","CR"))
rownames(anno_col) <- colnames(serum.expr)

pheatmap(s.expr, scale="row", cluster_rows = T, cluster_cols = F,
         show_rownames = T, show_colnames = T,
         color = colorRampPalette(c("#002f5f","#1565ab","#3894bf","#8cc5e2",
                                    "#d1e4ef", "#f9f9f9","#ffdbcd","#fba486","#dc5e51","#b91234","#6e001c"))(100),
         gaps_col = 5, 
         annotation_col = anno_col,
         annotation_colors = list(Treat=group_color))

########################### Figure S3A #############################
h.diff <- read_xlsx("./diff.metabo.heart.class.xlsx")
l.diff <- read_xlsx("./diff.metabo.liver.class.xlsx")
s.diff <- read_xlsx("./diff.metabo.serum.class.xlsx")

h.l <- intersect(h.diff$metabolites, l.diff$metabolites)
l.s <- intersect(l.diff$metabolites, s.diff$metabolites)
h.s <- intersect(h.diff$metabolites, s.diff$metabolites)

########################### Figure S3E-G #############################
# 1. Heart-Liver
h.l.diff <- h.diff %>% inner_join(l.diff, "metabolites") %>% 
  select(metabolites, log2FC.AL.CR..x, log2FC.AL.CR..y)
colnames(h.l.diff)[2:3] <- c("logfc.heart", "logfc.liver")
h.l.diff$logfc.heart = -h.l.diff$logfc.heart
h.l.diff$logfc.liver = -h.l.diff$logfc.liver

ggscatter(h.l.diff, x="logfc.heart", y="logfc.liver", add = "reg.line",
          add.params = list(color = "blue", fill = "lightgray"), 
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n"),
          ggtheme = theme_bw())+ggtitle("Heart-Liver")

# 2. Serum-Heart
s.h.diff <- s.diff %>% inner_join(h.diff, "metabolites") %>% 
  select(metabolites, log2FC.AL.CR..x, log2FC.AL.CR..y)
colnames(s.h.diff)[2:3] <- c("logfc.serum", "logfc.heart")
s.h.diff$logfc.serum = -s.h.diff$logfc.serum
s.h.diff$logfc.heart = -s.h.diff$logfc.heart

ggscatter(s.h.diff, x="logfc.serum", y="logfc.heart", add = "reg.line",
          add.params = list(color = "blue", fill = "lightgray"), 
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n"),
          ggtheme = theme_bw())+ggtitle("Serum-Heart")

# 3. Serum-Liver
s.l.diff <- s.diff %>% inner_join(l.diff, "metabolites") %>% 
  select(metabolites, log2FC.AL.CR..x, log2FC.AL.CR..y)
colnames(s.l.diff)[2:3] <- c("logfc.serum", "logfc.liver")
s.l.diff$logfc.serum = -s.l.diff$logfc.serum
s.l.diff$logfc.liver = -s.l.diff$logfc.liver

ggscatter(s.l.diff, x="logfc.serum", y="logfc.liver", add = "reg.line",
          add.params = list(color = "blue", fill = "lightgray"), 
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n"),
          ggtheme = theme_bw())+ggtitle("Serum-Liver")

########################### Figure S3H #############################

# pathway enrichment results from MetaboAnalyst5.0 (https://www.metaboanalyst.ca/)
h.l <- read.csv("./result/path.diff.h.l.csv") %>% filter(Raw.p<0.05)
l.s <- read.csv("./result/path.diff.l.s.csv") %>% filter(Raw.p<0.05)
h.s <- read.csv("./result/path.diff.h.s.csv") %>% filter(Raw.p<0.05)

kegg.res <- rbind(h.l, h.s, l.s)
kegg.res$Tissue <- c(rep("Heart-Liver", 3), rep("Heart-Serum", 3), rep("Liver-Serum", 2))
kegg.res$logp <- -log10(kegg.res$Raw.p)
kegg.res <- kegg.res[order(kegg.res$logp, decreasing = T), ]
colnames(kegg.res)[1] <- "Path"
kegg.res$Path <- factor(kegg.res$Path,
                        levels=rev(kegg.res$Path[!duplicated(kegg.res$Path)]))
kegg.res$Tissue <- factor(kegg.res$Tissue, levels = c("Heart-Liver", "Heart-Serum", "Liver-Serum"))
kegg.res$Ratio <- kegg.res$Hits/kegg.res$Total

library(RColorBrewer)
library(paletteer)
ggplot(kegg.res, aes(x=Tissue,y=Path))+
  geom_point(aes(size=Impact,color=logp))+
theme_bw()+
  labs(x="",y="")+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust = 0.5))+
  scale_color_gradient("logp",low = "#e5c9c8", high = "#e70022")+
  scale_size_continuous(name="RichFactor",
                        breaks = c(0.05, 0.1, 0.15, 0.2))

