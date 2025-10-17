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
serum.a <- read.csv("./expr.serum.angII.log2.csv") %>% column_to_rownames(var = "X")

## III. colors
treat_color <- c("#fe8b29", "#4fa8f8", "#c65911", "#305496")
names(treat_color) <- c("Ctrl_AL","Ctrl_CR","AngII_AL","AngII_CR")

tissue_color <- c("#BEBADA", "#8DD3C7", "#80B1D3")
names(tissue_color) <- c("Liver", "Heart", "Serum")

class_color <- c("#6282bd","#d66e2f","#dc8a8c","#e5b46f",
                 "#a4c3d8","#4a8f37","#1ad1c4","#d9bdd2","grey")
names(class_color) <- c("Carbohydrates","Lipids","Heterocyclic Compounds",
                        "Carboxylic Acids","Peptides","Vitamins and cofactors",
                        "Nucleotides","Others","Unknown")

########################### Figure S4A #############################
h=round(40/817*100, 2)
l=round(77/1003*100, 2)
s=round(65/658*100, 2)
data <- data.frame(tissue=c("Heart", "Liver", "Serum"), 
                   prop=c(h,l,s))
ggbarplot(data, x="tissue", y="prop", label = T, fill = "tissue")+
  scale_fill_manual(values = tissue_color)+
  xlab("")+ylab("Proportion of diff metabo/all metabo (%)")

########################### Figure 4A #############################
data <- t(as.matrix(serum.a))
coldata <- data.frame(Sample=rownames(data),
                      Treat=c(rep("Ctrl_AL",5),rep("Ctrl_CR",5),
                              rep("AngII_AL",5),rep("AngII_CR",5)))
coldata$Treat <- factor(coldata$Treat, levels = c("Ctrl_AL","Ctrl_CR","AngII_AL","AngII_CR"))

library(ropls)
data.plsda <- opls(data, coldata$Treat)
data.plot <- as.data.frame(data.plsda@scoreMN)
data.plot$group = coldata$Treat
data.plot$samples = rownames(data.plot)
x_lab <- data.plsda@modelDF[1, "R2X"] * 100
y_lab <- data.plsda@modelDF[2, "R2X"] * 100

ggplot(data.plot, aes(p1, p2, color = group)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5) + 
  geom_vline(xintercept = 0, linetype = 'dashed', size = 0.5) +
  geom_point() + 
  labs(x=paste0("P1 (",x_lab,"%)"), y=paste0("P2 (",y_lab,"%)"))+
  stat_ellipse(level = 0.95, linetype = 'solid',
               size = 0.5, show.legend = FALSE) +
  scale_color_manual(values = treat_color) +
  ggtitle("Serum_AngII")+
  theme_bw()

########################### Figure 4B #############################
diff.serum.ang <- read.csv("./diff.metabo.serum.ang.class.csv")
diff.expr <- serum.a[str_to_upper(rownames(serum.a)) %in% str_to_upper(diff.serum.ang$metabolites), ]
diff.expr <- select(diff.expr, serum_AL_A_1:serum_CR_A_5)

anno_col <- data.frame(Treat=c(rep("AngII_AL",5),rep("AngII_CR",5)))
anno_col$Treat <- factor(anno_col$Treat, levels = c("AngII_AL","AngII_CR"))
rownames(anno_col) <- colnames(diff.expr)

anno_row <- data.frame(metabolites=rownames(diff.expr))
anno_row$metabolites <- str_to_title(anno_row$metabolites)
diff.serum.ang$metabolites <- str_to_title(diff.serum.ang$metabolites)
anno_row <- left_join(anno_row, diff.serum.ang, "metabolites")

anno_row <- select(anno_row, metabolites, Define)
colnames(anno_row)[2] <- "Class"
anno_row <- column_to_rownames(anno_row, "metabolites")
anno_row$Class <- factor(anno_row$Class, levels = names(class_color))

pheatmap(diff.expr, scale = "row", show_rownames = T, show_colnames = F,
         cluster_rows = T, cluster_cols = F,
         color = colorRampPalette(c("#002f5f","#1565ab","#3894bf","#8cc5e2","#d1e4ef",
                                    "#f9f9f9","#ffdbcd","#fba486","#dc5e51","#b91234","#6e001c"))(100),
         cutree_rows = 2,
         gaps_col = 5,
         annotation_col = anno_col,
         annotation_row = anno_row,
         annotation_colors = list(Treat=treat_color[3:4], Class=class_color))

########################### Figure 4D #############################

# pathway enrichment results from MetaboAnalyst5.0 (https://www.metaboanalyst.ca/)
path <- read.csv("./result/serum_ang_pathway_results.csv")
path <- path %>% select(X, logp, Impact)
colnames(path)[1] <- "Path"
path$group <- "Serum+AngII"

path2 <- read_xlsx("./result/Serum_pathway_results.xlsx") %>% select(Path, logp, Impact)
path2$group <- "Serum"

res <- rbind(path2, path)
res$Path <- factor(res$Path, levels = c("Lysine degradation", 
                                        "Glycine, serine and threonine metabolism",
                                        "Glycerolipid metabolism", "Sphingolipid metabolism",
                                        "Histidine metabolism","Glycerophospholipid metabolism"))
ggplot(res, aes(x=group,y=Path))+
  geom_point(aes(size=Impact,color=logp))+
  theme_bw()+
  labs(x="",y="")+
scale_color_gradient("logp",low = "#e5c9c8", high = "#e70022")+
  scale_size_continuous(name="RichFactor",
                        breaks = c(0.1,0.2,0.3,0.4,0.5))

########################### Figure 4E #############################
inter <- intersect(str_to_upper(diff.serum$metabolites), str_to_upper(diff.serum.ang$metabolites))
diff.expr <- serum.a[str_to_upper(rownames(serum.a)) %in% inter, ]

diff.expr.serum <- as.data.frame(t(scale(t(diff.expr[,1:10]))))
diff.expr.ang <- as.data.frame(t(scale(t(diff.expr[,11:20]))))
diff.expr <- cbind(diff.expr.serum, diff.expr.ang)

anno_col <- data.frame(Treat=c(rep("Ctrl_AL",5),rep("Ctrl_CR",5),
                               rep("AngII_AL",5),rep("AngII_CR",5)))
anno_col$Treat <- factor(anno_col$Treat, levels = names(treat_color))
rownames(anno_col) <- colnames(diff.expr)

anno_row <- data.frame(metabolites=rownames(diff.expr))
anno_row$metabolites <- str_to_title(anno_row$metabolites)
diff.serum.ang$metabolites <- str_to_title(diff.serum.ang$metabolites)
anno_row <- left_join(anno_row, diff.serum.ang, "metabolites")

anno_row <- select(anno_row, metabolites, Define)
colnames(anno_row)[2] <- "Class"
anno_row <- column_to_rownames(anno_row, "metabolites")
anno_row$Class <- factor(anno_row$Class, levels = names(class_color))

pheatmap(diff.expr, show_rownames = T, show_colnames = F,
         cluster_rows = T, cluster_cols = F,
         color = colorRampPalette(c("#002f5f","#1565ab","#3894bf","#8cc5e2","#d1e4ef",
                                    "#f9f9f9","#ffdbcd","#fba486","#dc5e51","#b91234","#6e001c"))(100),
         main = "DEM of Serum AL vs CR overlap with Serum_AngII AL vs CR",
         gaps_col = c(5,10,15,20),
         annotation_col = anno_col,
         annotation_row = anno_row,
         annotation_colors = list(Treat=treat_color, Class=class_color))

########################### Figure 4G & S4C #############################
l.expr <- serum.a["Lpe(16:0/0:0)", ]
rownames(l.expr) <- str_to_upper(rownames(l.expr))
l.expr <- gather(l.expr, key="Sample", value = "Expr")
l.expr$Group <- substr(l.expr$Sample, start = 7, stop = 10)
names(treat_color) <- c("AL_O","CR_O", "AL_A", "CR_A")
compare <- list(c("AL_O","CR_O"), c("AL_O", "AL_A"), c("AL_A", "CR_A"))
ggboxplot(l.expr, x="Group", y="Expr", color = "Group", add = "jitter")+
  scale_color_manual(values=treat_color)+
  stat_compare_means(comparisons = compare, label = "p.signif")+
  xlab("")+ylab("Expression of LPE(16:0/0:0)")

l.expr <- serum.a["Citrulline", ]
rownames(l.expr) <- str_to_upper(rownames(l.expr))
l.expr <- gather(l.expr, key="Sample", value = "Expr")
l.expr$Group <- substr(l.expr$Sample, start = 7, stop = 10)
names(treat_color) <- c("AL_O","CR_O", "AL_A", "CR_A")
compare <- list(c("AL_O","CR_O"), c("AL_O", "AL_A"), c("AL_A", "CR_A"))
ggboxplot(l.expr, x="Group", y="Expr", color = "Group", add = "jitter")+
  scale_color_manual(values=treat_color)+
  stat_compare_means(comparisons = compare, label = "p.signif")+
  xlab("")+ylab("Expression of Citrulline")


