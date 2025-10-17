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

class_color <- c("#6282bd","#d66e2f","#dc8a8c","#e5b46f",
                 "#a4c3d8","#4a8f37","#1ad1c4","#d9bdd2","grey")
names(class_color) <- c("Carbohydrates","Lipids","Heterocyclic Compounds",
                        "Carboxylic Acids","Peptides","Vitamins and cofactors",
                        "Nucleotides","Others","Unknown")

########################### Figure 2A-C #############################
# 1. serum
data <- t(as.matrix(tissues.expr$serum))
coldata <- data.frame(Sample=rownames(data), Treat=c(rep("AL",5),rep("CR",5)))
coldata$Treat <- factor(coldata$Treat, levels = c("AL","CR"))

# 2. heart
data <- t(as.matrix(tissues.expr$heart))
coldata <- data.frame(Sample=rownames(data), Treat=c(rep("AL",4),rep("CR",5)))
coldata$Treat <- factor(coldata$Treat, levels = c("AL","CR"))

# 3. liver
data <- t(as.matrix(tissues.expr$liver))
coldata <- data.frame(Sample=rownames(data), Treat=c(rep("AL",4),rep("CR",6)))
coldata$Treat <- factor(coldata$Treat, levels = c("AL","CR"))

library(ropls)

data.oplsda <- opls(data, coldata$Treat, predI = 1, orthoI = NA)

sample.score = data.oplsda@scoreMN %>%
  as.data.frame() %>%
  mutate(Treat=coldata$Treat,
         o1 = data.oplsda@orthoScoreMN[,1])
head(sample.score)

ggplot(sample.score, aes(p1, o1, color = Treat)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5) + 
  geom_vline(xintercept = 0, linetype = 'dashed', size = 0.5) +
  geom_point() + 
  labs(x = 'Component 1',y = 'Component 2') +
  stat_ellipse(level = 0.95, linetype = 'solid',
               size = 0.5, show.legend = FALSE) +
  scale_color_manual(values = group_color) +
  ggtitle("Liver")+
  #ggtitle("Heart")+
  #ggtitle("Serum")+
  theme_bw()

############################# Figure 2F-H #############################
heart.class <- read_xlsx("./diff.metabo.heart.class.xlsx") %>% select(metabolites, Define)
liver.class <- read_xlsx("./diff.metabo.liver.class.xlsx") %>% select(metabolites, Define)
serum.class <- read_xlsx("./diff.metabo.serum.class.xlsx") %>% select(metabolites, Define)

h.diff <- heart.expr[str_to_upper(rownames(heart.expr)) %in% 
                       str_to_upper(heart.class$metabolites), ]
l.diff <- liver.expr[str_to_upper(rownames(liver.expr)) %in% 
                       str_to_upper(liver.class$metabolites), ]
s.diff <- serum.expr[str_to_upper(rownames(serum.expr)) %in% 
                       str_to_upper(serum.class$metabolites), ]

# 1. liver
anno_col <- data.frame(Treat=c(rep("AL",5),rep("CR",6)))
anno_col$Treat <- factor(anno_col$Treat, levels = c("AL","CR"))
rownames(anno_col) <- colnames(liver.expr)

tmp=rownames_to_column(l.diff, "metabolites")
tmp$metabolites=str_to_upper(tmp$metabolites)
liver.class$metabolites=str_to_upper(liver.class$metabolites)
anno_row <- tmp %>% left_join(., liver.class, by="metabolites")
anno_row$metabolites <- str_to_title(anno_row$metabolites)
anno_row <- select(anno_row, metabolites, Define)
colnames(anno_row)[2] <- "Class"
anno_row <- column_to_rownames(anno_row, "metabolites")
anno_row$Class <- factor(anno_row$Class, levels = names(class_color))

pheatmap(l.diff, scale = "row", show_rownames = T, show_colnames = F,
         cluster_rows = T, cluster_cols = T,
         color = colorRampPalette(c("#002f5f","#1565ab","#3894bf","#8cc5e2","#d1e4ef",
                                    "#f9f9f9","#ffdbcd","#fba486","#dc5e51","#b91234","#6e001c"))(100),
         cutree_rows = 2,
         cutree_cols = 2,
         annotation_col = anno_col,
         annotation_row = anno_row,
         annotation_colors = list(Treat=group_color, Class=class_color))

# 2. heart
anno_col <- data.frame(Treat=c(rep("AL",4),rep("CR",5)))
anno_col$Treat <- factor(anno_col$Treat, levels = c("AL","CR"))
rownames(anno_col) <- colnames(heart.expr)

tmp=rownames_to_column(h.diff, "metabolites")
tmp$metabolites=str_to_upper(tmp$metabolites)
heart.class$metabolites=str_to_upper(heart.class$metabolites)
anno_row <- tmp %>% left_join(., heart.class, by="metabolites")
anno_row$metabolites <- str_to_title(anno_row$metabolites)
anno_row <- select(anno_row, metabolites, Define)
colnames(anno_row)[2] <- "Class"
anno_row <- column_to_rownames(anno_row, "metabolites")
anno_row$Class <- factor(anno_row$Class, levels = names(class_color))

pheatmap(h.diff, scale = "row", show_rownames = T, show_colnames = F,
         cluster_rows = T, cluster_cols = F,
         color = colorRampPalette(c("#002f5f","#1565ab","#3894bf","#8cc5e2","#d1e4ef",
                                    "#f9f9f9","#ffdbcd","#fba486","#dc5e51","#b91234","#6e001c"))(100),
         cutree_rows = 2,
         gaps_col = 4,
         annotation_col = anno_col,
         annotation_row = anno_row,
         annotation_colors = list(Treat=group_color, Class=class_color))

# 3. serum
anno_col <- data.frame(Treat=c(rep("AL",5),rep("CR",5)))
anno_col$Treat <- factor(anno_col$Treat, levels = c("AL","CR"))
rownames(anno_col) <- colnames(serum.expr)

tmp=rownames_to_column(s.diff, "metabolites")
tmp$metabolites=str_to_upper(tmp$metabolites)
serum.class$metabolites=str_to_upper(serum.class$metabolites)
anno_row <- tmp %>% left_join(., serum.class, by="metabolites")
anno_row$metabolites <- str_to_title(anno_row$metabolites)
anno_row <- select(anno_row, metabolites, Define)
colnames(anno_row)[2] <- "Class"
anno_row <- column_to_rownames(anno_row, "metabolites")
anno_row$Class <- factor(anno_row$Class, levels = names(class_color))

pheatmap(s.diff, scale = "row", show_rownames = T, show_colnames = F,
         cluster_rows = T, cluster_cols = T,
         color = colorRampPalette(c("#002f5f","#1565ab","#3894bf","#8cc5e2","#d1e4ef",
                                    "#f9f9f9","#ffdbcd","#fba486","#dc5e51","#b91234","#6e001c"))(100),
         cutree_rows = 2,
         cutree_cols = 2,
         annotation_col = anno_col,
         annotation_row = anno_row,
         annotation_colors = list(Treat=group_color, Class=class_color))

############################# Figure S2A #############################
data <- data.frame(Heart=length(rownames(tissues.expr$heart)),
                   Liver=length(rownames(tissues.expr$liver)), 
                   Serum=length(rownames(tissues.expr$serum)))
data <- gather(data, key = "Tissue", value = "Num")

ggbarplot(data, x="Tissue", y="Num", label = T, fill = "Tissue")+
  scale_fill_manual(values = tissue_color)+
  xlab("")+ylab("Number of metabolites")

############################# Figure S2B-D #############################

# 1. Liver
liver <- tissues.expr$liver  
anno_col <- data.frame(Treat=c(rep("AL",5),rep("CR",6)))
anno_col$Treat <- factor(anno_col$Treat, levels = c("AL", "CR"))
rownames(anno_col) <- colnames(liver)
pheatmap(liver, scale = "row", show_rownames = F, show_colnames = T,
         cluster_rows = T, cluster_cols = F,
         color = colorRampPalette(c("#002f5f","#1565ab","#3894bf","#8cc5e2","#d1e4ef",
                                    "#f9f9f9","#ffdbcd","#fba486","#dc5e51","#b91234","#6e001c"))(100),
         annotation_col = anno_col,
         annotation_colors = list(Treat=group_color))
# 2. Heart
heart <- tissues.expr$heart
anno_col <- data.frame(Treat=c(rep("AL",4),rep("CR",5)))
rownames(anno_col) <- colnames(heart)
pheatmap(heart, scale = "row", show_rownames = F, show_colnames = T,
         cluster_rows = T, cluster_cols = F,
         color = colorRampPalette(c("#002f5f","#1565ab","#3894bf","#8cc5e2","#d1e4ef",
                                    "#f9f9f9","#ffdbcd","#fba486","#dc5e51","#b91234","#6e001c"))(100),
         annotation_col = anno_col,
         annotation_colors = list(Treat=group_color))
# 3. Serum
serum <- tissues.expr$serum
anno_col <- data.frame(Treat=c(rep("AL",5),rep("CR",5)))
rownames(anno_col) <- colnames(serum)
pheatmap(serum, scale = "row", show_rownames = F, show_colnames = T,
         cluster_rows = T, cluster_cols = F,
         color = colorRampPalette(c("#002f5f","#1565ab","#3894bf","#8cc5e2","#d1e4ef",
                                    "#f9f9f9","#ffdbcd","#fba486","#dc5e51","#b91234","#6e001c"))(100),
         annotation_col = anno_col,
         annotation_colors = list(Treat=group_color))
