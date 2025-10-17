## I. load packages
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pROC)
library(glmnet)
library(lme4)
library(lmerTest)
library(GSEABase)

########################### Figure 5A #############################
expr <- fread("./public data/GSE140947_bulkrna_24_samp_genecounts.csv.gz",data.table = F)
expr %>% mutate(Symbol = str_remove(pattern = ".*_",gene_id)) ->expr
aggregate(.~Symbol, data =expr[,-1], mean)->GSE140947.expr.annot
GSE140947.expr.annot<- column_to_rownames(GSE140947.expr.annot,
                                          var = "Symbol")
GSE140947.expr.annot=GSE140947.expr.annot[rowSums(GSE140947.expr.annot)>20,]

data.plot$class <- factor(data.plot$class, levels = c("Control", "AA"))
mylogit <- glm(class ~ expr, 
               data = data.plot, family = "binomial")
summary(mylogit)
prob <- predict(mylogit,type='response')
data.plot$prob<- prob

g <- roc(class ~ prob,data = data.plot)
g
plot(g)

pROC::ggroc(g,legacy.axes = T,color="brown",
            alpha = 0.9)+
  ggpubr::theme_pubr()+
geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="darkgrey", linetype="dashed")+
  annotate("text", x = .8, y = .4, 
           label = paste0("GSE140947 AUC: ",round(g$auc,3)))


########################### Figure 5B #############################
data <- read.csv("./public data/GSE57691_data_matrix.csv") %>% column_to_rownames(var = "Name")

if(T){
  data1 <- gly.data %>% 
    filter(!grepl("AOD", sample)) %>% 
    filter(!grepl("large.AAA", sample)) %>% 
    mutate(class2 = if_else(grepl(pattern = "small.AAA",sample), "small.AAA", "Control"))
  data1$class2 <- factor(data1$class2, levels = c("Control","small.AAA"))
  
  mylogit <- glm(class2 ~ expr, 
                 data = data1, family = "binomial")
  summary(mylogit)
  prob <- predict(mylogit,type='response')
  data1$prob<- prob
  g1 <- roc(class2 ~ prob,data = data1)
  g1
  plot(g1)
}

#large.AAA
if(T){
  data2 <- gly.data %>% 
    filter(!grepl("AOD", sample)) %>% 
    filter(!grepl("small.AAA", sample)) %>% 
    mutate(class3 = if_else(grepl(pattern = "large.AAA",sample), "large.AAA", "Control"))
  data2$class3 <- factor(data2$class3, levels = c("Control","large.AAA"))
  
  mylogit <- glm(class3 ~ expr, 
                 data = data2, family = "binomial")
  summary(mylogit)
  prob <- predict(mylogit,type='response')
  data2$prob<- prob
  
  g2 <- roc(class3 ~ prob,data = data2)
  g2
  plot(g2)
}

#AOD
if(T){
  data3 <- gly.data %>% 
    filter(!grepl("large.AAA", sample)) %>% 
    filter(!grepl("small.AAA", sample)) %>% 
    mutate(class4 = if_else(grepl(pattern = "AOD",sample), "AOD", "Control"))
  data3$class4 <- factor(data3$class4, levels = c("Control","AOD"))
  
  mylogit <- glm(class4 ~ expr, 
                 data = data3, family = "binomial")
  summary(mylogit)
  prob <- predict(mylogit,type='response')
  data3$prob<- prob
  g3 <- roc(class4 ~ prob,data = data3)
  g3
  plot(g3)
}

pROC::ggroc(list(small.AAA=g1, large.AAA=g2, AOD=g3),
            legacy.axes = T, #color="brown", 
            alpha = 0.9)+
  ggpubr::theme_pubr()+
  #scale_colour_jco()+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="darkgrey", linetype="dashed")+
  annotate("text", x = .8, y = .5, 
           label = paste0("AUC (GSE57691-small.AAA): ",
                          round(as.numeric(g1$auc),3)))+
  annotate("text", x = .8, y = .4, 
           label = paste0("AUC (GSE57691-large.AAA): ",
                          round(as.numeric(g2$auc),3)))+
  annotate("text", x = .8, y = .3, 
           label = paste0("AUC (GSE57691-AOD): ",
                          round(as.numeric(g3$auc),3)))

########################### Figure S5A #############################

# load data
expr <- read_xlsx("./China BAPE-20190925.xlsx", sheet = 2) %>% 
  column_to_rownames(var = "Compound")

meta <- read_xlsx("./info.xlsx")
res <- left_join(meta, expr, by="Subject_ID")

expr <- as.data.frame(t(expr))
expr <- rownames_to_column(expr, var = "Subject_ID")

cor.res <- list()
for(i in chr){
  f <- formula(paste("xAge~",i,"+ (1|xsubject) + xGender + xBMI + xCotinine + xEducation + xIncome"))
  model <- lmer(f, data=res)
  
  lm.tmp.sum=summary(model)
  tmp2=lm.tmp.sum$coefficients[-1,]
  tmp2=as.data.frame(t(tmp2))
  tmp2$metabo <- i
  cor.res[[i]] <- tmp2
}

cor.res_combined <- bind_rows(cor.res)
age.res <- na.omit(cor.res_combined)
  
gly <- read_xlsx("./glys.xlsx")
age.res$class <- if_else(age.res$metabo %in% gly, "Glycerophospholipids", "Others")

age.res$logp <- -log10(age.res$`Pr(>|t|)`) 
data.plot <- age.res[order(age.res$logp, decreasing = T), ]
data.plot <- data.plot %>%
  arrange(desc(`logp`)) %>%
  mutate(rank = rank(`logp`, ties.method = "first"))
data.plot$label <- if_else(data.plot$`Pr(>|t|)`>0.05 | data.plot$class == "Others", "no", "yes")
data.plot <- data.plot %>% mutate(size = if_else(label == "no", 1.5, 2),
                                  alpha = if_else(label == "no", 0.3, 1))

ggplot(data.plot, aes(x=rank,y=logp,colour=label))+
  geom_point(size=data.plot$size, alpha=data.plot$alpha)+
  scale_colour_manual(values = c("gray","royalblue4"))+
  theme_pubr()+
  scale_size_continuous(range = c(1,3))+
  geom_hline(yintercept = 1.3,linetype="dotted")

########################### Figure 5E-F #############################

meta2 <- read.csv("./methyage.csv")
res <- left_join(res, meta2, by="Subject_ID")

expr <- as.data.frame(t(expr))
expr <- rownames_to_column(expr, var = "Subject_ID")

cor.res <- list()
for(i in chr){
  f <- formula(paste("xDNAmAge~",i,"+ (1|xsubject) + xGender + xBMI + xCotinine + xEducation + xIncome"))
  f <- formula(paste("xDNAmPhenoAge~",i,"+ (1|xsubject) + xGender + xBMI + xCotinine + xEducation + xIncome"))
  f <- formula(paste("xDNAmGrimAge~",i,"+ (1|xsubject) + xGender + xBMI + xCotinine + xEducation + xIncome"))
  model <- lmer(f, data=res)
  
  lm.tmp.sum=summary(model)
  tmp2=lm.tmp.sum$coefficients[-1,]
  tmp2=as.data.frame(t(tmp2))
  tmp2$metabo <- i
  cor.res[[i]] <- tmp2
}

cor.res_combined <- bind_rows(cor.res)
age.res <- na.omit(cor.res_combined)

gly <- read_xlsx("./glys.xlsx")
age.res$class <- if_else(age.res$metabo %in% gly, "Glycerophospholipids", "Others")

age.res$logp <- -log10(age.res$`Pr(>|t|)`) 
data.plot <- age.res[order(age.res$logp, decreasing = T), ]
data.plot <- data.plot %>%
  arrange(desc(`logp`)) %>%
  mutate(rank = rank(`logp`, ties.method = "first"))
data.plot$label <- if_else(data.plot$`Pr(>|t|)`>0.05 | data.plot$class == "Others", "no", "yes")
data.plot <- data.plot %>% mutate(size = if_else(label == "no", 1.5, 2),
                                  alpha = if_else(label == "no", 0.3, 1))

ggplot(data.plot, aes(x=rank,y=logp,colour=label))+
  geom_point(size=data.plot$size, alpha=data.plot$alpha)+
  scale_colour_manual(values = c("gray","royalblue4"))+
  theme_pubr()+
  scale_size_continuous(range = c(1,3))+
  geom_hline(yintercept = 1.3,linetype="dotted")

########################### Figure S5B #############################

meta3 <- read.csv("./rDNAS.csv")
res <- left_join(res, meta3, by="Subject_ID")
res <- na.omit(res)

expr <- as.data.frame(t(expr))
expr <- rownames_to_column(expr, var = "Subject_ID")

cor.res <- list()
for(i in chr){
  f <- formula(paste("xX18S~",i,"+(1|xsubject) + xGender + xBMI + xCotinine + xEducation + xIncome"))
  f <- formula(paste("xX5.8S~",i,"+(1|xsubject) + xGender + xBMI + xCotinine + xEducation + xIncome"))
  f <- formula(paste("xX5S~",i,"+(1|xsubject) + xGender + xBMI + xCotinine + xEducation + xIncome"))
  model <- lmer(f, data=res)
  
  lm.tmp.sum=summary(model)
  tmp2=lm.tmp.sum$coefficients[-1,]
  tmp2=as.data.frame(t(tmp2))
  tmp2$metabo <- i
  cor.res[[i]] <- tmp2
}

cor.res_combined <- bind_rows(cor.res)
age.res <- na.omit(cor.res_combined)

gly <- read_xlsx("./glys.xlsx")
age.res$class <- if_else(age.res$metabo %in% gly, "Glycerophospholipids", "Others")

age.res$logp <- -log10(age.res$`Pr(>|t|)`) 
data.plot <- age.res[order(age.res$logp, decreasing = T), ]
data.plot <- data.plot %>%
  arrange(desc(`logp`)) %>%
  mutate(rank = rank(`logp`, ties.method = "first"))
data.plot$label <- if_else(data.plot$`Pr(>|t|)`>0.05 | data.plot$class == "Others", "no", "yes")
data.plot <- data.plot %>% mutate(size = if_else(label == "no", 1.5, 2),
                                  alpha = if_else(label == "no", 0.3, 1))

ggplot(data.plot, aes(x=rank,y=logp,colour=label))+
  geom_point(size=data.plot$size, alpha=data.plot$alpha)+
  scale_colour_manual(values = c("gray","royalblue4"))+
  theme_pubr()+
  scale_size_continuous(range = c(1,3))+
  geom_hline(yintercept = 1.3,linetype="dotted")

########################### Figure 5G #############################
meta4 <- read.csv("./mtDNAS.csv")
res <- left_join(res, meta4, by="Subject_ID")
res <- na.omit(res)

expr <- as.data.frame(t(expr))
expr <- rownames_to_column(expr, var = "Subject_ID")

cor.res <- list()
for(i in chr){
  f <- formula(paste("xMT1.S~",i,"+(1|xsubject) + xGender + xBMI + xCotinine + xEducation + xIncome"))
  model <- lmer(f, data=res)
  
  lm.tmp.sum=summary(model)
  tmp2=lm.tmp.sum$coefficients[-1,]
  tmp2=as.data.frame(t(tmp2))
  tmp2$metabo <- i
  cor.res[[i]] <- tmp2
}

cor.res_combined <- bind_rows(cor.res)
age.res <- na.omit(cor.res_combined)

gly <- read_xlsx("./glys.xlsx")
age.res$class <- if_else(age.res$metabo %in% gly, "Glycerophospholipids", "Others")

age.res$logp <- -log10(age.res$`Pr(>|t|)`) 
data.plot <- age.res[order(age.res$logp, decreasing = T), ]
data.plot <- data.plot %>%
  arrange(desc(`logp`)) %>%
  mutate(rank = rank(`logp`, ties.method = "first"))
data.plot$label <- if_else(data.plot$`Pr(>|t|)`>0.05 | data.plot$class == "Others", "no", "yes")
data.plot <- data.plot %>% mutate(size = if_else(label == "no", 1.5, 2),
                                  alpha = if_else(label == "no", 0.3, 1))

ggplot(data.plot, aes(x=rank,y=logp,colour=label))+
  geom_point(size=data.plot$size, alpha=data.plot$alpha)+
  scale_colour_manual(values = c("gray","royalblue4"))+
  theme_pubr()+
  scale_size_continuous(range = c(1,3))+
  geom_hline(yintercept = 1.3,linetype="dotted")

########################### Figure 5I #############################
tmp = read.csv("./wzy/result/BAPE.gly.csscore.csv")
tmp %>%
  mutate(CSscore = case_when(score>quantile(tmp$score,.5)~"High",
                             score<quantile(tmp$score,.5)~"Low",
                             TRUE~"median")) %>% 
  filter(CSscore!="median") %>% 
  ggboxplot(x="CSscore",y="KEGG_GLYCEROPHOSPHOLIPID_METABOLISM",
            color = "CSscore",add = "jitter", palette = "npg")+
  scale_color_manual(values = c("Low" = "#6cb8d2", "High" = "#d55740"))+
  stat_compare_means()+xlab("")+ylab("Glycerophospholipid-metabolism pathway activity")

########################### Figure 5J #############################
gly <- read_xlsx("./China BAPE-20190925.xlsx", sheet = 3) %>% column_to_rownames(var = "Compound")
LPE16 <- t(gly)
rownames(LPE16) <- paste0("X", rownames(LPE16))
LPE16 <- rownames_to_column(as.data.frame(LPE16), var = "samples")

data <- left_join(tmp, LPE16, by="samples")

data %>%
mutate(CSscore = case_when(score>quantile(data$score,.5)~"High",
                             score<quantile(data$score,.5)~"Low",
                             TRUE~"median"
  )) %>% 
  filter(CSscore!="median") %>% 
  ggboxplot(x="CSscore",y="`LPE(16:0)+H*pos`",
            color = "CSscore",add = "jitter",palette = "npg")+
  scale_color_manual(values = c("Low" = "#6cb8d2", "High" = "#d55740"))+
  ylab("LPE(16:0)")+xlab("")+
  stat_compare_means()

########################### Figure S5C #############################
data <- read.csv("./GSE57691_data_matrix.csv") %>% column_to_rownames(var = "Name")

# cs score
aging_sg = readRDS("/Users/heping/Documents/文稿 - 和平的MacBook Pro/projects/aging/cancer aging/code/aging_signature.rds")
cs.score=GSEAScoring(aging_sg = aging_sg, method = "Gaussian", geneMat = data)
GSEAScoring<- function(aging_sg,geneMat,
                       method=c("Gaussian","Poisson"),parallel.sz=1){
  require(GSVA)
  pathways = list(up_path = names(aging_sg[aging_sg>0]),
                  down_path = names(aging_sg[aging_sg<0]))
  # count data: Poisson; tpm or fpkm data: Gaussian
  gsvapar <- gsvaParam(as.matrix(geneMat), kcdf = c("Gaussian"), geneSets=pathways, maxDiff=TRUE)
  tmp.score <- gsva(gsvapar)
  tmp.score=as.data.frame(t(tmp.score))
  tmp.score$score= tmp.score$up_path-tmp.score$down_path
  return(tmp.score)
}

my_bp<-getGmt("./KEGG_GLYCEROPHOSPHOLIPID_METABOLISM.v2024.1.Hs.gmt")
res <- gsvaParam(exprData = as.matrix(data), kcdf = c("Gaussian"), geneSets=my_bp, maxDiff=TRUE)
res <- gsva(res)

gly <- t(res)
gly <- rownames_to_column(as.data.frame(gly), var = "samples")
cs.score <- rownames_to_column(cs.score, var = "samples")
data <- left_join(cs.score, gly, by="samples")

data %>%
  mutate(CSscore = case_when(score>quantile(tmp$score,.5)~"High", 
                             score<quantile(tmp$score,.5)~"Low", 
                             TRUE~"median"
  ), CSscore = factor(CSscore, levels=c("Low","High"))) %>% 
  filter(CSscore != "median") %>% 
  ggboxplot(x="CSscore", y="KEGG_GLYCEROPHOSPHOLIPID_METABOLISM",
            color = "CSscore",add = "jitter",palette = "npg")+
  scale_color_manual(values = c("Low" = "#6cb8d2", "High" = "#d55740"))+
  xlab("")+ylab("Glycerophospholipid-metabolism pathway activity")+
  stat_compare_means()


