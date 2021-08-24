setwd("C:/Users/postu003/Dropbox/Wageningen/epiGBS/finalbenchmarkingepigbs/")
library("ggplot2")
library("ggpubr")
library("ggrepel")
#args <- commandArgs(trailingOnly=T)

# epi-platypus
data1a <- "results/snp_calling/denovo/QUAL_Filt_squash/weighted_roc.tsv.gz"
data1a <- read.table(data1a,header=F)
data1a$tool <- "denovo-filtered"

# epi-platypus
data1b <- "results/snp_calling/denovo/QUAL_noFilt_squash/weighted_roc.tsv.gz"
data1b <- read.table(data1b,header=F)
data1b$tool <- "denovo-no filter"

# epi-freebayes
data2a <- "results/snp_calling/reference/QUAL_Filt_squash/weighted_roc.tsv.gz"
data2a <- read.table(data2a,header=F)
data2a$tool <- "reference-filtered"

# epi-freebayes
data2b <- "results/snp_calling/reference/QUAL_noFilt_squash/weighted_roc.tsv.gz"
data2b <- read.table(data2b,header=F)
data2b$tool <- "reference-no filter"

# epi-gatk
data3a <- "results/snp_calling/denovo/QUAL_Filt_squash_Col0/weighted_roc.tsv.gz"
data3a <- read.table(data3a,header=F)
data3a$tool <- "Col0 denovo-filtered"
# epi-gatk
data3b <- "results/snp_calling/denovo/QUAL_noFilt_squash_Col0/weighted_roc.tsv.gz"
data3b <- read.table(data3b,header=F)
data3b$tool <- "Col0 denovo-no filter"
# biscuit

denovo <- rbind(data1a,data1b)
colnames(denovo) <- c("score","tpb","fp","tpc","fn","precision","recall","fscore","tool")
df1 <- as.data.frame(denovo)
#df$tool <- as.factor(df$tool)
df1$tool <- as.factor(sub("^.*-","",df1$tool))
hist(df1$score)
  
reference <- rbind(data2a,data2b)
colnames(reference) <- c("score","tpb","fp","tpc","fn","precision","recall","fscore","tool")
df2 <- as.data.frame(reference)
#df$tool <- as.factor(df$tool)
df2$tool <- as.factor(sub("^.*-","",df2$tool))


denovo_col <- rbind(data3a,data3b)
colnames(denovo_col) <- c("score","tpb","fp","tpc","fn","precision","recall","fscore","tool")
df3 <- as.data.frame(denovo_col)
#df$tool <- as.factor(df$tool)
df3$tool <- as.factor(sub("^.*-","",df3$tool))



str(df1)
head(df1)
str(df2)
head(df2)

# manipulation
df1 <- df1[df1$recall>=0.01&df1$score>=10,]
df2 <- df2[df2$recall>=0.01&df2$score>=10,]

scoredf1<-df1[c(round(seq(1,length(df1[df1$tool=="filtered",1]),length.out=6)),round(seq(1,length(df1[df1$tool=="no filter",1]),length.out=6))+length(df1[df1$tool=="filtered",1])),]

scoredf2<-df2[c(seq(1,length(df2[df2$tool=="filtered",1]),length.out=6),round(seq(1,length(df2[df2$tool=="no filter",1]),length.out=6))+length(df2[df2$tool=="filtered",1])),]

# ggplot2
library(ggplot2)
library(ggpubr)

#cbPalette <- c("#56B4E9","#009E73","#CC79A7","#F0E422","#000000","#D55E00","#0072B2","#E69F00")
#cbPalette <- c("#CC79A7","#999999","#D55E00","#0072B2","#000000","#E69F00","#009E73")
cbPalette <- c("Black","Blue","#000000","#0072B2","#E69F00")

pr1 <- ggplot(df1,aes(x=recall,y=precision,color=tool)) +
  geom_step(size=0.9) +
  #geom_text_repel(data = scoredf1,aes(x=recall,y=precision,color=tool,label=round(score)),size=2)+
  ylab("Precision") +
  xlab("Sensitivity") +
  scale_y_continuous(breaks=seq(0, 1, 0.2))+
  scale_x_continuous(breaks=seq(0, 1, 0.2))+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  scale_color_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  #scale_linetype_manual(values=c(6,3,2,1,2,3,6,1,4))+
  scale_linetype_manual(values=c(3,4,5,6,1,1,1))+
  theme_bw() +
  theme(legend.title=element_blank(), legend.text = element_text(size = 11), legend.key.width = unit(1,"cm")) +
  guides(color=guide_legend(nrow=1),linetype=guide_legend(nrow=1,override.aes = list(size=0.8)))
  
pr1

pr2 <- ggplot(df2,aes(x=recall,y=precision,color=tool)) +
  geom_step(size=0.9) +
  #geom_text_repel(data = scoredf2,aes(x=recall,y=precision,color=tool,label=round(score)),size=2)+
  ylab("Precision") +
  xlab("Sensitivity") +
  scale_y_continuous(breaks=seq(0, 1, 0.2))+
  scale_x_continuous(breaks=seq(0, 1, 0.2))+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  scale_color_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  #scale_linetype_manual(values=c(6,3,2,1,2,3,6,1,4))+
  scale_linetype_manual(values=c(3,4,5,6,1,1,1))+
  theme_bw() +
  theme(legend.position="none")
pr2

pr3 <- ggplot(df3[df3$score>10,],aes(x=recall,y=precision,color=tool)) +
  geom_step(size=0.9) +
  ylab("Precision") +
  xlab("Sensitivity") +
  scale_y_continuous(breaks=seq(0, 1, 0.2))+
  scale_x_continuous(breaks=seq(0, 1, 0.2))+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  scale_color_manual("",values=cbPalette)+
  scale_fill_manual("",values=cbPalette)+
  #scale_linetype_manual(values=c(6,3,2,1,2,3,6,1,4))+
  scale_linetype_manual(values=c(3,4,5,6,1,1,1))+
  theme_bw()+
  guides(color=guide_legend(nrow=1),linetype=guide_legend(nrow=1,override.aes = list(size=0.8)))
  
pr3
ggsave("")

g <- ggarrange(pr2,pr1, ncol=2, nrow=1, widths=c(6,6), align="h", labels="AUTO", common.legend = TRUE, legend="top", legend.grob=get_legend(pr3))
g
ggsave("../epiGBS_paper/Figures/pdf-Figures/Figure7.pdf",g,height=4,width=9,dpi="retina")
?ggarrange

p<- ggarrange(pr2,pr1,pr3, ncol=3, nrow=1, widths=c(6,6), align="h", labels="AUTO", common.legend = TRUE, legend="top", legend.grob=get_legend(pr3))
p

df2[c(20660:20668),]


######################




pr12 <- ggplot(df2,aes(x=tpc,y=fp,color=tool)) +
  geom_step(size=0.9) +
  geom_text_repel(data = scoredf2,aes(x=tpc,y=fp,color=tool,label=round(score)),size=2)+
  ylab("False Positives") +
  xlab("True Positives") +
  scale_color_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  #scale_linetype_manual(values=c(6,3,2,1,2,3,6,1,4))+
  scale_linetype_manual(values=c(3,4,5,6,1,1,1))+
  theme_bw() +
  theme(legend.position="none")
pr12


pr11 <- ggplot(df1,aes(x=tpc,y=fp,color=tool)) +
  geom_step(size=0.6) +
  geom_text_repel(data = scoredf1,aes(x=tpc,y=fp,color=tool,label=round(score)),size=2)+
  ylab("False Positives") +
  xlab("True Positives") +
  scale_color_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  #scale_linetype_manual(values=c(6,3,2,1,2,3,6,1,4))+
  scale_linetype_manual(values=c(3,4,5,6,1,1,1))+
  theme_bw() +
  theme(legend.position="none")
pr11


g <- ggarrange(pr12,pr11, ncol=2, nrow=1, widths=c(6,6), align="h", labels="AUTO", common.legend = TRUE, legend="top", legend.grob=get_legend(pr3))
g
ggsave("../epiGBS_paper/Figures/SNP_calling_FP.tiff",g,height=4,width=9,dpi="retina")




    