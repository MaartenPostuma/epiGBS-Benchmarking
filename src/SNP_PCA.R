rm(list=ls())
library(adegenet)
library(vcfR)
library(ggplot2)
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

                       
vcf.fnDenovo<-"results/snp_calling/denovo/testDenovo.vcf.gz"
vcfDenovo<-read.vcfR(vcf.fnDenovo)
genlightDenovo<-vcfR2genlight(vcfDenovo)


pcaDenovo<-glPca(genlightDenovo,nf = 10)



pcaPlotDenovo <- data.frame(pcaDenovo$scores,sample=row.names(pcaDenovo$scores),
                            pop=sub("_.*$","",row.names(pcaDenovo$scores)),ID=row.names(pcaDenovo$scores))

ggplot(pcaPlotDenovo,aes(x=PC1,y=PC2,col=pop,label=ID))+theme_light()+
  xlab("PCA 1")+ylab("PCA 2")+geom_point(size=3)+
  scale_colour_manual(values = colorBlindBlack8)+ggtitle("B (denovo)")


distanceDenovo<-dist(genlightDenovo,method = "euc")
clustersDenovo<-hclust(distanceDenovo,method = "ward.D")


dend_datadenovo <- dendro_data(clustersDenovo, type = "rectangle")
dend_datadenovo$labels$col<-as.character(sub('_[^_]*$', '', dend_datadenovo$labels$label))
dend_datadenovo$labels$pop<-factor(dend_datadenovo$labels$col,
                                   c("Cvi0","C24","Ler0","Gu0","Col0","Ei2"),c("Cvi-0","C24","Ler-0","Gu-0","Col-0","Ei-2"))
dend_datadenovo2<-dend_datadenovo$labels

colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pr2<-ggplot(dend_datadenovo$segments)+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+theme_light()+xlab("")+ylab("genetic distance")+
  geom_text(data = dend_datadenovo2, aes(x, y,col=pop,label=label),size =4,angle=90,nudge_y = -0.05)+
  geom_point(data = dend_datadenovo2,aes(x,y,col=pop),size=1.5)+
  scale_colour_manual(values = colorBlindBlack8)
pr2



vcf.fnref<-"results/snp_calling/reference//testRef.vcf.gz"
vcfref<-read.vcfR(vcf.fnref)
genlightref<-vcfR2genlight(vcfref)


pcaref<-glPca(genlightref,nf = 10)



pcaPlotref <- data.frame(pcaref$scores,sample=row.names(pcaref$scores),
                            pop=sub("_.*$","",row.names(pcaref$scores)),ID=row.names(pcaref$scores))

ggplot(pcaPlotref,aes(x=PC1,y=PC2,col=pop,label=ID))+theme_light()+
  xlab("PCA 1")+ylab("PCA 2")+geom_point(size=3)+
  scale_colour_manual(values = colorBlindBlack8)+ggtitle("A (reference)")

distanceref<-dist(genlightref,method = "euc")
clustersref<-hclust(distanceref,method = "ward.D")


dend_dataref <- dendro_data(clustersref, type = "rectangle")
dend_dataref$labels$col<-as.character(sub('_[^_]*$', '', dend_dataref$labels$label))
dend_dataref$labels$pop<-factor(dend_dataref$labels$col,
                                   c("Cvi0","C24","Ler0","Gu0","Col0","Ei2"),c("Cvi-0","C24","Ler-0","Gu-0","Col-0","Ei-2"))
dend_dataref2<-dend_dataref$labels

pr1<-ggplot(dend_dataref$segments)+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+theme_light()+xlab("")+ylab("genetic distance")+
  geom_point(data = dend_dataref2,aes(x,y,col=pop),size=1.5)+
  scale_colour_manual("accession",values = colorBlindBlack8)

pr2<-ggplot(dend_datadenovo$segments)+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+theme_light()+xlab("")+ylab("genetic distance")+
  geom_point(data = dend_datadenovo2,aes(x,y,col=pop),size=1.5)+
  scale_colour_manual(values = colorBlindBlack8)


g <- ggarrange(pr1,pr2, ncol=2, nrow=1, widths=c(6,6), align="h", labels="AUTO", common.legend = TRUE, legend="right", legend.grob=get_legend(pr1))
g

ggsave("SNP_Tree.tiff",g,height=4,width=9,dpi="retina")
