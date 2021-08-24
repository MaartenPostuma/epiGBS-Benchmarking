rm(list=ls())
library(reshape2)
library(ggplot2)
##Is now included in the paramTest mode.
setwd("../epiGBS/finalbenchmarkingepigbs/")
parameterTest<-read.csv('results/parameterTest/parameterTest.csv')
parameterTest$clustersMappedOnRef<-parameterTest$clustersMappedOnRef/2
parameterTest$nClusters<-parameterTest$nClusters/2

colnames(parameterTest)[1]<-"min.depth"

meltedData<-melt(parameterTest,id.vars = c("min.depth","clusterID"))
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


meltedData$variable<-factor(meltedData$variable,labels = c("N clusters","Mapping % denovo on ref","No. Clusters map on Ref","Mapping % reads on denovo","% multimapped","Average depth","CG methylation","CHG methylation","CHH methylation"))

Figure5<-meltedData[meltedData$variable%in%c("N clusters","Mapping % denovo on ref","Average depth"),]
Figure5$letter<-rep(c("A","B","C"),each=9)
ggplot(Figure5,aes(x=min.depth,y=value,col=as.factor(clusterID)))+geom_line(size=0.5)+facet_wrap(.~variable,scales = "free")+theme_light()+xlab("minimal depth")+geom_point(size=2)+
  scale_colour_manual("Clustering ID %",values = colorBlindBlack8)+    
  geom_text(aes(x = -Inf, y = Inf, label = letter, group = letter),size = 5,hjust = -1,vjust = 1.4,
            inherit.aes = FALSE)
ggsave("../epiGBS_paper/Figures/Figure5.pdf",height=9,width=9,dpi="retina")


FigureS4<-meltedData[-which(meltedData$variable%in%c("N clusters","Mapping % denovo on ref","Average depth")),]
FigureS4$letter<-rep(c("A","B","C","D","E","F"),each=9)
ggplot(FigureS4,aes(x=min.depth,y=value,col=as.factor(clusterID)))+geom_line(size=0.5)+facet_wrap(.~variable,scales = "free")+theme_light()+xlab("minimal depth")+geom_point(size=2)+
  scale_colour_manual("Clustering ID %",values = colorBlindBlack8)+    
  geom_text(aes(x = -Inf, y = Inf, label = letter, group = letter),size = 5,hjust = -1,vjust = 1.4,
            inherit.aes = FALSE)

ggsave("../epiGBS_paper/Figures/Supplemental/FigureS4.pdf",height=9,width=9,dpi="retina")
