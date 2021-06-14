rm(list=ls())
library(ggpmisc)
library(ggplot2)

combineFilesRef<-function(bisFiles,indNames){  #Read and combine individual bismark files
  for(i in 1:length(bisFiles)){
    bisFiles[[i]]$indName<-indNames[i]
  }
  out<-do.call(rbind,bisFiles)
  colnames(out)<-c("chr","pos","pos2","fracMeth","nMeth","nNonMeth","indNames")
  out$loci<-paste(out$chr,out$pos,sep="_")
  out$total<-out$nMeth+out$nNonMeth
  out$fracMeth<-out$fracMeth/100
  return(out)  
}


####################
#Works the same as the denovo branch with out the lift over data.
ascencsions<-c("Col","C24","Cvi","Ei","Gu")
refAll<-data.frame()
for(i in ascencsions){
ascencsionID<-i
#load truth dataset and filter for 10x coverage
Truth<-read.table(grep(ascencsionID,list.files("/mnt/nfs/bioinfdata/home/NIOO/maartenp/finalBenchmarking/data/methylation-calls/",full.names = T),value=T),h=T) 
Truth<-Truth[Truth$total_bases>9,]
Truth$lociRef<-paste(Truth$chrom,Truth$pos,sep="_")

#Load the individual methylation data from epiGBS
allIndsBismarkRef<-list.files("/mnt/nfs/bioinfdata/home/NIOO/maartenp/epiGBS-ref/output/methylation_calling/",full.names = T)
bisFilesRef<-lapply(grep(paste(ascencsionID,".*cov.gz$",sep=""),allIndsBismarkRef,value = T),read.table)
indNamesRef<-sub("^.*//","",sub("_bismark.*$","",grep(ascencsionID,grep("cov.gz$",allIndsBismarkRef,value = T),value=T)))
bismarkRef<-combineFilesRef(bisFilesRef,indNamesRef)
bismarkRef$lociRef<-bismarkRef$loci

#Pool the data
pooledRef<-aggregate(.~lociRef,bismarkRef[,c("nMeth","total","lociRef")],sum)
pooledRef$chr<-as.numeric(sub("_.*$","",pooledRef$lociRef))
pooledRef$pos<-as.numeric(sub("^.*_","",pooledRef$lociRef))
pooledRef$lociRef<-paste(pooledRef$chr,pooledRef$pos-1,sep="_")

#Combine truth and pooled data
refVsPooledTruth<-merge(pooledRef,Truth,by='lociRef')
#Calculate methylation percentages
refVsPooledTruth$fracRef<-refVsPooledTruth$nMeth/refVsPooledTruth$total
refVsPooledTruth$fracTrue<-refVsPooledTruth$methylated_bases/refVsPooledTruth$total_bases
refVsPooledTruth10x<-refVsPooledTruth[refVsPooledTruth$total>length(bisFilesRef)*10,]
refVsPooledTruth10x$context<-"CHH" #Create context based on the truth dataset
refVsPooledTruth10x$context[grep("CG",refVsPooledTruth10x$mc_class)]<-"CG"
refVsPooledTruth10x$context[grep("C.G",refVsPooledTruth10x$mc_class)]<-"CHG"

refSave<-data.frame(fracPooled=refVsPooledTruth10x$fracRef,fracTrue=refVsPooledTruth10x$fracTrue,context=refVsPooledTruth10x$context,ascencsion=ascencsionID)
refAll<-rbind(refAll,refSave)
}


head(refSave)
colnames(refAll)<-c("fracPooled","fracTrue","context","ascencsion")
#Make the plot
ggplot(refAll,aes(x=fracPooled,y=fracTrue))+geom_point(alpha=0.05)+xlab("Fraction methylation Reference branch")+ylab("Fraction methylation truth")+
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..rr.label..,sep="~",parse=T)),col="blue", 
               parse = TRUE)+ggtitle("epiGBS reference branch vs Truth")+
  facet_grid(ascencsion~context)



# refVsTruth<-merge(bismarkRef,Truth,by='lociRef')

# refVsTruth$fracRef<-refVsTruth$nMeth/refVsTruth$total
# refVsTruth$fracTrue<-refVsTruth$methylated_bases/refVsTruth$total_bases
# refVsTruth10x<-refVsTruth[refVsTruth$total>9,]
# refVsTruth10x$context<-"CHH" #Create context based on the truth dataset
# refVsTruth10x$context[grep("CG",refVsTruth10x$mc_class)]<-"CG"
# refVsTruth10x$context[grep("C.G",refVsTruth10x$mc_class)]<-"CHG"


# ggplot(refVsTruth10x,aes(x=fracRef,y=fracTrue))+geom_point(alpha=0.05)+xlab("Fraction methylation Reference branch")+ylab("Fraction methylation truth")+
#   stat_poly_eq(formula = y~x,
#                aes(label = paste(..rr.label..,sep="~",parse=T)),col="blue",
#                 parse = TRUE)+ggtitle("epiGBS reference branch vs Truth")+
#   facet_grid(.~context)

# ggplot(refVsTruth10x,aes(x=fracRef,y=fracTrue))+geom_point(alpha=0.05)+xlab("Fraction methylation Reference branch")+ylab("Fraction methylation truth")+
#   stat_poly_eq(formula = y~x,
#                aes(label = paste(..rr.label..,sep="~",parse=T)),col="blue",
#                parse = TRUE)+ggtitle("epiGBS reference branch vs Truth")+
#   facet_grid(indNames~context)
