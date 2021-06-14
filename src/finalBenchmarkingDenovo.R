rm(list=ls())
library(ggpmisc)
library(ggplot2)

combineFilesDenovo<-function(bisFiles,indNames){  #Read and combine individual bismark files
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

ascencsions<-c("Col","C24","Cvi","Ei","Gu")
refAll<-data.frame()
liftOver<-read.table("/mnt/nfs/bioinfdata/home/NIOO/maartenp/finalBenchmarking/results/liftover.depth") #Load lift over coordinates
liftOver$lociDenovo<-paste(liftOver$V6,liftOver$V8,sep="_") #Make the same columns in the lift over file
liftOver$lociRef<-paste(liftOver$V1,liftOver$V2,sep="_") #To move the epiGBS contigs to the reference loci


for(i in ascencsions){ #Run for loop over all ascensions
  ascencsionID<-i
  #Load "Truth dataset", format loci to be similar only use sites with 10x coverage or more 
  Truth<-read.table(grep(ascencsionID,list.files("/mnt/nfs/bioinfdata/home/NIOO/maartenp/finalBenchmarking/data/methylation-calls/",full.names = T),value=T),h=T) 
  Truth$lociRef<-paste(Truth$chrom,Truth$pos,sep="_")
  Truth<-Truth[Truth$lociRef%in%liftOver$lociRef&Truth$total_bases>9,]
  

  #Load all methylation files, and merge with the lift over coordinates
  allIndsBismarkRef<-list.files("/mnt/nfs/bioinfdata/home/NIOO/maartenp/epiGBS-denovo/output/methylation_calling/",full.names = T)
  bisFilesRef<-lapply(grep(paste(ascencsionID,".*cov.gz$",sep=""),allIndsBismarkRef,value = T),read.table)
  indNamesRef<-sub("^.*//","",sub("_bismark.*$","",grep(ascencsionID,grep("cov.gz$",allIndsBismarkRef,value = T),value=T)))
  bismarkRef<-combineFilesRef(bisFilesRef,indNamesRef)
  bismarkRef$lociDenovo<-bismarkRef$loci
  bismarkLifted<-merge(bismarkRef,liftOver,by="lociDenovo")
  
  #Pool all individuals
  pooledRef<-aggregate(.~lociRef,bismarkLifted[,c("nMeth","total","lociRef")],sum)
  #Combine the pooled results with the "truth dataset", calculate methylation fractions and save them
  refVsPooledTruth<-merge(pooledRef,Truth,by='lociRef')
  refVsPooledTruth$fracRef<-refVsPooledTruth$nMeth/refVsPooledTruth$total
  refVsPooledTruth$fracTrue<-refVsPooledTruth$methylated_bases/refVsPooledTruth$total_bases
  refVsPooledTruth10x<-refVsPooledTruth[refVsPooledTruth$total>length(bisFilesRef)*10,]
  refVsPooledTruth10x$context<-"CHH" #Create context based on the truth dataset
  refVsPooledTruth10x$context[grep("CG",refVsPooledTruth10x$mc_class)]<-"CG"
  refVsPooledTruth10x$context[grep("C.G",refVsPooledTruth10x$mc_class)]<-"CHG"
  
  refSave<-data.frame(fracPooled=refVsPooledTruth10x$fracRef,fracTrue=refVsPooledTruth10x$fracTrue,context=refVsPooledTruth10x$context,ascencsion=ascencsionID)
  refAll<-rbind(refAll,refSave)
}

#Plot the figure.
colnames(refAll)<-c("fracPooled","fracTrue","context","accession")
ggplot(refAll,aes(x=fracPooled,y=fracTrue))+geom_point(alpha=0.05)+xlab("Fraction methylation Denovo branch")+ylab("Fraction methylation truth")+
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..rr.label..,sep="~",parse=T)),col="blue", 
               parse = TRUE)+ggtitle("epiGBS denovo branch vs Truth")+
  facet_grid(ascencsion~context)


# 
# refVsTruth<-merge(bismarkRef,Truth,by='lociRef')
# 
# refVsTruth$fracRef<-refVsTruth$nMeth/refVsTruth$total
# refVsTruth$fracTrue<-refVsTruth$methylated_bases/refVsTruth$total_bases
# refVsTruth10x<-refVsTruth[refVsTruth$total>9,]
# refVsTruth10x$context<-"CHH" #Create context based on the truth dataset
# refVsTruth10x$context[grep("CG",refVsTruth10x$mc_class)]<-"CG"
# refVsTruth10x$context[grep("C.G",refVsTruth10x$mc_class)]<-"CHG"
# 
# 
# ggplot(refVsTruth10x,aes(x=fracRef,y=fracTrue))+geom_point(alpha=0.05)+xlab("Fraction methylation Reference branch")+ylab("Fraction methylation truth")+
#   stat_poly_eq(formula = y~x,
#                aes(label = paste(..rr.label..,sep="~",parse=T)),col="blue",
#                parse = TRUE)+ggtitle("epiGBS reference branch vs Truth")+
#   facet_grid(.~context)
# 
# ggplot(refVsTruth10x,aes(x=fracRef,y=fracTrue))+geom_point(alpha=0.05)+xlab("Fraction methylation Reference branch")+ylab("Fraction methylation truth")+
#   stat_poly_eq(formula = y~x,
#                aes(label = paste(..rr.label..,sep="~",parse=T)),col="blue",
#                parse = TRUE)+ggtitle("epiGBS reference branch vs Truth")+
#   facet_grid(indNames~context)
