rm(list=ls())
library(ggpmisc)
library(ggplot2)
library(dplyr)
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

ascencsions<-c("Col","C24","Cvi","Ei","Gu")
refAll<-data.frame()
for(i in ascencsions){
ascencsionID<-i

Truth<-read.table(grep(ascencsionID,list.files("/mnt/nfs//home//maartenp/finalBenchmarking/data/methylation-calls/",full.names = T),value=T),h=T) 
Truth<-Truth[Truth$total_bases>9,]
Truth$lociRef<-paste(Truth$chrom,Truth$pos,sep="_")


allIndsBismarkRef<-list.files("/mnt/nfs//home//maartenp/adam/adam-ref/output/methylation_calling/",full.names = T)
bisFilesRef<-lapply(grep(paste(ascencsionID,".*cov.gz$",sep=""),allIndsBismarkRef,value = T),read.table)
indNamesRef<-sub("^.*//","",sub("_bismark.*$","",grep(ascencsionID,grep("cov.gz$",allIndsBismarkRef,value = T),value=T)))
bismarkRef<-combineFilesRef(bisFilesRef,indNamesRef)
bismarkRef$lociRef<-bismarkRef$loci
pooledRef<-aggregate(.~lociRef,bismarkRef[,c("nMeth","total","lociRef")],sum)
pooledRef$chr<-as.numeric(sub("_.*$","",pooledRef$lociRef))
pooledRef$pos<-as.numeric(sub("^.*_","",pooledRef$lociRef))
pooledRef$lociRef<-paste(pooledRef$chr,pooledRef$pos-1,sep="_")
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

colnames(refAll)<-c("fracPooled","fracTrue","context","accession")
A<-ggplot(refAll,aes(x=fracPooled,y=fracTrue))+geom_point(alpha=0.05)+xlab("fraction methylation reference branch")+ylab("fraction methylation baseline")+
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..rr.label..,sep="~",parse=T)),col="blue", 
               parse = TRUE)+#ggtitle("epiGBS reference branch vs baseline")+
  stat_smooth(method = "lm",formula=y~x,col="blue",se = F)+
  theme_minimal(base_size = 15)+  
  theme(strip.background =element_rect(fill="black"))+
  theme(strip.text = element_text(colour = 'white'))+
  facet_grid(accession~context)+
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(0,"line"))+
  geom_vline(data=tibble(f=2, x=c(-1,1)*Inf), aes(xintercept=x), col="grey30",size=0.5)+
  geom_hline(data=tibble(f=2, x=c(-1,1)*Inf), aes(yintercept=x), col="grey30",size=1)
  ggsave("results/Figures/referenceVsBaseline.pdf",height=10,width=6)

  
  rm(list=ls()[ls()!="A"])
  library(ggpmisc)
  library(ggplot2)
  library(dplyr)
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
  liftOver<-read.table("/mnt/nfs//home//maartenp/finalBenchmarking/results/liftover.depth")
  liftOver$lociDenovo<-paste(liftOver$V6,liftOver$V8,sep="_") #Make the same columns in the lift over file
  liftOver$lociRef<-paste(liftOver$V1,liftOver$V2,sep="_") #To move the epiGBS contigs to the reference loci
  
  
  for(i in ascencsions){
    ascencsionID<-i
    
    Truth<-read.table(grep(ascencsionID,list.files("/mnt/nfs//home//maartenp/finalBenchmarking/data/methylation-calls/",full.names = T),value=T),h=T) 
    Truth$lociRef<-paste(Truth$chrom,Truth$pos,sep="_")
    Truth<-Truth[Truth$lociRef%in%liftOver$lociRef&Truth$total_bases>9,]
    
    allIndsBismarkRef<-list.files("/mnt/nfs//home//maartenp/epiGBS-denovo/output/methylation_calling/",full.names = T)
    bisFilesRef<-lapply(grep(paste(ascencsionID,".*cov.gz$",sep=""),allIndsBismarkRef,value = T),read.table)
    indNamesRef<-sub("^.*//","",sub("_bismark.*$","",grep(ascencsionID,grep("cov.gz$",allIndsBismarkRef,value = T),value=T)))
    bismarkRef<-combineFilesDenovo(bisFilesRef,indNamesRef)
    bismarkRef$lociDenovo<-bismarkRef$loci
    bismarkLifted<-merge(bismarkRef,liftOver,by="lociDenovo")
    
    
    pooledRef<-aggregate(.~lociRef,bismarkLifted[,c("nMeth","total","lociRef")],sum)
    #pooledRef$chr<-as.numeric(sub("_.*$","",pooledRef$lociRef))
    #pooledRef$pos<-as.numeric(sub("^.*_","",pooledRef$lociRef))
    #pooledRef$lociRef<-paste(pooledRef$chr,pooledRef$pos-1,sep="_")
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
  
  
  colnames(refAll)<-c("fracPooled","fracTrue","context","ascencsion")
  B<-ggplot(refAll,aes(x=fracPooled,y=fracTrue))+geom_point(alpha=0.05)+xlab("fraction methylation denovo branch")+ylab("fraction methylation baseline")+
    stat_poly_eq(formula = y~x, 
                 aes(label = paste(..rr.label..,sep="~",parse=T)),col="blue", 
                 parse = TRUE)+#ggtitle("epiGBS reference branch vs baseline")+
    stat_smooth(method = "lm",formula=y~x,col="blue",se = F)+
    theme_minimal(base_size = 15)+  
    theme(strip.background =element_rect(fill="black"))+
    theme(strip.text = element_text(colour = 'white'))+
    facet_grid(ascencsion~context)+
    theme(panel.border = element_blank(),
          panel.background = element_blank(),
          panel.spacing.x = unit(0,"line"))+
    geom_vline(data=tibble(f=2, x=c(-1,1)*Inf), aes(xintercept=x), col="grey30",size=0.5)+
    geom_hline(data=tibble(f=2, x=c(-1,1)*Inf), aes(yintercept=x), col="grey30",size=1)
  ggsave("results/Figures/denovoVsBaseline.pdf",height=10,width=6)
