################################################################################
## Personalized molecular profile scoring benchmark
## R version 4.3.1 (2024-04-10)
## 
################################################################################
#' Scores and Clinical associations


# devtools::install_github("jordimartorell/pathMED")
##······································································· Step 0
## Set environment

setwd("C:/Users/danie/Desktop/WORK/BENCHMARK_25")
set.seed(1234)

load(paste0(getwd(),"/RData/SCORES.RData"))

source(paste0(getwd(),"/code/utils.R"))

check.packages(c("vctrs","stats","stringr","caret","singscore","GSVA","parallel",
                 "lsa","metrica","BiocParallel","pbapply","reshape","ggplot2",
                 "ggpubr","psych","jaccard","decoupleR","reshape2","qvalue",
                 "ggbreak","pheatmap","scales","dplyr","DExMA"))

## Load pathMED locally
# devtools::install_local(paste0(getwd(),"/pathMED-main"))
library("pathMED")

colors <- c("GSVA"="#F4A460",
            "M-Scores"="#FF6347",
            "Z-score"="#A0522D",
            "ssGSEA"="burlywood4",
            "singscore"="rosybrown3",
            "Plage"="#20B2AA",
            "AUCell"="#A8A8A8",
            "MDT"="#4169E1",       
            "MLM"="#6495ED",
            "ORA"="#AFEECC",
            "UDT"="lightskyblue3",
            "ULM"="#1E90CC",
            "FGSEA"="#f0c35e",
            "norm_FGSEA"="gold3",
            "WMEAN"="#32CD32",
            "norm_WMEAN"="#6B8E23",
            "corr_WMEAN"= "#228B22",
            "WSUM"="orchid4",
            "norm_WSUM"="plum",
            "corr_WSUM"="#FFB6C1")

# methods<-methods[!grepl("corr",methods)]

clin.ln1$group<-ifelse(clin.ln1$NeprBX=="Proliferative","pLN","NoLN")
clin.ln2$group<-ifelse(clin.ln2$LN=="YES","pLN","NoLN")
clin.ln3$group<-ifelse(clin.ln3$pLN=="YES","pLN","NoLN")
clin.ln4$group<-ifelse(clin.ln4$pLN=="YES","pLN","NoLN")

colnames(SCORES$prec$`M-Scores`)<-gsub("_HC","",colnames(SCORES$prec$`M-Scores`))

##······································································· Step 1
## Q3: Are the scores capturing the disease clinical heterogeneity

PREC<-SCORES$prec[methods]

remPaths<-unique(c(unlist(lapply(PREC,function(dat){
  res<-rownames(dat)[(apply(dat, 1, function(row) any(is.na(row) | is.infinite(row)))) | 
                       (apply(dat, 1, function(row) all(row == 0, na.rm = TRUE)))]
  return(res)
}))))

PREC<-lapply(PREC,function(dat){
  dat<-dat[!rownames(dat) %in% remPaths,]
  return(dat)
})


## Make score-based clustering
RESULTS.CL<-lapply(PREC,function(tmp.dat){
  d = sweep(tmp.dat,1, apply(tmp.dat,1,median,na.rm=T))
  clusters = ConsensusClusterPlus(as.matrix(d),maxK=10,reps=250,pItem=0.8,verbose=T,
                                  pFeature=1,clusterAlg="km",distance="euclidean",
                                  innerLinkage = "complete",seed=12345,plot=NULL)
  clusters<-do.call("cbind",lapply(2:10,function(k){
    return(clusters[[k]]$consensusClass)
  }))
  colnames(clusters)<-paste0("k",2:10)
  return(clusters)
})
names(RESULTS.CL)<-names(PREC)




## Get enrichment with clinical variables and clusters
clin.hc<-as.data.frame(matrix(data=NA,ncol=ncol(clin.prec$symptoms),nrow=ncol(HC.prec)))
colnames(clin.hc)<-colnames(clin.prec$symptoms)
rownames(clin.hc)<-colnames(HC.prec)
clin.hc$Diagnosis_Arm<-"Healthy"
clin.prec$symptoms<-rbind(clin.prec$symptoms,clin.hc)


CLIN.table<-lapply(clin.prec,function(clin){
  
  res.met<-do.call("cbind",lapply(RESULTS.CL,function(met){
    pats<-intersect(rownames(met),rownames(clin))
    
    
    allK<-do.call("cbind",lapply(c("k3","k4","k5","k6","k7"),function(k){
      
      metadata<-cbind(met[pats,k],clin[pats,])
      colnames(metadata)[1]<-"group"
      metadata$group<-paste0("CL",metadata$group)
      
      res<-testAssoc(metadata=metadata,p.adj = "none", remove.vars=NULL,colors=NULL,
                     normality=NULL,positiveClass=NULL)
      
      res.k<-unlist(lapply(1:length(res),function(i){
        x<-res[[i]]$general.test
        switch(class(x),
               NULL = { 
                 return(NA) 
               },
               numeric = { 
                 return(x) 
               },
               data.frame = {
                 return(as.numeric(x$pvalue))
               })
      }))
      names(res.k)<-names(res)
      return(res.k)
    }))
    colnames(allK)<-c("k3","k4","k5","k6","k7")
    allK <- allK[, order(-colSums(allK <= 0.05))]
    
    return(allK[,1,drop=FALSE])
  })) #met
  colnames(res.met)<-names(RESULTS.CL)
  
  return(res.met)
  
})
names(CLIN.table)<-names(clin.prec)



## Plot Heatmap ··············
DAT<-NULL
ANN<-NULL

for(i in 1:3){
  
  dat<-CLIN.table[[i]]
  dat<-dat[!is.na(rownames(dat)),]
  dat[is.na(dat)]<-1
  dat<- (-log10(dat))
  dat<-dat[ifelse(apply(dat,1,max)<(-log10(0.05)),F,T),]
  
  DAT<-rbind(DAT,dat)
  ann<-data.frame("ids"=rownames(dat),"type"=names(CLIN.table)[i])
  ANN<-rbind(ANN,ann)

}
rownames(ANN)<-ANN$ids


breaks <- c(seq(0,-log10(0.05), length.out = 100),
            seq(-log10(0.05)+0.001,-log10(0.01), length.out = 100),
            seq(-log10(0.01)+0.001, 4, length.out = 100),
            seq(4.001, max(DAT), length.out = 100))


anncolors<-list("type"=c("autoantibodies"="#cd5554","symptoms"="#0198b1","cytokines"="#8fbc8f"))

colors <- c(colorRampPalette(c("white", "#FFCC99"))(100),
            colorRampPalette(c("#FFCC99", "coral3"))(100),
            colorRampPalette(c("coral3", "coral4"))(100),
            colorRampPalette(c("coral4", "#660033"))(100))

pheatmap(t(DAT),scale="none",show_colnames = T,cluster_cols = T,
         cluster_rows = T, show_rownames = T, border_color = "black", fontsize = 7,
         breaks=breaks,annotation_col = ANN[,"type",drop=FALSE],
         color = colors,annotation_colors = anncolors)


save.image(paste0(getwd(),"/RData/ClinicalAssocc.RData"))

