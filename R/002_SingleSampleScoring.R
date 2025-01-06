################################################################################
## Personalized molecular profile scoring benchmark
## R version 4.3.1 (2024-04-10)
## 
################################################################################
## Calculate single sample scores for LN and PRECISESADS cohorts:


# devtools::install_github("jordimartorell/pathMED")
##······································································· Step 0
## Set environment

setwd("D:/Work/BMS/PMPS/PMSPbenchmark")
set.seed(1234)

source(paste0(getwd(),"/code/utils.R"))

check.packages(c("vctrs","stats","stringr","GSVA","singscore","BiocParallel",
                 "car","pbapply","parallel","limma","decoupleR","reshape2",
                 "pathMED","doMC"))

## Load Datasets
load("D:/Work/BMS/PMPS/PMSPbenchmark/RData/PRECISESADS.RData")
load("D:/Work/BMS/PMPS/PMSPbenchmark/RData/LNDatasets.RData")


##······································································· Step 1
## Preparing data for scoring

methods<-c("GSVA","M-score","ssGSEA","singscore","Z-score","Plage",
           "AUCell","MDT","MLM","ORA","UDT","ULM",
           "FGSEA","norm_FGSEA","WMEAN","norm_WMEAN","corr_WMEAN",
           "WSUM","norm_WSUM","corr_WSUM")


## list of gene-expression from different studies
ref.data<-list("ln1"=list("expData"=cbind(HC.ln1,SLE.ln1),
                          "labels"=c(rep("Healthy",ncol(HC.ln1)),rep("SLE",ncol(SLE.ln1)))),
               "ln2"=list("expData"=cbind(HC.ln2,SLE.ln2),
                          "labels"=c(rep("Healthy",ncol(HC.ln2)),rep("SLE",ncol(SLE.ln2)))),
               "ln3"=list("expData"=cbind(HC.ln3,SLE.ln3),
                          "labels"=c(rep("Healthy",ncol(HC.ln3)),rep("SLE",ncol(SLE.ln3)))),
               "ln4"=list("expData"=cbind(HC.ln4,SLE.ln4),
                          "labels"=c(rep("Healthy",ncol(HC.ln4)),rep("SLE",ncol(SLE.ln4)))),
               "prec"=list("expData"=cbind(HC.prec,SLE.prec),
                           "labels"=c(rep("Healthy",ncol(HC.prec)),rep("SLE",ncol(SLE.prec)))) )

## Using SLEDiseasome 
load("D:/Work/BMS/PMPS/PMSPbenchmark/RData/SLEDiseasome.RData")


##······································································· Step 2
## Get Scores

SCORES<-lapply(ref.data,function(set){
  
  res<-lapply(methods,function(method){
    cat(paste0("\nRuning ",method,"\n"))
    dat<-set$expData
    labs<-set$labels
    
    if(method=="M-score"){ ## Get Mscores also for HC
      hc<-dat[,labs=="Healthy"]
      colnames(hc)<-paste0(colnames(hc),"_HC")
      dat<-cbind(dat,hc)
      labs<-c(labs,rep("SLE",ncol(hc)))
    }
    
    ncores<-ifelse(method %in% c("M-scores", "MDT", "UDT", "ULM"), 12, 1)
    
    tmp.scores<-getScores(inputData=dat, geneSets=geneset.DB, method = method,
                          labels = labs,cores = ncores, minSize=3,
                          times =ifelse(grepl("norm",method),100,
                                        ifelse(grepl("corr",method),200,2)))
    
    return(tmp.scores)
  })
  
  names(res)<-methods
  return(res)
})
names(SCORES)<-names(ref.data)

## Remove all not needed objects and save
save.image("D:/Work/BMS/PMPS/PMSPbenchmark/RData/SCORES.RData")

