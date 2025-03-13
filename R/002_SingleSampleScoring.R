################################################################################
## Personalized molecular profile scoring benchmark
## R version 4.3.1 (2024-04-10)
## 
################################################################################
## Calculate single sample scores for LN and PRECISESADS cohorts:


# devtools::install_github("jordimartorell/pathMED")
##······································································· Step 0
## Set environment

setwd("C:/Users/danie/Desktop/WORK/BENCHMARK_25")
set.seed(1234)

source(paste0(getwd(),"/code/utils.R"))

check.packages(c("vctrs","stats","stringr","GSVA","singscore","BiocParallel",
                 "car","pbapply","parallel","limma","decoupleR","reshape2",
                 #"pathMED",
                 "doMC"))

## Load pathMED locally
# devtools::install_local(paste0(getwd(),"/pathMED-main"))
library("pathMED")

## Load Datasets
load(paste0(getwd(),"/RData/PRECISESADS.RData"))
load(paste0(getwd(),"/RData/LNDatasets.RData"))


##······································································· Step 1
## Preparing data for scoring

methods<-c("GSVA","M-Scores","ssGSEA","singscore","Z-score","Plage","AUCell",
           "MDT","MLM","ORA","UDT","ULM","FGSEA","norm_FGSEA","WMEAN",
           "norm_WMEAN","WSUM","norm_WSUM")


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
geneset.DB<-readRDS(paste0(getwd(),"/RData/SLEDiseasome_db.rds"))


##······································································· Step 2
## Get Scores


# res <- readRDS("C:/Users/danie/Desktop/WORK/BENCHMARK_25/RData/tmpRes.rds")
# SCORES <- readRDS("C:/Users/danie/Desktop/WORK/BENCHMARK_25/RData/tmpScores.rds")

SCORES<-list()
for(i in 1:length(ref.data)){
  set<-ref.data[[i]]
  
  res<-list()
  for(method in methods){
    
    cat(paste0("\nRuning ",method,"\n"))
    dat<-set$expData
    labs<-set$labels
    dat<-dat[!is.na(rownames(dat)) & !rownames(dat)=="",]
    
    ncores<-ifelse(method %in% c("M-Scores", "MDT", "UDT", "ULM"), 14, 1)
    
    tmp.scores<-getScores(inputData=dat, geneSets=geneset.DB, method = method,
                          labels = labs, returnHC = TRUE, cores = ncores, 
                           minSize=3, times =ifelse(grepl("norm",method),100,2))
    res[[method]]<-tmp.scores
    saveRDS(res,"C:/Users/danie/Desktop/WORK/BENCHMARK_25/RData/tmpRes.rds")
  }
  
  SCORES[[names(ref.data)[i]]]<-res
  saveRDS(SCORES,"C:/Users/danie/Desktop/WORK/BENCHMARK_25/RData/tmpScores.rds")
  
}

# SCORES<-lapply(ref.data,function(set){
#   
#   res<-lapply(methods,function(method){
#     cat(paste0("\nRuning ",method,"\n"))
#     dat<-set$expData
#     labs<-set$labels
#     
#     if(method=="M-score"){ ## Get Mscores also for HC
#       hc<-dat[,labs=="Healthy"]
#       colnames(hc)<-paste0(colnames(hc),"_HC")
#       dat<-cbind(dat,hc)
#       labs<-c(labs,rep("SLE",ncol(hc)))
#     }
#     
#     ncores<-ifelse(method %in% c("M-Scores", "MDT", "UDT", "ULM"), 12, 1)
#     
#     tmp.scores<-getScores(inputData=dat, geneSets=geneset.DB, method = method,
#                           labels = labs,cores = ncores, minSize=3,
#                           times =ifelse(grepl("norm",method),100,
#                                         ifelse(grepl("corr",method),200,2)))
#     
#     return(tmp.scores)
#   })
#   
#   names(res)<-methods
#   return(res)
# })
# names(SCORES)<-names(ref.data)

## Remove all not needed objects and save


save.image(paste0(getwd(),"/RData/SCORES.RData"))

