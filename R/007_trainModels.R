################################################################################
## Personalized molecular profile scoring benchmark
## R version 4.3.1 (2024-04-10)
## 
################################################################################
#' Get ML-based models results


# devtools::install_github("jordimartorell/pathMED")
##······································································· Step 0
## Set environment

setwd("D:/Work/BMS/PMPS/PMSPbenchmark")
set.seed(1234)

load(paste0(getwd(),"/RData/SCORES.RData"))

source(paste0(getwd(),"/code/utils.R"))

check.packages(c("vctrs","stats","stringr","caret","singscore","GSVA","parallel",
                 "lsa","metrica","BiocParallel","pbapply","reshape","ggplot2",
                 "ggpubr","psych","jaccard","decoupleR","reshape2","qvalue",
                 "ggbreak","pheatmap","pathMED","scales","dplyr","DExMA"))


colors <- c("GSVA"="#F4A460",
            "M-score"="#FF6347",
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

methods<-methods[!grepl("corr",methods)]

clin.ln1$group<-ifelse(clin.ln1$NeprBX=="Proliferative","pLN","NoLN")
clin.ln2$group<-ifelse(clin.ln2$LN=="YES","pLN","NoLN")
clin.ln3$group<-ifelse(clin.ln3$pLN=="YES","pLN","NoLN")
clin.ln4$group<-ifelse(clin.ln4$pLN=="YES","pLN","NoLN")

## Fix colnames
colnames(SCORES$ln1$`M-score`)<-gsub("_HC","",colnames(SCORES$ln1$`M-score`))
colnames(SCORES$ln2$`M-score`)<-gsub("_HC","",colnames(SCORES$ln2$`M-score`))
colnames(SCORES$ln3$`M-score`)<-gsub("_HC","",colnames(SCORES$ln3$`M-score`))
colnames(SCORES$ln4$`M-score`)<-gsub("_HC","",colnames(SCORES$ln4$`M-score`))


##······································································· Step 1
## Prepare Train and Test sets for Run getML (gene expression)

## Same genes are needed for ML
features<-Reduce(intersect, list(rownames(SLE.ln1), rownames(SLE.ln2), 
                                 rownames(SLE.ln3), rownames(SLE.ln4)))
SLE.ln1<-SLE.ln1[features,]
HC.ln1<-HC.ln1[features,]
SLE.ln2<-SLE.ln2[features,]
HC.ln2<-HC.ln2[features,]
SLE.ln3<-SLE.ln3[features,]
HC.ln3<-HC.ln3[features,]
SLE.ln4<-SLE.ln4[features,]
HC.ln4<-HC.ln4[features,]
## Train: ln1 + ln3, Test: ln2 + ln4

## Feature selection based on Meta-analysis using training datasets
obj.MA<-list("ln1"=list("exp"=as.matrix(cbind(HC.ln1,SLE.ln1)),
                        "cls"=c(rep(0,ncol(HC.ln1)),rep(1,ncol(SLE.ln1)))),
             "ln3"=list("exp"=as.matrix(cbind(HC.ln3,SLE.ln3)),
                        "cls"=c(rep(0,ncol(HC.ln3)),rep(1,ncol(SLE.ln3)))) )


## Test heterogeneity 
heterogeneityTest(obj.MA)

resultsMA <- metaAnalysisDE(obj.MA, typeMethod="REM",
                            missAllow=0.3, proportionData=1)

selGenes<-rownames(resultsMA[resultsMA$FDR<=0.05 & abs(resultsMA$AveFC)>1,])
print(length(selGenes))

## Select samples for train and test 
train<-list()
test<-list()

CLIN.train<-data.frame("group"=c(clin.ln1$group,clin.ln3$group),
                       "dataset"=c(rep("ln1",ncol(SLE.ln1)),rep("ln3",ncol(SLE.ln3))),
                       "patientID"=c(clin.ln1$Patient,clin.ln3$PatientID))
rownames(CLIN.train)<-c(rownames(clin.ln1),rownames(clin.ln3))

CLIN.test<-data.frame("group"=c(clin.ln2$group,clin.ln4$group),
                      "dataset"=c(rep("ln2",ncol(SLE.ln2)),rep("ln4",ncol(SLE.ln4))))
rownames(CLIN.test)<-c(rownames(clin.ln2),rownames(clin.ln4))

## Gene expression data
train[["geneExp"]]<-cbind(SLE.ln1[selGenes,rownames(clin.ln1)],
                          SLE.ln3[selGenes,rownames(clin.ln3)])
test[["geneExp"]]<-cbind(SLE.ln2[selGenes,rownames(clin.ln2)],
                         SLE.ln4[selGenes,rownames(clin.ln4)])

## Gene expression normalized by zscore by samples
train[["geneExp_norm"]] <- apply(train$geneExp, 2, function(x) (x - mean(x)) / sd(x))

test[["geneExp_norm"]] <- apply(test$geneExp, 2, function(x) (x - mean(x)) / sd(x))


##······································································· Step 2
## Prepare Train and Test sets for Run getML (scores)
#fc<-c(0.75,1.5,0.1,0.1,3.2,0.125,0.075,3,0.25,1,0,1.25,0.15,0.5,0.75,2,15,2)

## Select 150 most significant features (FDR + abs(fc))

for(m in 1:length(methods)){
  
  print(methods[m])
  
  dat<-list("ln1"=SCORES$ln1[[methods[m]]],
            "ln2"=SCORES$ln2[[methods[m]]],
            "ln3"=SCORES$ln3[[methods[m]]],
            "ln4"=SCORES$ln4[[methods[m]]])
  
  dat <- lapply(dat, function(mat){
    mat<-as.matrix(mat)
    return(mat[!apply(is.infinite(mat), 1, any), ])
  }) 
  sharedPaths <- Reduce(intersect, lapply(dat, rownames))

  obj.MA<-list("ln1"=list("exp"=as.matrix(cbind(dat$ln1[sharedPaths,colnames(HC.ln1)],
                                                dat$ln1[sharedPaths,colnames(SLE.ln1)])),
                          "cls"=c(rep(0,ncol(HC.ln1)),rep(1,ncol(SLE.ln1)))),
               "ln3"=list("exp"=as.matrix(cbind(dat$ln3[sharedPaths,colnames(HC.ln3)],
                                                dat$ln3[sharedPaths,colnames(SLE.ln3)])),
                          "cls"=c(rep(0,ncol(HC.ln3)),rep(1,ncol(SLE.ln3)))) )
  
  resultsMA <- metaAnalysisDE(obj.MA, typeMethod="REM",
                              missAllow=0.3, proportionData=1)
  
  resultsMA<-resultsMA[!is.na(resultsMA$AveFC),]
  
  resultsMA<-resultsMA[resultsMA$FDR<=0.05 & resultsMA$propDataset==1,]
  resultsMA$abs<-abs(resultsMA$AveFC)
  
  resultsMA<-resultsMA[order(resultsMA$abs,decreasing = T),]
  
  if(nrow(resultsMA)>150){ resultsMA<-resultsMA[1:150,] }
  
  #plot(density(resultsMA$AveFC),main = methods[m])
  selPaths<-rownames(resultsMA); print(length(selPaths))
  
  train[[methods[m]]]<-cbind(SCORES$ln1[[methods[m]]][selPaths,rownames(clin.ln1)],
                             SCORES$ln3[[methods[m]]][selPaths,rownames(clin.ln3)])
  
  test[[methods[m]]]<-cbind(SCORES$ln2[[methods[m]]][selPaths,rownames(clin.ln2)],
                            SCORES$ln4[[methods[m]]][selPaths,rownames(clin.ln4)])
}

rm(dat,m,obj.MA,resultsMA,features,selGenes,selPaths,sharedPaths)
gc()


##······································································· Step 3
## Run GetML for all training sets

modelList<-lapply(train,function(tmpTrain){
  tmpFit<-pathMED::trainModel(expData=tmpTrain,
                              metadata=CLIN.train,
                              models=methodsML(outcomeClass="character",
                                               algorithms="all"),
                              var2predict="group",
                              positiveClass="pLN",
                              paired ="patientID",
                              Koutter=4,
                              Kinner=4,
                              repeatsCV=2,
                              continue_on_fail = TRUE,
                              saveLogFile = NULL)
  print(tmpFit$stats)
  gc()
  return(tmpFit)
})


save.image((paste0(getwd(),"/RData/modelsLN.RData")))


##······································································· Step 4
## Validation in external datasets

## Select subset of each dataset, because test are imbalanced in classes
set.seed(1234)

baTest<-do.call("rbind",lapply(unique(CLIN.test$dataset),function(ds){
  tmp<-CLIN.test[CLIN.test$dataset==ds,]
  mg<-names(which.min(table(tmp$group)))
  res<-tmp[tmp$group==mg,]
  other<-tmp[tmp$group!=mg,]
  
  size<-round(nrow(res)+ (nrow(res)*0.25),digits = 0)
  if(size>nrow(other)){size=nrow(other)}
  other<-other[sample(1:nrow(other),size),]
  res<-rbind(res,other)
  return(res)
}))

##

orig.data<-baTest$group
names(orig.data)<-rownames(baTest)

df<-do.call("cbind",lapply(1:length(modelList),function(i){
  
  tmp.test<-predictExternal(testData = t(test[[i]][,rownames(baTest)]),
                            model = modelList[[i]]$model,
                            realValues = orig.data,
                            positiveClass = "pLN")$stats
  x<-tmp.test$Score
  names(x)<-tmp.test$Metric
  return(x)
}))
colnames(df)<-names(modelList)

## 

save.image((paste0(getwd(),"/RData/modelsLN_val.RData")))


df<-df[c("mcc","balacc","accuracy","precision","recall","specificity","fscore","npv"),
       order(df["mcc",],decreasing = T)]
df<-round(df, 3)

rownames(df)<-c("MCC","BAcc","Acc","Prec","Recall","Spec","FScore","NPV")

write.table(df,paste0(getwd(),"/Results/ModelValidation.txt"),sep="\t",quote = F)





  

  

