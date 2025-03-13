################################################################################
## Personalized molecular profile scoring benchmark
## R version 4.3.1 (2024-04-10)
## 
################################################################################
#' Get ML-based models results


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

## Fix colnames
#colnames(SCORES$ln1$`M-Scores`)<-gsub("_HC","",colnames(SCORES$ln1$`M-Scores`))
#colnames(SCORES$ln2$`M-Scores`)<-gsub("_HC","",colnames(SCORES$ln2$`M-Scores`))
#colnames(SCORES$ln3$`M-Scores`)<-gsub("_HC","",colnames(SCORES$ln3$`M-Scores`))
#colnames(SCORES$ln4$`M-Scores`)<-gsub("_HC","",colnames(SCORES$ln4$`M-Scores`))


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

# modelList<-lapply(train,function(tmpTrain){
#   set.seed(1234)
#   tmpFit<-pathMED::trainModel(inputData =tmpTrain,
#                               metadata=CLIN.train,
#                               models=methodsML(outcomeClass="character",
#                                                algorithms="all"),
#                               var2predict="group",
#                               positiveClass="pLN",
#                               pairingColumn ="patientID",
#                               Koutter=4,
#                               Kinner=4,
#                               repeatsCV=2,
#                               continue_on_fail = TRUE,
#                               saveLogFile = NULL)
#   
#   print(tmpFit$stats)
#   gc()
#   return(tmpFit)
# })
# 
# 
# save.image((paste0(getwd(),"/RData/modelsLN.RData")))

library(pathMED)
library(tidyverse)
library(parallel)

# load("forML.RData")

progress_bar <- function(current, total, width = 50) {
  progress <- round(current / total * width)
  bar <- paste0(
    "[",
    paste0(rep("=", progress), collapse = ""),
    paste0(rep(" ", width - progress), collapse = ""),
    "] ",
    paste0(current, "/", total, " "),
    sprintf("%d%%", round(current / total * 100))
  )
  message("\r", bar)
  if (current == total) cat("\n")
  flush.console()
}

start_timer <- function(){
  start_time <- proc.time()
  return(start_time)
}

stop_timer <- function(start_time) {
  end_time <- proc.time()
  elapsed_time <- end_time - start_time
  elapsed_time <- elapsed_time[3]
  names(elapsed_time) <- NULL
  return(elapsed_time)
}

# Format the elapsed time nicely
time2seconds <- function(elapsed_time) {
  seconds <- as.numeric(elapsed_time, units = "secs")
  secs <- round(seconds, 2)  # Get remaining seconds and round
  return(secs)
}

seconds2mins <- function(secs, r_secs){
  mins <- secs - r_secs
  mins <- mins / 60
}

mins2hours <- function(mins, r_mins){
  hours <- mins - r_mins
  hours <- hours / 60
}

hours2days <- function(hours, r_hours){
  days <- hours - r_hours
  days <- days / 24
}

seconds2timer <- function(secs){
  if (secs < 60){
    time_string <- paste0(secs, " secs")
  } else{
    r_secs <- round(secs %% 60,2)
    time_string <- paste0(r_secs, " secs")
    mins <- seconds2mins(secs, r_secs)
    if (mins < 60){
      time_string <- paste0(mins, " mins ", time_string)
    } else{
      r_mins <- round(mins %% 60, 2)
      time_string <- paste0(r_mins, " mins ", time_string)
      hours <- mins2hours(mins, r_mins)
      if (hours < 24){
        time_string <- paste0(hours, " hours ", time_string)
      } else{
        r_hours <- round(hours %% 24, 2)
        time_string <- paste0(r_hours, " hours ", time_string)
        days <- hours2days(hours, r_hours)
        time_string <- paste0(days, " days ", time_string)
      }
    }
  }
  return(time_string)
}

message_time <- function(start_time, total){
  elapsed_time <- stop_timer(start_time)  # Stop the timer
  formatted_time <- time2seconds(elapsed_time)
  avg_time <- formatted_time / total
  message(paste("Total time:", seconds2timer(formatted_time)))
  message(paste("Average time per iteration:", seconds2timer(avg_time)))
}

mclapply_pb <- function(X, FUN, ..., mc.cores = 2, width = 50) {
  fun_name <- substitute(FUN)
  message(paste0("Running ", fun_name, " with ", mc.cores, " cores"))
  total <- length(X)
  batch_size <- mc.cores
  batches <- split(seq_along(X), ceiling(seq_along(X) / batch_size))
  total_batches <- length(batches)
  results <- vector("list", total)
  start_time <- start_timer()
  progress_bar(0, total, width)
  for (batch_indices in batches) {
    results[batch_indices] <- mclapply(X[batch_indices], FUN, ..., mc.cores = mc.cores)
    progress_bar(batch_indices[length(batch_indices)], total, width)
  }
  message_time(start_time, total)
  return(results)
}

trainModelSeed <- function(x, metadata, models, var2predict, positiveClass, pairingColumn, Koutter, Kinner, repeatsCV, continue_on_fail, saveLogFile, seed) {
  set.seed(seed)
  trainModel(
    x,
    metadata = metadata,
    models = models,
    var2predict = var2predict,
    positiveClass = positiveClass,
    pairingColumn = pairingColumn,
    Koutter = Koutter,
    Kinner = Kinner,
    repeatsCV = repeatsCV,
    continue_on_fail = continue_on_fail,
    saveLogFile = saveLogFile
  )
}

seeds <- c(1234, 12345, 123456, 1234567, 12345678)
modelLists <- lapply(seeds, function(seed) {
  modelList <- mclapply_pb(train, trainModelSeed, metadata = CLIN.train, 
                           models = methodsML(outcomeClass = "character", algorithms = "all"), 
                           var2predict = "group", positiveClass = "pLN", 
                           pairingColumn = "patientID", Koutter = 4, Kinner = 4, 
                           repeatsCV = 2, continue_on_fail = TRUE, saveLogFile = NULL,
                           seed = seed, mc.cores = 20)
  names(modelList) <- names(train)
  return(modelList)
})

names(modelLists) <- as.character(seeds)
saveRDS(modelLists, "C:/Users/danie/Desktop/WORK/BENCHMARK_25/RData/modelListsBenchmark.rds")



##······································································· Step 4
## Validation in external datasets


setwd("C:/Users/danie/Desktop/WORK/BENCHMARK_25")
set.seed(1234)

load(paste0(getwd(),"/RData/forML.RData"))

source(paste0(getwd(),"/code/utils.R"))

check.packages(c("vctrs","stats","stringr","caret","singscore","GSVA","parallel",
                 "lsa","metrica","BiocParallel","pbapply","reshape","ggplot2",
                 "ggpubr","psych","jaccard","decoupleR","reshape2","qvalue",
                 "ggbreak","pheatmap","scales","dplyr","DExMA"))

# devtools::install_local(paste0(getwd(),"/pathMED-main"))
library("pathMED")



modelList<-readRDS(paste0(getwd(),"/RData/modelListsBenchmark_5Seeds.rds"))


methods<-names(modelList[[1]])


TEST.tables<-lapply(methods,function(meth){
  
  tmpModels <- lapply(modelList, function(listML.seed) {
    listML.seed[names(listML.seed) == meth]
  })
  
  df<-do.call("cbind",lapply(1:length(tmpModels),function(iter){
    
    res<-do.call("cbind",lapply(unique(CLIN.test$dataset),function(dt){
      tmp.clin.test<-CLIN.test[CLIN.test$dataset==dt,]
      orig.data<-tmp.clin.test$group
      tmp.test<-test[[meth]][,rownames(tmp.clin.test)]
      
      tmp.test<-predictExternal(testData = tmp.test,
                                model = tmpModels[[iter]][[1]]$model,
                                realValues = orig.data,
                                positiveClass = "pLN")$stats
      x<-tmp.test$Score
      names(x)<-tmp.test$Metric
      return(x)
    }))
    
    res[is.na(res)]<-0
    
    #res<-rowMeans(res)
    return(res)
  }))
  colnames(df)<-paste0(meth,rep(c("_test1", "_test2"), times = 5), ".iter", rep(1:length(tmpModels), each = 2))
  return(df)
  
})

names(TEST.tables)<-methods



## Get results

table.result <- as.data.frame(sapply(TEST.tables, function(df) {
  # rowMeans_df <- rowMeans(df, na.rm = TRUE)  # Median
  # sd_df <- apply(df, 1, sd, na.rm = TRUE)    # SD
  # return(paste0(round(rowMeans_df, 2), " (±", round(sd_df, 2),")"))
  
  rowMedians_df <- apply(df,1,function(x){median(x,na.rm = TRUE)})
  sd_df <- apply(df, 1, sd, na.rm = TRUE)    # SD
  IQR<-apply(df,1,function(x){
    res<-quantile(x,c(0.75,0.25))
    return(res[1]-res[2])
  })
  #print(rowMedians_df+IQR)
  
  return(paste0(round(rowMedians_df, 2), " (±", round(sd_df, 2),")"))
}))
rownames(table.result)<-rownames(TEST.tables[[1]])


save.image((paste0(getwd(),"/RData/modelsLN_val.RData")))



## Save results in table format
table.result<-table.result[c("mcc","balacc","accuracy","precision","recall","specificity","fscore","npv"),]

rownames(table.result)<-c("MCC","BAcc","Acc","Prec","Recall","Spec","FScore","NPV")


means <- sapply(TEST.tables, function(df) {
  rowMedians_df <- apply(df,1,function(x){median(x,na.rm = TRUE)})
})[c("balacc","mcc"),]
order(apply(means,2,mean),decreasing = T)


table.result<-table.result[,order(as.numeric(means["mcc",]),decreasing = T)]


write.table(t(table.result),paste0(getwd(),"/Results/ModelValidation.txt"),sep="\t",quote = F)


