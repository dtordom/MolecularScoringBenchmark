#'··············································································
#' Benchmarking single-sample gene set scoring methods for application in 
#' precision medicine  
#' 
#' R version 4.5.1 (2025-06-13 ucrt)
#' contact: danieltorodominguez@gmail.com
#' 
#' Example of how to get single sample scores using pathMED
#'··············································································

# ······································································· Step 1
## Set environment

set.seed(12345)

## Install pathMED and dependencies
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pathMED")


## Load libraries
library("GEOquery")
library("parallel")
library("doParallel")
library("pathMED")
library("pbapply")
library("BiocParallel")
library("dplyr")
library("matrixStats")
library("caret")
library("stringr")
library("tidyr")

## Load functions
## Merge expression by gene. Internaly used by annotateGenes()
#' @param gene Gene identifier to annotate
#' @param genome Table with gene-probe set connections
#' @param expressionMatrix Gene expression data.frame
#' @param method Method used to merge expression of several probes pointing to 
#' the same gene. "median", "mean", "max", "sum"
mergeExpression = function(gene, 
                           genome,
                           expressionMatrix,
                           method){ 
  
  probe_ids = genome[genome$toGenes==gene,"fromGenes"] ## Select probes for each gene
  if (length(probe_ids)>1){
    method<-ifelse(method %in% c("median","mean","sum","max"),method,"median")
    switch (method,
            median = { res =  colMedians(as.matrix(expressionMatrix[probe_ids,]))},
            mean = { res =  colMeans2(as.matrix(expressionMatrix[probe_ids,]))},
            max = { res<-apply(as.matrix(expressionMatrix[probe_ids,]),2,function(x){
              max(x,na.rm = TRUE)})},
            sum = { res<-apply(as.matrix(expressionMatrix[probe_ids,]),2,function(x){
              sum(x,na.rm = TRUE)})})
  } else{
    res <- as.numeric(expressionMatrix[probe_ids,])  
  }
  return(res)
}

## Log Normalization
#' @param data Gene expression data.frame
norm.log<-function(data){
  qx <- as.numeric(quantile(data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) { data[which(data <= 0)] <- NaN
  data <- log2(data) }
  return(data)
}

## Summary of gene variance
#' @param data Gene expression data.frame
summ.var<-function(data){ ## Gene expression matrix
  require(caret)
  nzv <- nearZeroVar(t(data), saveMetrics= TRUE)
  return(nzv)
}

# ······································································· Step 2
## Load example data from public repository (NCBI GEO)

## Download. data from NCBI GEO
#' GSE99967: Dataset with 42 samples from Systemic Lupus Erythematosus and 17
#' Healthy samples
gset <- getGEO("GSE99967", GSEMatrix = TRUE)[[1]] 
data <- exprs(gset) ## Get gene-expression matrix
data <- norm.log(data) ## Log2 transformation

nonVar.genes <- summ.var(data = data) ## Remove genes with near to zero variance
data <- data[!nonVar.genes$nzv, ]

## Annotate Probes to GeneSymbol
## Download GPL21970 db from NCBI GEO
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL21970
genome <- getGEO("GPL21970")@dataTable@table[, c("GeneSymbol", "ID")]
colnames(genome) <- c("toGenes", "fromGenes")
genome <- genome[genome$fromGenes != "", ] ## Filtering non anotategd genes
genome <- genome[genome$fromGenes %in% rownames(data), ]
genome <- na.omit(genome)
data <- data[genome$fromGenes, ]
finalGenes <- unique(genome$toGenes)

temp <- mclapply(finalGenes, mergeExpression,
                 genome = genome,
                 expressionMatrix = data,
                 method ="median", mc.cores = 1)
temp <- as.data.frame(do.call("rbind", temp))
rownames(temp) <- finalGenes
colnames(temp) <- colnames(data)
data.exp <- temp

rm(temp,data,genome,finalGenes)

## Preparing metadata
metadata <- phenoData(gset) ## Download clinical data from NCBI GEO
metadata <- pData(metadata)[, c("title", "disease state:ch1")]
colnames(metadata) <- c("title", "diagnosis")
metadata$diagnosis<-ifelse(metadata$diagnosis=="SLE","SLE","Healthy")
data.exp<-data.exp[,rownames(metadata)]

## Remove non-variable genes in Healthy or SLE samples
keepGenes<-unique(unlist(lapply(names(table(metadata$diagnosis)),function(group){
  tmp<-data.exp[,rownames(metadata[metadata$diagnosis==group,])]
  nonVar.genes<-summ.var(data = tmp)
  tmp<-tmp[!nonVar.genes$nzv,]
  return(rownames(tmp))
})))
keepGenes<-keepGenes[keepGenes!="" & !is.na(keepGenes)]

data.exp<-data.exp[keepGenes,]

# ······································································· Step 3
## Get single sample scores

#' @param inputData Matrix, data frame, ExpressionSet or SummarizedExperiment
#' with omics data. Feature names must match the gene sets nomenclature. To use
#' preloaded databases, they must be gene symbols.
#' @param geneSets A named list with each gene set,
#' or the name of one preloaded database (go_bp, go_cc, go_mf, kegg, reactome,
#' pharmgkb, lincs, ctd, disgenet, hpo, wikipathways, tmod)
#' or a GeneSetCollection. For using network methods,
#' a data frame including columns:
#' "source","target","weight" and "mor" (optional).
#' @param method Scoring method: M-Scores, GSVA, ssGSEA, singscore, Plage,
#' Z-score, AUCell, MDT, MLM, ORA, UDT, ULM, FGSEA, norm_FGSEA, WMEAN,
#' norm_WMEAN, corr_WMEAN, WSUM, norm_WSUM or corr_WSUM.
#' @param labels (Only for M-Scores) Vector with the samples class labels (0 or
#' "Healthy" for control samples). Optional.
#' @param cores Number of cores to be used.
#' @param use.assay If SummarizedExperiments are used, the number of the assay 
#' to extract the data.
#' @param ... Additional parameters for the scoring functions.
#'

labels<-metadata$diagnosis
names(labels)<-rownames(metadata)

score.data<-getScores(inputData =data.exp,
                      geneSets = "reactome",
                      method = "GSVA",
                      labels = labels, # Only for M-scores
                      cores = 1) 

#' score.data is a data.frame with samples in columns, gene sets in rows, and
#' scores as matrix entries

