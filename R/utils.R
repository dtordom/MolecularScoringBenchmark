##############################
## PMPS Bechmark 
## R version 4.3.3 (2024-02-29 ucrt)
## 
##############################
## utils.R

##····························································· check.packages()
## Check installed packages and load packages
#' @param pkg Character vector contains package names
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)>0){
    cat(paste0("\nInstalling ", paste(new.pkg,collapse = ", ")))
    if (!requireNamespace("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")
    }
    bioconductor_packages <- BiocManager::available()
    bioconductor_packages <- bioconductor_packages[bioconductor_packages %in% new.pkg]
    if (length(bioconductor_packages) > 0){
      BiocManager::install(bioconductor_packages, dependencies = TRUE,ask=FALSE)
    }
    new.pkg = new.pkg[!(new.pkg %in% bioconductor_packages)]
    if (length(new.pkg)>0){
      install.packages(new.pkg, dependencies = TRUE)
    }
  }
  res <- lapply(pkg,function(x){
     cat(paste0("\nLoading ",x))
     suppressMessages(require(x,character.only = T))
  })
} 
##··············································································

##····································································norm.log()
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
##··············································································

##····································································summ.var()
## Summary of gene variance
#' @param data Gene expression data.frame
summ.var<-function(data){ ## Gene expression matrix
  require(caret)
  nzv <- nearZeroVar(t(data), saveMetrics= TRUE)
  return(nzv)
}
##··············································································

##·····························································mergeExpression()
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
##··············································································

##·······························································annotateGenes()
## Function to annotate genes from a gene-expression matrix 
#' @param data Gene expression data.frame
#' @param toGenes Identifier of db to annotate genes
#' @param fromGenes Identifier of db of probes/genes used in data
#' @param method Method used to merge expression of several probes pointing to 
#' the same gene. "median", "mean", "max", "sum"
#' @param ensembl ensembl object from biomaRt package
annotateGenes<-function(data,
                        toGenes='external_gene_name',
                        fromGenes='ensembl_gene_id',
                        method = "median",
                        ensembl = NULL){
  require("biomaRt")
  require("parallel")
  require("tidyr")
  
  ## Get Gene annotation database
  if(is.null(ensembl)){
    
    ensembl = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  }
  genome = getBM(attributes=c(toGenes,fromGenes),mart = ensembl)
  genome <- genome %>% `colnames<-`(c("toGenes","fromGenes")) %>% 
    replace(.=="",NA) %>% drop_na() %>% filter(fromGenes %in% rownames(data))
    
  data = data[genome$fromGenes,]
  finalGenes = unique(genome$toGenes)
  nCores<-ifelse(as.character(Sys.info()["sysname"])=="Windows",1,detectCores())

  temp = as.data.frame(do.call("rbind", mclapply(finalGenes, mergeExpression,
                                                 genome = genome,
                                                 expressionMatrix = data,
                                                 method = method,
                                                 mc.cores = nCores)))
  rownames(temp) = finalGenes
  colnames(temp) = colnames(data)
  
  return(temp)
}
##··············································································

##································································data_summary()
## Get data summary (mean and sd) for metadata variables
#' @param data data table
#' @param varname variable to extract statistics
#' @param groupnames variable to group the data
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
##··············································································

##···································································applyGSVA()
## Get GSVA scores by sample
#' @param data Gene expression data.frame
#' @param exMatrix Gene expression matrix
#' @param geneSets list of pathways
#' @param minSize minimum number of genes to select a pathway as suitable
#' @param kcdf: "Gaussian" for normalized data and "Poisson" for counts
#' @param internal_n_cores Cores used to parellelize the process
applyGSVA <- function(exMatrix, geneSets, minSize = 3, kcdf = "Gaussian", internal_n_cores = 1){
  paramMatrix <- gsvaParam(exMatrix, geneSets, minSize =  minSize, kcdf = kcdf)
  gsvaMatrix <- gsva(paramMatrix, BPPARAM = MulticoreParam(workers = internal_n_cores))
  return(gsvaMatrix)
}
##··············································································

##··································································applyPlage()
## Get Plage scores by sample
#' @param data Gene expression data.frame
#' @param exMatrix Gene expression matrix
#' @param geneSets list of pathways
#' @param minSize minimum number of genes to select a pathway as suitable
#' @param kcdf: "Gaussian" for normalized data and "Poisson" for counts
#' @param internal_n_cores Cores used to parellelize the process
applyPlage <- function(exMatrix, geneSets, minSize = 3, kcdf = "Gaussian", internal_n_cores = 1){
  paramMatrix <- plageParam(exMatrix, geneSets, minSize =  minSize)
  plageMatrix <- gsva(paramMatrix, BPPARAM = MulticoreParam(workers = internal_n_cores))
  return(plageMatrix)
}
##··············································································

##·································································applyZscore()
## Get Zscores by sample
#' @param exMatrix Gene expression matrix
#' @param geneSets list of pathways
#' @param minSize minimum number of genes to select a pathway as suitable
#' @param internal_n_cores Cores used to parellelize the process
applyZscore <- function(exMatrix, geneSets, minSize = 3, internal_n_cores = 1){
  paramMatrix <- zscoreParam(exMatrix, geneSets, minSize =  minSize)
  ZscoreMatrix <- gsva(paramMatrix, BPPARAM = MulticoreParam(workers = internal_n_cores))
  #Scale data between -1 and 1
  #ZscoreMatrix <- ZscoreMatrix / max(abs(ZscoreMatrix))
  return(ZscoreMatrix)
}
##··············································································

##·································································applyssGSEA()
## Get ssGSEA scores by sample
#' @param exMatrix Gene expression matrix
#' @param geneSets list of pathways
#' @param minSize minimum number of genes to select a pathway as suitable
#' @param internal_n_cores Cores used to parellelize the process
applyssGSEA <- function(exMatrix, geneSets, minSize = 3, internal_n_cores = 1){
  paramMatrix <- ssgseaParam(exMatrix, geneSets, minSize =  minSize)
  ssGSEAMatrix <- gsva(paramMatrix, BPPARAM = MulticoreParam(workers = internal_n_cores))
  #Scale data between -1 and 1
  #ssGSEAMatrix <- ssGSEAMatrix / max(abs(ssGSEAMatrix))
  return(ssGSEAMatrix)
}
##··············································································

##······························································applySingscore()
## Get Singscores by sample
#' @param exMatrix Gene expression matrix
#' @param geneSets list of pathways
#' @param minSize minimum number of genes to select a pathway as suitable
#' @param internal_n_cores Cores used to parellelize the process
applySingscore <- function(exMatrix, geneSets, minSize = 3, internal_n_cores = 1){
  exMatrixgenes <- rownames(exMatrix)
  pathways_filtered <- lapply(geneSets, function(pathway) {
    genes_com <- intersect(pathway, exMatrixgenes)
    if (length(genes_com) >= minSize) {
      return(pathway)
    } else {
      return(NULL)
    }
  })
  pathways_filtered <- pathways_filtered[!sapply(pathways_filtered, is.null)]
  geneSets <- pathways_filtered
  rankMatrix <- rankGenes(exMatrix, tiesMethod = "average")
  # Save currect options
  op <- pboptions()
  # New barr type
  pboptions(type = "txt", style = 3, char = "=")
  print("Estimating singscore values")
  listSign <- suppressWarnings(pblapply(geneSets,
                                        function(x) simpleScore(rankData = rankMatrix, upSet = x)))
  # Back to old optiops
  pboptions(op)
  listScores <- sapply(listSign, function(x) x$TotalScore)
  if(is.list(listScores)){
    pathMatrix <- do.call(rbind, listScores)}
  else{
    pathMatrix <- t(listScores)
  }
  colnames(pathMatrix) <- colnames(exMatrix)
  #Escale data between -1 and 1
  #pathMatrix <- pathMatrix / 0.5
  return(pathMatrix)
}
##··············································································

##·····································································varTest()
## Variance test between groups
#' @param data Gene expression matrix
#' @param group vector of groups
#' @param ref reference id (character) of group
#' @param norm TRUE/FALSE, are the data normal?
#' @param adj.pvalue method to correct p values
varTest<-function(data,group,ref="HC",norm=FALSE,adj.pvalue = "fdr"){
  
  RES<-as.data.frame(matrix(data=0,ncol=5,nrow=nrow(data)))
  colnames(RES)<-c("Ref_var","Case_var","Pvalue","adj.Pvalue","Difference")
  RES[,"Ref_var"]<-apply(data[,group==ref],1,var)
  RES[,"Case_var"]<-apply(data[,group!=ref],1,var)
  RES[,"Difference"]<-RES$Case_var - RES$Ref_var
  
  for(gene in 1:nrow(data)){
    
    if(sd(as.numeric(data[gene,]),na.rm = T)!=0){
      if(norm){
        res<-bartlett.test(list(as.numeric(data[gene,group==ref]),
                                as.numeric(data[gene,group!=ref])))
        res<-res$p.value
      }else{
        res<-leveneTest(y = as.numeric(data[gene,]),
                        group = as.factor(unname(group)),center = "median")
        res<-res$`Pr(>F)`[1]
      }
      RES[gene,"Pvalue"]<-res
    }else{
      RES[gene,"Pvalue"]<-1
    }
  }
  RES$adj.Pvalue<-p.adjust(as.numeric(RES$Pvalue), method = "bonferroni",n = nrow(RES))
  
  RES$color<-rep("#99999999",nrow(RES))
  RES$color[RES$adj.Pvalue<=0.05 & RES$Difference>0]<-"lightgreen"
  RES$color[RES$adj.Pvalue<=0.05 & RES$Difference<0]<-"lightblue"
  rownames(RES)<-rownames(data)
  
  return(RES)
}
##··············································································

##······························································vectSimilarity()
## Return similarity measurements between two vectors
#' @param x numeric or categorical vector 
#' @param y numeric or categorical vector 
#' @param method one or several metrics of similarity: "euclidean","maximum",
#' "canberra","binary","minkowski","pearson","spearman","kendall","cosine",
#' "jaccard", "cohenk","adjrand"
vectSimilarity<-function(x,y,method="euclidean"){
  
  require("jaccard")
  require("psych")
  require("lsa")
  require("mclust")
  
  ## Build for other similarity metrics
  results<-NULL
  
  if(any(c("euclidean","maximum","canberra","binary","minkowski") %in% method)){
    tmpMethods<-method[method %in% c("euclidean","maximum","canberra","binary","minkowski")]
    res<-unlist(lapply(tmpMethods,function(met){
      dist(rbind(x,y),method = met)}))
    names(res)<-tmpMethods
    results<-c(results,res)
  }
  
  if(any(c("pearson","spearman","kendall") %in% method)){
    tmpMethods<-method[method %in% c("pearson","spearman","kendall")]
    res<-unlist(lapply(tmpMethods,function(met){
      cor.test(x,y,method = met)$estimate}))
    names(res)<-tmpMethods
    results<-c(results,res)
  }
  
  if("cosine" %in% method){
    results<-c(results,cosine(x,y))
    names(results)[length(results)]<-"cosine"
  }
  
  if("jaccard" %in% method){
    results<-c(results,jaccard(x,y))
    names(results)[length(results)]<-"jaccard"
  }
  if("cohenk" %in% method){
    results<-c(results,cohen.kappa(cbind(x,y))$kappa)
    names(results)[length(results)]<-"cohenk"
  }
  if("adjrand" %in% method){
    results<-c(results,adjustedRandIndex(x,y))
    names(results)[length(results)]<-"adjrand"
  }
  
  return(results)
}


##··············································································

##······································································getGap() 
## Get Gaps for heatmaps to separe clusters 
#' @vect numeric ordered vector
getGap<-function(vect){
  isum<-0
  res<-NULL
  for(i in 1:(length(vect)-1)){
    isum<-isum+vect[i]
    res<-c(res,isum)
  }
  return(res)
}
##··············································································

##···································································testAssoc()
#' @param metadata data.frame with variables in columns, First column: "Group"
#' is used to subgroup the data. Samples/ observations in rows  
#' @param p.adj pvalue adjustment:"holm", "hochberg", "hommel", "bonferroni", 
#' "BH", "BY","fdr", "none")
#' @param remove.vars Variable value (categorical) to remove for the analysis
#' @param colors named Vector specifying colors for variable ids
#' i,e: colors = c("Case"="red","Healthy"="blue","cl1"="grey",etc)
#' @param normality bool to indicate normality on data, if it is NULL, 
#' shaphiro.test is run internally to define normality
testAssoc<-function(metadata,
                    p.adj = "none",
                    remove.vars=NULL,
                    colors=NULL,
                    normality=NULL,
                    positiveClass=NULL){
  
  require("ggplot2")
  require("scales")
  require("dplyr")
  require("reshape2")
  require("rstatix")
  
  theme_ggplot<-theme(axis.title.x = element_text(size = 8),
                      axis.title.y = element_text(size = 8),
                      axis.text.y = element_text(size = 6),
                      axis.text.x = element_text(size = 6),
                      axis.ticks = element_line(size = 0.25),
                      axis.line=element_line(linewidth = 0.25),
                      plot.title = element_text(size = 8))
  
  
  RESULTS<-lapply(2:ncol(metadata),function(i){
    
    #print(i)
    met<-metadata[,c(1,i)]
    colnames(met)<-c("group","variable")
    if(!is.null(remove.vars)){ met<-met[!met$variable %in% remove.vars,] }
    met<-met[!is.na(met[,2]),]
    
    if(nrow(met)<3 | dim(table(met$group))<2 | dim(table(met$variable))<2){
      return(NULL)
    }else{
      
      if(is.numeric(met$variable)){ # Numerical variables ·························
        
        met.stats<- as.data.frame(met %>% group_by(group) %>% 
                                    summarise(mean = mean(variable,na.rm = T),
                                              sd = sd(variable,na.rm = T)))
        ## Check normality in data
        if(is.null(normality)){
          normality<-ifelse(shapiro.test(met[,2])$p.value<=0.05,F,T)
        }
        
        coords<-NULL
        
        if(normality){
          info<-c("type"="num","general.test"="Anova","posthosc"="Pairwise t-test")
          
          ## ANOVA test
          all.pvalue <- aov(met$variable ~ met$group)
          all.pvalue<-do.call("cbind",summary(all.pvalue))
          colnames(all.pvalue)[colnames(all.pvalue)=="Pr(>F)"]<-"pvalue"
          
          ## T-test by groups
          ph.pvalue<-as.data.frame(pairwise_t_test(met,variable ~ group, p.adjust.method = p.adj))
          
          ## Plot parameters (coords)
          if(min(ph.pvalue$p)<=0.05){
            coords<-as.data.frame(ph.pvalue[ph.pvalue$p<=0.05,c("group1","group2","p")])
            colnames(coords)[3]<-"pvalue"
            coords$x1<-gsub("[^0-9]","",coords$group1)
            coords$x2<-gsub("[^0-9]","",coords$group2)
            coords$group<-coords$group1
            coords$ypos<-seq(from=max(met$variable),
                             by=round(max(met$variable)*0.05,digits = 0),
                             length.out=nrow(coords))
            coords$p<-format(coords$p,scientific = TRUE,digits = 3)
          }
          
        }else{
          info<-c("type"="num","general.test"="Kruskal-wallis","posthosc"="Wilcox-test")
          
          ## Kruskal Wallis test
          all.pvalue<-kruskal.test(variable ~ group, data = met)
          all.pvalue<-data.frame("statistic"=all.pvalue$statistic,
                                 "df"=all.pvalue$parameter,
                                 "pvalue"=all.pvalue$p.value)
          
          ## Wilcox test
          ph.pvalue<-pairwise.wilcox.test(met$variable, met$group,p.adjust.method = p.adj)
          ph.pvalue<-as.data.frame(ph.pvalue$p.value)
          
          ## Plot parameters (coords)
          if(min(ph.pvalue,na.rm=T)<=0.05){
            
            coords<-data.frame("group1"=rep(colnames(ph.pvalue),each=ncol(ph.pvalue)),
                               "group2"=rep(rownames(ph.pvalue),nrow(ph.pvalue)),
                               "pvalue"=as.numeric(unlist(as.vector(ph.pvalue))))
            coords <- coords %>% filter(!is.na(pvalue)) %>% filter(pvalue<=0.05)
            coords$x1<-gsub("[^0-9]","",coords$group1)
            coords$x2<-gsub("[^0-9]","",coords$group2)
            coords$group<-coords$group1
            coords$ypos<-seq(from=max(met$variable),
                             by=round(max(met$variable)*0.05,digits = 0),
                             length.out=nrow(coords))
            
            coords$pvalue<-format(coords$pvalue,scientific = TRUE,digits = 3)
          }
          
        }
        
        
        ## Plotting results
        p1<-NULL
        
        pvaltext<-ifelse(!is.na(all.pvalue$pvalue[1]),
                         format(all.pvalue$pvalue[1],scientific = TRUE,digits = 3),"")
        
        p1<-ggplot(met,aes(x=group,y=variable,fill=group))+theme_classic()+
          geom_jitter(size=1,alpha=0.8,shape=21,color="black",stroke=0.2,width = 0.15)+
          theme_ggplot+
          geom_boxplot(lwd=0.2,alpha=0.8,outlier.colour = NA,width = 0.5,alpha=0.8)+ 
          labs(x="Groups",
               y=colnames(metadata)[i],
               title=paste0(colnames(metadata)[i],": (",info[2]," = ",pvaltext,")"))
        
        if(is.null(colors)){
          p1<-p1 + scale_fill_brewer(palette="Set2")
        }else{
          p1<-p1 + scale_fill_manual(values=colors)
        }
        
        ## Add pvalues for pair comparisons
        if(!is.null(coords)){
          for(p in 1:nrow(coords)){
            p1 <- p1 + geom_segment(data = coords[p, ],aes(x =group1,
                                                           xend =group2,y = ypos, yend = ypos),
                                    size =0.1, color ="black") +
              geom_text(data = coords[p,],aes(x = (as.numeric(x1)+as.numeric(x2))/2, y= ypos + (max(met$variable)*0.025),
                                              label = pvalue),size=2)
          }
        }
        
        res<-list(met.stats,all.pvalue,ph.pvalue,info,p1)
        names(res)<-c("table","general.test","posthosc","info","plot")
        
      }else{ # Categorical variables ·························
        met$variable<-as.character(met$variable)
        
        #' Function can be modify to add chisq.test when sample size >5
        #' But results are similar and more stringent for fisher.test
        #' Only computational time is reduced using chisq in big data 
        
        fisher_table <- table(met)
        fisher_table[is.na(fisher_table)]<-0
        
        info<-c("type"="cat","general.test"="fisher","posthosc"="fisher")
        
        ## General pvalue
        all.pvalue<-min(fisher.test(fisher_table,alternative = "greater",simulate.p.value=TRUE)$p.value,
                        fisher.test(fisher_table,alternative = "less",simulate.p.value=TRUE)$p.value)
        
        ## Pvalue One vs Others
        if(length(unique(met$variable))>2){
          
          pvals<-unlist(lapply(unique(met$variable),function(g){
            tmp.met<-met
            tmp.met$variable<-ifelse(tmp.met$variable==g,g,"other")
            tmp.fisher_table <- table(tmp.met)
            
            pvalue<-min(fisher.test(tmp.fisher_table,alternative = "greater")$p.value,
                        fisher.test(tmp.fisher_table,alternative = "less")$p.value)
            return(pvalue)
          }))
          ph.pvalue<-data.frame("group"=unique(met$variable),"pvalue"=pvals)
        }else{
          ph.pvalue<-NULL
        }
        
        
        ## Formatting Fisher table
        fisher_table<-as.data.frame.matrix(fisher_table)
        if(is.null(positiveClass)){
          positiveClass<-colnames(fisher_table)[which.max(apply(fisher_table,2,function(cl){sum(cl,na.rm=T)}))]
        }else{
          if(!positiveClass %in% colnames(fisher_table)){
            
            positiveClass<-colnames(fisher_table)[which.max(apply(fisher_table,2,function(cl){sum(cl,na.rm=T)}))]
          }
        }
        
        ## Plotting results  
        p1<-NULL
        
        nSize<-rep(rowSums(fisher_table),ncol(fisher_table))
        
        fisher_pl<-cbind(melt(fisher_table),"Cluster"=rep(rownames(fisher_table),ncol(fisher_table)))
        fisher_pl$value<-(fisher_pl$value / nSize)*100
        colnames(fisher_pl)<-c("Group","Freq","Cluster")
        rownames(fisher_pl)<-NULL
        
        # fisher_pl<-data.frame("Cluster"=rownames(fisher_table),fisher_table,
        #                       "Freq"=as.numeric((fisher_table[,positiveClass]/
        #                                            (rowSums(fisher_table)))*100))
        # 
        # fisher_pl<-rbind(fisher_pl[,c("Cluster","Freq")],
        #                  data.frame("Cluster"=rownames(fisher_table),
        #                             "Freq"=as.numeric(rowSums(fisher_table[,!colnames(fisher_table) %in% positiveClass])/
        #                                                 (rowSums(fisher_table))*100)))
        # if(ncol(fisher_table)>2){
        #   fisher_pl$Group = c(rep(positiveClass,nrow(fisher_table)),
        #                       rep("Others",nrow(fisher_table)))
        # }else{
        #   fisher_pl$Group = c(rep(positiveClass,nrow(fisher_table)),
        #                       rep(colnames(fisher_table)[colnames(fisher_table)!=positiveClass],nrow(fisher_table)))
        # }
        # 
        
        
        fisher_pl <- fisher_pl %>% group_by(Cluster) %>%
          mutate(Prop = round((Freq / sum(Freq))*100,digits = 2))
        fisher_pl$Prop<- round(fisher_pl$Freq)
        
        p1<-ggplot(fisher_pl,aes(x=Cluster,fill=Group,weight=Freq,by=Cluster))+
          geom_bar(position="fill",color="black",linewidth=0.3) + theme_minimal()+
          theme_ggplot +
          scale_y_continuous(labels =label_percent(scale=100))+
          labs(x="Groups",
               y="Proportion",
               title=paste0(colnames(metadata)[i],": (pval: ",
                            format(all.pvalue,scientific = TRUE,digits = 3),")"))
        
        ## Add colors
        if(is.null(colors)){
          p1<-p1 + scale_fill_brewer(palette="Set2")
        }else{
          p1<-p1 + scale_fill_manual(values=colors)
        }
        
        
        ## Add text into the boxes
        p1 <- p1 + geom_text(
          aes(label = sprintf("%.2f", Freq), y = ..count.. / sum(..count..)),
          stat = "count",
          position = position_fill(vjust = 0.5),
          color = "black",
          size = 2
        )
        
        ## Return results
        res<-list(fisher_table,all.pvalue,ph.pvalue,info,p1)
        names(res)<-c("table","general.test","posthosc","info","plot")
        
      } # Categorical variables ·························
      
      ## lapply return
      return(res)
    }
    
  })
  names(RESULTS)<-colnames(metadata)[2:ncol(metadata)]
  
  #RESULTS<-Filter(Negate(is.null),RESULTS)
  return(RESULTS)
}

##··············································································

##································································orderSamples() 
## This function is used to sort patients for heatmap visualization
#' @data gene expression data
#' @groups vector with patient assignation to clusters (patient ids as names)
orderSamples<-function(data,groups){
 
  require("pheatmap")
  pats<-NULL
  group.id<-unique(groups)
  
  for(i in group.id){
    
    dat.tmp<-data[,names(groups)[groups==i]]
    p<-pheatmap(dat.tmp,scale="none",show_colnames = F,cluster_cols = T,
                cluster_rows = T, show_rownames = F, border_color = NA,
                breaks=seq(-2,2,length.out = 100), fontsize = 5,
                color = colorRampPalette(c("deepskyblue4","white","coral2"))(100))
    dev.off()
    pats<-c(pats,colnames(dat.tmp)[p$tree_col$order])
  }
    
  return(pats)
}
##··············································································

##···························································mergeListMatrices() 
## merge matrices with different dimensions from a list, keeping all rows
#' @mlist list with matrices
#' @value value to put in non-shared rows
mergeListMatrices<-function(mlist, value =NA){
  
  rowNames<-unique(c(unlist(sapply(mlist,function(x){rownames(x)}))))
  colNames<-c(unlist(sapply(mlist,function(x){colnames(x)})))
  newMatrix<-as.data.frame(matrix(value,nrow=length(rowNames),ncol=sum(sapply(mlist,ncol))))
  rownames(newMatrix)<-rowNames
  colnames(newMatrix)<-colNames
  
  for(i in 1:length(mlist)){
    matrix.i<-mlist[[i]]
    newMatrix[rownames(matrix.i),colnames(matrix.i)]<-matrix.i
  }
  return(newMatrix)
}

##··············································································

##······························································· ClusterStats()
## Get cluster stability
#' @param list.i 
#' @param Dist.list 
#' @param combinations
ClusterStats<-function(list.i, ## permutation list
                       Dist.list, ## List of distance matrices
                       combinations){ ## Matrix with parameter combinations
  require("SNFtool")
  require("NbClust")
  cat("Running ClusterStats")
  
  nDatasets<-as.numeric(unlist(list.i))
  Times = 15; 	# Number of Iterations, usually (10~20)
  
  ListaNB<-lapply(1:nrow(combinations),function(i){
    
    ## Parameter and dataset selection
    K<<-combinations[i,1]
    Alpha<<-combinations[i,2]
    sel.samples<-sample(1:length(Dist.list),nDatasets,replace=F)
    
    tmp.list<-list()
    for(l.i in 1:length(Dist.list)){
      data<-Dist.list[[l.i]]
      w<-affinityMatrix(data,K,Alpha)
      tmp.list[[l.i]]<-w
    }
    tmp.list<-tmp.list[sel.samples]
    
    W = SNF(tmp.list, K, Times)
    Nb<-NbClust(data = W, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 15,
                method = "complete", index = "ch", alphaBeale = 0.1)
    
    return(Nb$All.index)
  })
  
  ## Get summary
  nbclusts<-lapply(ListaNB,function(x){
    res<-as.numeric(names(x[order(x,decreasing=T)][1]))
  })
  nbclusts<-unlist(nbclusts)
  res<-table(nbclusts)
  
  res<-c(as.numeric(res["2"]),as.numeric(res["3"]),as.numeric(res["4"]),as.numeric(res["5"]),
         as.numeric(res["6"]),as.numeric(res["7"]),as.numeric(res["8"]),as.numeric(res["9"]),
         as.numeric(res["10"]),as.numeric(res["11"]),as.numeric(res["12"]),as.numeric(res["13"]),
         as.numeric(res["14"]),as.numeric(res["15"]))
  res[is.na(res)]<-0
  
  res<-data.frame(cbind(rep(nDatasets,14),2:15,res))
  colnames(res)<-c("datasets","clusters","nbindex")
  
  
  return(res)
  
}

##··············································································

##································································ concordance()
## Get concordance between cluster results across iteractions
#' @param df dataframe with samples in rows and cluster assigment in each 
#' iteraction (column). First column is used as reference 
concordance <- function(df) {
  df[,1]<-paste0("CL_",df[,1])
  for(j in 2:ncol(df)){
    
    cls<-table(df[,c(1,j)])
    cls<-apply(cls,2,function(x) names(which.max(x)))
    df[,j] <- cls[as.character(df[,j])]
  }
  return(df)
}

##··············································································
