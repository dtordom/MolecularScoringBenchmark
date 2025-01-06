################################################################################
## Personalized molecular profile scoring benchmark
## R version 4.3.1 (2024-04-10)
## 
################################################################################
#' Testing NA and sample size influence on single sample scoring


# devtools::install_github("jordimartorell/pathMED")
##······································································· Step 0
## Set environment

setwd("D:/Work/BMS/PMPS/PMSPbenchmark")
set.seed(1234)

load(paste0(getwd(),"/RData/SCORES.RData"))

source(paste0(getwd(),"/code/utils.R"))

check.packages(c("vctrs","stats","stringr","caret","singscore","GSVA","parallel","lsa",
                 "metrica","BiocParallel","pbapply","reshape","ggplot2","ggpubr","psych",
                 "jaccard","NbClust","ConsensusClusterPlus","mclust","decoupleR",
                 "reshape2","qvalue","ggbreak","pheatmap","pathMED","scales","dplyr"))


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

##······································································· Step 1
## Q1. How similar the score results are between them?  
#' Pearson correlation are calculated for each pair of methods and for each sample
#' Pearson correlation average across patients was computed as final similarity 
#' measure between pair of methods (using PRECISESADS cohort)

Cor.table<-as.data.frame(matrix(data=1,ncol=length(methods),nrow=length(methods)))
colnames(Cor.table)<-methods
rownames(Cor.table)<-methods

Dis.table<-as.data.frame(matrix(data=1,ncol=length(methods),nrow=length(methods)))
colnames(Dis.table)<-methods
rownames(Dis.table)<-methods

combinations<-combn(methods,2)

for (pat in 1:ncol(SCORES$prec$`M-score`)){ ## loop by patient
  for(metr in 1:ncol(combinations)){ ## compare metrics
    paths<-intersect(rownames(SCORES$prec[[combinations[1,metr]]]),
                     rownames(SCORES$prec[[combinations[2,metr]]]))
    x<-cbind(SCORES$prec[[combinations[1,metr]]][paths,pat],
             SCORES$prec[[combinations[2,metr]]][paths,pat])
    x[is.infinite(x)] <- NA
    
    res<-vectSimilarity(x=x[,1],y=x[,2],method = c("euclidean","pearson"))
    
    Cor.table[combinations[1,metr],combinations[2,metr]]<-res["pearson"]
    Cor.table[combinations[2,metr],combinations[1,metr]]<-res["pearson"]
    Dis.table[combinations[1,metr],combinations[2,metr]]<-res["euclidean"]
    Dis.table[combinations[2,metr],combinations[1,metr]]<-res["euclidean"]
  }
}

pheatmap(Cor.table,scale="none",show_colnames = T,cluster_cols = T,
         cluster_rows = T, show_rownames = T, border_color = "black", fontsize = 7,
         breaks=seq(-1,1,length.out = 100),#gaps_col =ncol(m2),
         color = colorRampPalette(c("darkblue","white","darkgreen"))(100))

rm(Dis.table,Cor.table,pat,metr,combinations,x,res,paths)

##······································································· Step 2
## Q1.1: Are the scores influenced by the cohort composition?
#' Cohort composition: number of patients and heterogeneity of the cohort
#' Using PRECISESADS data, test: 1,5,10,20,30,40,50,60,70,80,90


## Use database subsets of 250 pats (to reduce computational time)
tmp.geneset<-geneset.DB[sample(1:length(geneset.DB),250)]

## Build a small reference (250 paths instead 4.4k paths, to reduce computational time)
REFs<-lapply(methods,function(method){
    cat(paste0("\nRuning ",method,"\n"))
    tmp.SLE<-SLE.prec
    labs<-rep("SLE",ncol(tmp.SLE))
    
    if(method=="M-score"){ ## Get Mscores also for HC
      dat<-cbind(tmp.SLE,HC.prec)
      labs<-c(labs,rep("Healthy",ncol(HC.prec)))
    }else{
      dat<-tmp.SLE
    }
    ncores<-ifelse(method %in% c("M-scores", "MDT", "UDT", "ULM"), 14, 1)
    
    tmp.scores<-getScores(inputData=dat, geneSets=tmp.geneset, method = method,
                          labels = labs,cores = ncores, minSize=3,
                          times =ifelse(grepl("norm",method),100,
                                        ifelse(grepl("corr",method),200,2)))
    return(tmp.scores)
  })
names(REFs)<-methods



## Run influence of sample size on scores
sizes<-rep(trunc(ncol(SLE.prec)*(c(1,5,10,20,30,40,50,60,70,80,90)/100)),each = 10)

Results.cc<-as.data.frame(do.call("rbind",lapply(1:length(sizes),function(i){
  
  ## Select random patients within each iteration and cohort size)
  tmp.SLE<-SLE.prec[,sample(1:ncol(SLE.prec),sizes[i]),drop=FALSE]
  labs <- rep("SLE",ncol(tmp.SLE))
  
  res<-as.data.frame(do.call("rbind",lapply(methods,function(method){
    
    ncores<-ifelse(method %in% c("M-scores", "MDT", "UDT", "ULM"), 14, 1)
    
    if(method=="M-score"){ ## Get Mscores also for HC
      dat<-cbind(tmp.SLE,HC.prec)
      labs<-c(labs,rep("Healthy",ncol(HC.prec)))
    }else{
      dat<-tmp.SLE
    }
    
    tmp.sc<-getScores(inputData=dat, geneSets=tmp.geneset, method=method,
                      labels=labs, cores=ncores, minSize=3,
                      times = ifelse(grepl("norm",method),100, ifelse(grepl("corr",method),200,2)))
    
    ## Mean similarity across pair of patient comparisons
    sim.vals<-as.data.frame(do.call("rbind",lapply(1:ncol(tmp.sc),function(pat){
      sharedPaths<-intersect(rownames(REFs[[method]]),rownames(tmp.sc))
      x<-tmp.sc[sharedPaths,pat]
      y<-REFs[[method]][sharedPaths,colnames(tmp.sc)[pat]]
      df<-cbind(x,y)
      df <- df[apply(df, 1, function(row) all(is.finite(row))), ]
      sim.values.p<-vectSimilarity(x=df[,1],y=df[,2],method = c("euclidean","pearson"))
      sim.values.p[is.na(sim.values.p)]<-0
      return(sim.values.p)
    })))
    sim.vals<-apply(sim.vals,2,function(col){mean(col,na.rm=T)})
    return(sim.vals)
  })))
  
  res$method<-methods
  res$size<-sizes[i]
  print(res)
  return(res)
})))

save.image(paste0(getwd(),"/RData/Results_CC1.RData"))


## Plots ············
library(ggplot2)
library(dplyr)

Results.cc<-Results.cc[!grepl("corr",Results.cc$method),]

Results.cc$x<-as.numeric(unlist(lapply(1:nrow(Results.cc),function(i){
  round((Results.cc[i,"size"]*100)/ncol(SLE.prec),digits = 0)
})))

df <- Results.cc %>%  mutate(x = factor(x, levels = sort(unique(x)))) %>%
  mutate(size = factor(size, levels = sort(unique(size))))

df$method <- factor(df$method, levels = methods)

p.cc.eu<-ggplot(df, aes(x = x, y = euclidean,fill=method)) +
  geom_jitter(aes(color = "black",), shape = 21, color = "black", width = 0.2, size = 1, alpha = 0.8) +
  geom_boxplot(aes(color = "black"), outlier.shape = NA,alpha=1) + 
  #ylim(min(df$euclidean),100) +
  facet_wrap(~method, scale = "free_x",ncol=5,nrow=4) +  # free
  labs(x = "", y = "Euclidean distance") +
  theme_classic() +
  theme(legend.text = element_text(size=8),
        axis.text.x=element_text(size=7, color = "black",hjust=1, vjust=0.5,angle = 90), #,angle = 90
        axis.text.y = element_text(size=7, color = "black"),
        axis.title = element_text(size=8, color = "black", face = "bold"), 
        strip.text = element_text(size=8,face = "bold"), 
        #strip.background = element_rect(fill = "white", color = "black", size = 0.5),  
        strip.background = element_blank(),
        plot.title = element_text(size=10,hjust=0.5,face="bold"))+
  theme(legend.position = "none") +
  scale_color_manual(values="black")+
  scale_fill_manual(values=colors)
plot(p.cc.eu)

p.cc.cor<-ggplot(df, aes(x = x, y = pearson,fill=method)) +
  geom_jitter(aes(color = "black",), shape = 21, color = "black", width = 0.2, size = 1, alpha = 0.8) +
  geom_boxplot(aes(color = "black"), outlier.shape = NA,alpha=1) + 
  ylim(min(df$pearson),1) +
  facet_wrap(~method, scales = "free_y",ncol=5,nrow=4) + 
  labs(x = "", y = "Pearson correlation") +
  theme_classic() +
  theme(legend.text = element_text(size=8),
        axis.text.x=element_text(size=7, color = "black",hjust=1, vjust=0.5,angle = 90), #
        axis.text.y = element_text(size=7, color = "black"),
        axis.title = element_text(size=8, color = "black", face = "bold"), 
        strip.text = element_text(size=8,face = "bold"), 
        #strip.background = element_rect(fill = "white", color = "black", size = 0.5),  
        strip.background = element_blank(),
        plot.title = element_text(size=10,hjust=0.5,face="bold"))+
  theme(legend.position = "none") +
  scale_color_manual(values="black")+
  scale_fill_manual(values=colors)
plot(p.cc.cor)


meth.cor<-c("GSVA","Z-score","Plage","MDT")
meth.eu<-c("GSVA","ssGSEA","Z-score","Plage","AUCell","MDT","ORA")

## Plot only the variable methods
p1.eu<-ggplot(df[df$method %in% meth.eu,], aes(x = x, y = euclidean,fill=method)) +
  geom_jitter(aes(color = "black",), shape = 21, color = "black", width = 0.2, size = 1, alpha = 0.8) +
  geom_boxplot(aes(color = "black"), outlier.shape = NA,alpha=1) + 
  #ylim(min(df$pearson),1) +
  facet_wrap(~method, scales = "free",ncol=4,nrow=4) +  # Separar por método
  labs(x = "Sample size (%)", y = "Euclidean distance") +
  theme_classic() +
  theme(legend.text = element_text(size=8),
        axis.text.x=element_text(size=7, color = "black",hjust=1, vjust=0.5,angle = 90), #
        axis.text.y = element_text(size=7, color = "black"),
        axis.title = element_text(size=8, color = "black", face = "bold"), 
        strip.text = element_text(size=8,face = "bold"), 
        #strip.background = element_rect(fill = "white", color = "black", size = 0.5),  
        strip.background = element_blank(),
        plot.title = element_text(size=10,hjust=0.5,face="bold"))+
  theme(legend.position = "none") +
  scale_color_manual(values="black")+
  scale_fill_manual(values=colors)
plot(p1.eu)

p1.corr<-ggplot(df[df$method %in% meth.cor,], aes(x = x, y = pearson,fill=method)) +
  geom_jitter(aes(color = "black",), shape = 21, color = "black", width = 0.2, size = 1, alpha = 0.8) +
  geom_boxplot(aes(color = "black"), outlier.shape = NA,alpha=1) + 
  #ylim(min(df$pearson),1) +
  facet_wrap(~method, scales = "free_y",ncol=5,nrow=4) +  # Separar por método
  labs(x = "Sample size (%)", y = "Pearson correlation") +
  theme_classic() +
  theme(legend.text = element_text(size=8),
        axis.text.x=element_text(size=7, color = "black",hjust=1, vjust=0.5,angle = 90), #
        axis.text.y = element_text(size=7, color = "black"),
        axis.title = element_text(size=8, color = "black", face = "bold"), 
        strip.text = element_text(size=8,face = "bold"), 
        #strip.background = element_rect(fill = "white", color = "black", size = 0.5),  
        strip.background = element_blank(),
        plot.title = element_text(size=10,hjust=0.5,face="bold"))+
  theme(legend.position = "none") +
  scale_color_manual(values="black")+
  scale_fill_manual(values=colors)
plot(p1.corr)



#········································································ STEP 3
## Q1.2: Could these changes on scores affect to patient clustering?

## Test only methods affected by cohort composition
methods<-unique(c(meth.cor,meth.eu))

## Use database subsets of 250 pats (to reduce computational time)
tmp.geneset<-geneset.DB[sample(1:length(geneset.DB),250)]

## Build a small reference (250 paths instead 4.4k paths, to reduce computational time)
REFs<-lapply(methods,function(method){
  cat(paste0("\nRuning ",method,"\n"))
  tmp.SLE<-SLE.prec
  labs<-rep("SLE",ncol(tmp.SLE))
  
  if(method=="M-score"){ ## Get Mscores also for HC
    dat<-cbind(tmp.SLE,HC.prec)
    labs<-c(labs,rep("Healthy",ncol(HC.prec)))
  }else{
    dat<-tmp.SLE
  }
  ncores<-ifelse(method %in% c("M-scores", "MDT", "UDT", "ULM"), 14, 1)
  
  tmp.scores<-getScores(inputData=dat, geneSets=tmp.geneset, method = method,
                        labels = labs,cores = ncores, minSize=3,
                        times =ifelse(grepl("norm",method),100,
                                      ifelse(grepl("corr",method),200,2)))
  return(tmp.scores)
})
names(REFs)<-methods


## Select 30-50 samples to use as core of samples / 10% of the samples
sampleCl<-colnames(SLE.prec)[sample(1:ncol(SLE.prec),50)]

sizes<-rep(trunc(ncol(SLE.prec)*(c(1,5,10,20,30,40,50,60,70,80)/100)),each = 5)


Results.fRate<-as.data.frame(do.call("rbind",lapply(1:length(sizes),function(i){
  
  ## Select random patients within each iteration and cohort size)
  tmp.SLE<-SLE.prec[,!colnames(SLE.prec) %in% sampleCl]
  tmp.SLE<-tmp.SLE[,sample(1:ncol(tmp.SLE),sizes[i]),drop=FALSE]
  labs <- rep("SLE",ncol(tmp.SLE))
  
  res<-as.numeric(unlist(lapply(methods,function(method){
    
    ncores<-ifelse(method %in% c("M-scores", "MDT", "UDT", "ULM"), 12, 1)
    
    if(method=="M-score"){
      tmp.SLE.m<-cbind(tmp.SLE,HC.prec)
      labs<-c(labs,rep("Healthy",ncol(HC.prec)))
    }else{
      tmp.SLE.m<-tmp.SLE
    }
    
    tmp.sc<-getScores(inputData=tmp.SLE.m, geneSets=tmp.geneset, method=method,
                      labels=labs, cores=ncores, minSize=3,
                      times = ifelse(grepl("norm",method),100, ifelse(grepl("corr",method),200,2)))
    
    sharedPaths<-intersect(rownames(tmp.sc),rownames(REFs[[method]]))
    
    dat1<-REFs[[method]][sharedPaths,c(sampleCl,intersect(colnames(tmp.sc),colnames(SLE.prec)))]
    dat2<-cbind(REFs[[method]][sharedPaths,sampleCl], tmp.sc[sharedPaths,])
    
    valid_rows <- apply(dat1, 1, function(row) all(is.finite(row))) & 
      apply(dat2, 1, function(row) all(is.finite(row)))
    dat1<-dat1[valid_rows,]
    dat2<-dat2[valid_rows,]
    
    ## Clustering (k = 10 max)
    cluster.res<-lapply(list(dat1,dat2),function(tmp.dat){
      
      d = sweep(tmp.dat,1, apply(tmp.dat,1,median,na.rm=T))
      clusters = ConsensusClusterPlus(as.matrix(d),maxK=10,reps=250,pItem=0.8,
                                      pFeature=1,clusterAlg="km",distance="euclidean",
                                      innerLinkage = "complete",seed=12345,plot=NULL)
      clusters<-do.call("cbind",lapply(2:10,function(k){
        return(clusters[[k]]$consensusClass)
      }))
      colnames(clusters)<-paste0("k",2:10)
      return(clusters)
    })
    
    ## Get failure rate (using different K)
    fRate<-c(unlist(lapply(1:ncol(cluster.res[[1]]),function(k){
      df<-as.data.frame(cbind(cluster.res[[1]][,k],cluster.res[[2]][,k]))
      df<-concordance(df)
      df<-df[colnames(tmp.sc),]
      frate <- (sum(df[,1] != df[,2])/nrow(df))*100
      return(frate)
    })))
    fRate<-round(mean(fRate,na.rm = T),digits = 3) #; print(fRate)
    
    return(fRate)
  })))
  
  res<-data.frame("method"=methods,"size"=sizes[i],"fRate"=res)
  print(res)
  return(res)
})))

save.image((paste0(getwd(),"/RData/results_frate.RData")))


## PLOT 
df <- Results.fRate

df <- df %>%
  group_by(size, method) %>%
  summarise(
    mean_frate = mean(fRate),
    sd_frate = sd(fRate),
    .groups = "drop"
  )

p.frate<-ggplot(df, aes(x = factor(size, levels = rev(unique(size))), y = mean_frate, fill = factor(method))) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8, color = "black", width = 0.8) +
  geom_errorbar(aes(ymin = mean_frate - sd_frate, ymax = mean_frate + sd_frate), 
                width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~method, ncol = 1) +
  labs(
    title = "Failure rate",
    x = "Sample size (%)",
    y = "Average of failure rate",
    fill = "Size"
  ) +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 8),
    axis.text.x = element_text(size = 7, color = "black", hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8, color = "black", face = "bold"),
    strip.text = element_text(size = 8, face = "bold"),
    strip.background = element_blank(),
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold")
  ) +
  theme(legend.position = "none") +
  coord_flip() + # Intercambiar los ejes
  scale_fill_manual(values = colors)
plot(p.frate)

rm(list=ls())

#········································································ STEP 4
## Q2: Are the scores influenced by NA in genes? 
## Run first 'Step 0', Set Environment, 7, 10, 11, 14

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


## Select a subset of 50 samples (to reduce computational time)
sampleRd<-colnames(SLE.prec)[sample(1:ncol(SLE.prec),50)]

## Subsets
subsets<-rep(seq(0.1, 0.9, by = 0.1),each = 5)
subsets<-subsets[length(subsets):1]


Results.NA<-as.data.frame(do.call("rbind",lapply(1:length(subsets),function(i){
  
  print(subsets[i])
  #' Select random patients within each iteration and cohort size)
  #' Remove different percentage of genes 
  tmp.SLE<-SLE.prec[-sample(1:nrow(SLE.prec),trunc(nrow(SLE.prec)*subsets[i])),sampleRd]
  
  
  res<-as.data.frame(do.call("rbind",lapply(methods,function(method){
  
    print(method)
    ## Select 150 genesets with at least 3 shared genes with tmp.SLE 
    ## To reduce computational time and to keep alwways the same length to compare results
    tmp.geneset<-geneset.DB[rownames(PREC[[method]])]
    sel.geneset<-unlist(lapply(tmp.geneset,function(x){
      if(length(intersect(rownames(tmp.SLE),x))>=3){
        return(TRUE)
      }else{
        return(FALSE)
      }
    })) #; print(table(sel.geneset))
    tmp.geneset<-tmp.geneset[sel.geneset]
    
    if(method=="MLM"){
      net <- do.call("rbind", lapply(seq_len(length(tmp.geneset)), function (x) {
        res <- data.frame("source"=rep(names(tmp.geneset)[[x]],
                                       length(tmp.geneset[[x]])),
                          "target"=as.character(tmp.geneset[[x]]),
                          "weight"=rep(1,length(tmp.geneset[[x]])),
                          "mor"=rep(1,length(tmp.geneset[[x]])))
        return(res)
      }))
      
      co.lin <- as.data.frame(decoupleR::check_corr(net,
                                                    .source="source",
                                                    .target = "target",
                                                    .mor = "mor"))
      net <- net[!net$source %in% as.character(co.lin[
        co.lin$correlation > 0.4,"source"]),]
      
      tmp.geneset<-tmp.geneset[names(tmp.geneset) %in% co.lin$source]
    }
    print(length(tmp.geneset))
    tmp.geneset<-tmp.geneset[sample(1:length(tmp.geneset),150)]
    
      
    #print(method)
    ncores<-ifelse(method %in% c("M-scores", "MDT", "UDT", "ULM"), 12, 1)
    
    labs <- rep("SLE",ncol(tmp.SLE))
    if(method=="M-score"){
      tmp.SLE.m<-cbind(tmp.SLE,HC.prec[rownames(tmp.SLE),])
      labs<-c(labs,rep("Healthy",ncol(HC.prec)))
    }else{
      tmp.SLE.m<-tmp.SLE
    }
    
    tmp.sc<-getScores(inputData=tmp.SLE.m, geneSets=tmp.geneset, method=method,
                      labels=labs, cores=ncores, minSize=3,
                      times = ifelse(grepl("norm",method),100, ifelse(grepl("corr",method),200,2)))
    
    ## Check valid paths, get Loss pathways
    lossPaths<-rownames(tmp.sc)[(apply(tmp.sc, 1, function(row) any(is.na(row) | is.infinite(row)))) | 
                           (apply(tmp.sc, 1, function(row) all(row == 0, na.rm = TRUE)))]
    tmp.sc<-tmp.sc[!rownames(tmp.sc) %in% lossPaths,]
    
    sharedPaths<-intersect(rownames(tmp.sc),
                           rownames(PREC[[method]][names(tmp.geneset),]))
    
    tmp.sc<-tmp.sc[sharedPaths,]
    ref<-PREC[[method]][sharedPaths,colnames(tmp.sc)]
    
    ## Measure pearson and euclidean distance
    sim.vals<-as.data.frame(do.call("rbind",lapply(1:ncol(tmp.sc),function(pat){
      x<-tmp.sc[,pat]
      y<-ref[,pat]
      df<-cbind(x,y)
      #df <- df[apply(df, 1, function(row) all(is.finite(row))), ]
      sim.values.p<-vectSimilarity(x=df[,1],y=df[,2],method = c("euclidean","pearson"))
      sim.values.p[is.na(sim.values.p)]<-0
      return(sim.values.p)
    })))
    sim.vals<-apply(sim.vals,2,function(col){mean(col,na.rm=T)})
    
    return(c(sim.vals,"lossPathways"=length(lossPaths)))
  })))
  
  res<-data.frame("method"=methods,"size"=subsets[i],"pearson"=res$pearson,
                  "euclidean"=res$euclidean,"lossPathways"=res$lossPathways)
  print(res)
  return(res)
})))

save.image((paste0(getwd(),"/RData/results_NA.RData")))

## Plot

df<-Results.NA
df$size<-(1-df$size)*100

df <- df %>%
  group_by(size, method) %>%
  summarise(
    mean_pearson = mean(pearson),
    mean_euclidean = mean(euclidean),
    .groups = "drop"
  )


p1<-ggplot(df,aes(x=size,y=mean_pearson,color=method))+
  geom_point(size=1.2) +
  geom_line(aes(group=method),linewidth=0.75) +
  scale_x_reverse()+
  labs(x="Gene retention (%)",
       y="Pearson corelation",
       title = "")+ theme_bw()+
  theme(legend.text = element_text(size=8),
        axis.text.x=element_text(size=7, color = "black"),
        axis.text.y = element_text(size=7, color = "black"),
        axis.title = element_text(size=8, color = "black", face = "bold"), 
        strip.text = element_text(size=7), 
        plot.title = element_text(size=10,hjust=0.5,face="bold"))+
  scale_color_manual(values=colors) +
  scale_x_continuous(breaks = c(90,80,70,60,50,40,30,20,10),
                     labels = c(90,80,70,60,50,40,30,20,10))
plot(p1)

p2<-ggplot(df,aes(x=size,y=mean_euclidean,color=method))+
  geom_point(size=1.2) +
  geom_line(aes(group=method),linewidth=0.75) +
  scale_x_reverse()+
  labs(x="Gene retention (%)",
       y="Euclidean distance",
       title = "")+ theme_bw()+
  theme(legend.text = element_text(size=8),
        axis.text.x=element_text(size=7, color = "black"),
        axis.text.y = element_text(size=7, color = "black"),
        axis.title = element_text(size=8, color = "black", face = "bold"), 
        strip.text = element_text(size=7), 
        plot.title = element_text(size=10,hjust=0.5,face="bold"))+
  scale_color_manual(values=colors) +
  scale_x_continuous(breaks = c(90,80,70,60,50,40,30,20,10),
                     labels = c(90,80,70,60,50,40,30,20,10))
plot(p2)


ggarrange(p1,p2,ncol = 2,nrow = 1,common.legend = T,legend = "bottom")


