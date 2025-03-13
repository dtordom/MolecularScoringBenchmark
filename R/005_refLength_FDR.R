################################################################################
## Personalized molecular profile scoring benchmark
## R version 4.3.1 (2024-04-10)
## 
################################################################################
#' Testing impact of number of genesets in the reference on score reproducibility
#' and false rate discovery across methods


# devtools::install_github("jordimartorell/pathMED")
##······································································· Step 0
## Set environment

setwd("C:/Users/danie/Desktop/WORK/BENCHMARK_25")
set.seed(1234)

load(paste0(getwd(),"/RData/SCORES.RData"))

source(paste0(getwd(),"/code/utils.R"))

check.packages(c("vctrs","stats","stringr","caret","singscore","GSVA","parallel","lsa",
                 "metrica","BiocParallel","pbapply","reshape","ggplot2","ggpubr","psych",
                 "jaccard","NbClust","ConsensusClusterPlus","mclust","decoupleR",
                 "reshape2","qvalue","ggbreak","pheatmap","scales","dplyr"))

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

##······································································· Step 1
## Q3: Are the scores influenced by the number of genesets tested? 

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

## Different lengths for geneset database
sizes<-c(1,2,3,4,5,10,20,30,50,75,100,150,200,400,750,1000)
sizes.list <- replicate(10, sizes, simplify = FALSE)

## Select 50 random samples (to reduce computational time)
tmp.SLE<-SLE.prec[,sample(1:ncol(SLE.prec),50),drop=FALSE]


Results.gl<-as.data.frame(do.call("rbind",lapply(methods,function(method){
  
  print(method)
  res.perm<-as.data.frame(do.call("rbind",lapply(sizes.list,function(perm){
    
    ## Suffle genesets
    tmp.geneset<-geneset.DB[rownames(PREC[[method]])]
    tmp.geneset<-tmp.geneset[sample(1:length(tmp.geneset),length(tmp.geneset))]
    
    res.met<-do.call("rbind",lapply(perm,function(i){
      
      ## Select 'i' random genesets
      tmp.geneset.i<-tmp.geneset[1:i]
      
      ncores<-ifelse(method %in% c("M-Scores", "MDT", "UDT", "ULM"), detectCores()-2, 1)
      labs <- rep("SLE",ncol(tmp.SLE))
      if(method=="M-Scores"){ ## HC is also needed for M-Scores
        dat<-cbind(tmp.SLE,HC.prec)
        labs<-c(labs,rep("Healthy",ncol(HC.prec)))
      }else{
        dat<-tmp.SLE
      }
      
      tmp.sc<-getScores(inputData=dat, geneSets=tmp.geneset.i, method=method,
                        labels=labs, returnHC = FALSE, cores=ncores, minSize=3,
                        times = ifelse(grepl("norm",method),100,2))
      
      return(as.numeric(tmp.sc[names(tmp.geneset.i)[1],]))
    }))
    
    combs<-as.data.frame(t(combn(1:nrow(res.met),2)))
    
    sim.vals<-as.data.frame(do.call("rbind",lapply(1:nrow(combs),function(p){
      return(vectSimilarity(as.numeric(res.met[combs[p,1],]),
                            as.numeric(res.met[combs[p,2],]),
                            method = c("euclidean","pearson"))) 
    })))
    sim.vals.mean<-unlist(lapply(1:ncol(sim.vals),function(p){
      mean(as.numeric(sim.vals[,p]),na.rm = T)}))
    names(sim.vals.mean)<-colnames(sim.vals)
    
    return(sim.vals.mean)
  }))) # perm
  
  res.perm$method<-method

  return(res.perm)
})))

save.image(paste0(getwd(),"/RData/Results_GL.RData"))


## Plot
df<-Results.gl

custom_trans <- trans_new(
  name = "custom",
  transform = function(x) ifelse(x <= 30, x * 2 / 3, 20 + (x - 30) / 470 * (1 / 3) * 80),
  inverse = function(x) ifelse(x <= 20, x * 3 / 2, 30 + (x - 20) * 470 / (1 / 3 * 80))
)


p1<-ggplot(df,aes(x=method,y=euclidean,fill=method))+theme_classic()+
  geom_jitter(size=1,alpha=0.8,shape=21,color="black",stroke=0.2)+
  theme_classic() +
  theme(legend.text = element_text(size=8),
        axis.text.x=element_text(size=7, color = "black",hjust=1, vjust=0.5,angle = 90), #,angle = 90
        axis.text.y = element_text(size=7, color = "black"),
        axis.title = element_text(size=8, color = "black", face = "bold"), 
        strip.text = element_text(size=8,face = "bold"), 
        #strip.background = element_rect(fill = "white", color = "black", size = 0.5),  
        strip.background = element_blank(),
        plot.title = element_text(size=10,hjust=0.5,face="bold"))+
  theme(legend.position = "none")+
  labs(x = "", y = "Euclidean distance") +
  scale_y_continuous(
    trans = custom_trans,
    breaks = c(0, 10, 20, 30, 200, 400, 600),  # Posiciones de las etiquetas
    labels = c(0, 10, 20, 30, 200, 400, 600)   # Textos de las etiquetas
  )+
  geom_boxplot(lwd=0.2,alpha=0.8,outlier.colour = NA,width = 0.8,alpha=0.8)+ 
  scale_fill_manual(values=colors)+scale_x_discrete(limits=methods)
plot(p1)


p2<-ggplot(df,aes(x=method,y=pearson,fill=method))+theme_classic()+
  geom_jitter(size=1,alpha=0.8,shape=21,color="black",stroke=0.2)+
  theme_classic() +
  theme(legend.text = element_text(size=8),
        axis.text.x=element_text(size=7, color = "black",hjust=1, vjust=0.5,angle = 90), #,angle = 90
        axis.text.y = element_text(size=7, color = "black"),
        axis.title = element_text(size=8, color = "black", face = "bold"), 
        strip.text = element_text(size=8,face = "bold"), 
        #strip.background = element_rect(fill = "white", color = "black", size = 0.5),  
        strip.background = element_blank(),
        plot.title = element_text(size=10,hjust=0.5,face="bold"))+
  theme(legend.position = "none")+
  labs(x = "", y = "Pearson correlation") +
  ylim(0.4,1)+
  geom_boxplot(lwd=0.2,alpha=0.8,outlier.colour = NA,width = 0.8,alpha=0.8)+ 
  scale_fill_manual(values=colors)+scale_x_discrete(limits=methods)
plot(p2)

ggarrange(p2,p1,ncol = 2,nrow = 1,common.legend = T,legend = "bottom")



##······································································· Step 1
## Q4: FDR
# Run Set Environment first

## Create a reference of random genesets with different lengths
genes<-unique(unlist(geneset.DB))
setSizes<-c(5,10,25,50,100)

geneset.RD<-lapply(setSizes,function(sz){
  gs.x<-lapply(1:1000,function(i){
    return(sample(genes,sz))
  })
  names(gs.x)<-paste0("rd.",1:1000,"_",sz)
  return(gs.x)
})


## subset of 100 pats to reduce computational size
sle.tmp<-SLE.prec[,sample(1:ncol(SLE.prec),100)]
hc.tmp<-HC.prec[,sample(1:ncol(HC.prec),100)]

## Get scores
RD.scores<-lapply(methods,function(method){
  cat(paste0("\nRuning ",method,"\n"))
  
  dat<-cbind(sle.tmp,hc.tmp)
  labs<-c(rep("SLE",ncol(sle.tmp)),
          rep("Healthy",ncol(hc.tmp)))
  
  ncores<-ifelse(method %in% c("M-Scores", "MDT", "UDT", "ULM"), detectCores()-2, 1)

  tmp.scores<-do.call("rbind",lapply(geneset.RD,function(tmp.geneset){
    tmp.sc<-getScores(inputData=dat, returnHC = TRUE, geneSets=tmp.geneset,
                      method = method, labels = labs,cores = ncores, minSize=3,
                     times =ifelse(grepl("norm",method),100,2))
    return(tmp.sc)
  }))
  
 
  return(tmp.scores)
})
names(RD.scores)<-methods


save.image(paste0(getwd(),"/RData/Results_fdr.RData"))


## Get significance ·················

PL<-lapply(1:length(RD.scores),function(i){
  
  print(names(RD.scores)[i])
  dat<-RD.scores[[i]]
  
  ## Fix Inf values
  dat[dat == Inf] <- NA
  
  nhv<-dat[,colnames(dat) %in% colnames(HC.prec)]
  nhv<-data.frame("mean"=apply(nhv,1,function(r){mean(r,na.rm=T)}),
                  "sd"=apply(nhv,1,function(r){sd(r,na.rm=T)}))
  
  DAT<-dat[,colnames(dat) %in% colnames(SLE.prec)]
  
  if(names(RD.scores)[i]!="M-Scores"){ 
    ## scale with respect to Healthy controls
    tmp<-do.call("cbind",lapply(1:ncol(DAT),function(x){
      res <- (DAT[,x] - nhv[, "mean"]) / nhv[, "sd"]
      return(res)}))
    colnames(tmp)<-colnames(DAT)
    rownames(tmp)<-rownames(DAT)
    DAT<-tmp
  }
  
  ## Fix Inf values and 0 values
  DAT[DAT == Inf | is.na(DAT)] <- 0
  
  SIGN<-unlist(lapply(1:nrow(DAT),function(row){
    return(sum(abs(DAT[row,]) > 1.65))
  }))
  #plot(density(SIGN),main = names(RD.scores)[i])
  
  df<-data.frame("fdr"=SIGN,"ID"=paste0("gs",sub(".*_", "", rownames(DAT))))
  
  p1<-ggplot(df, aes(x = fdr, color = ID)) +
    geom_density(linewidth=1) + xlim(0,60)+
    theme(legend.text = element_text(size=7),
          axis.text.x=element_text(size=7, color = "black"), #,angle = 90
          axis.text.y = element_text(size=7, color = "black"),
          axis.title = element_blank(), 
          strip.text = element_text(size=7,face = "bold"), 
          #strip.background = element_rect(fill = "white", color = "black", size = 0.5),  
          strip.background = element_blank(),
          plot.title = element_text(size=7,hjust=0.5,face="bold"))+
    labs(title = names(RD.scores)[i],
         x = "False positive rate",
         y = "Density") + theme_bw() +
    scale_color_manual(values=c("gs100"="lightblue","gs50"="lightskyblue3","gs25"="lightskyblue4",
                                "gs10"="steelblue3","gs5"="steelblue4"))
  
  plot(p1)
  return(p1)
})
names(PL)<-names(RD.scores)


ggarrange(PL[[1]],PL[[2]],PL[[3]],PL[[4]],PL[[5]],PL[[6]],
          PL[[7]],PL[[8]],PL[[9]],PL[[10]],PL[[11]],PL[[12]],
          PL[[13]],PL[[14]],PL[[15]],PL[[16]],PL[[17]],PL[[18]],
          ncol = 6,nrow = 3,common.legend = TRUE,legend = "bottom")

