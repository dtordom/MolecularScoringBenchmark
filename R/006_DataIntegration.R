################################################################################
## Personalized molecular profile scoring benchmark
## R version 4.3.1 (2024-04-10)
## 
################################################################################
#' Data integration 


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
                 "reshape2","qvalue","ggbreak","pheatmap","pathMED","scales","dplyr",
                 "cluster"))


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


##······································································· Step 0
## Check data integration

DATA<-list()

MET<-rbind(data.frame("SampleID"=c(rownames(clin.ln1),colnames(HC.ln1)),
                      "Diagnosis"=c(rep("SLE",ncol(SLE.ln1)),rep("NHV",ncol(HC.ln1))),
                      "pLN"=c(clin.ln1$group,rep("NHV",ncol(HC.ln1))),
                      "Dataset"=rep("ln1",ncol(SLE.ln1)+ncol(HC.ln1))),
           
           data.frame("SampleID"=c(rownames(clin.ln2),colnames(HC.ln2)),
                      "Diagnosis"=c(rep("SLE",ncol(SLE.ln2)),rep("NHV",ncol(HC.ln2))),
                      "pLN"=c(clin.ln2$group,rep("NHV",ncol(HC.ln2))),
                      "Dataset"=rep("ln2",ncol(SLE.ln2)+ncol(HC.ln2))),
           
           data.frame("SampleID"=c(rownames(clin.ln3),colnames(HC.ln3)),
                      "Diagnosis"=c(rep("SLE",ncol(SLE.ln3)),rep("NHV",ncol(HC.ln3))),
                      "pLN"=c(clin.ln3$group,rep("NHV",ncol(HC.ln3))),
                      "Dataset"=rep("ln3",ncol(SLE.ln3)+ncol(HC.ln3))),
           
           data.frame("SampleID"=c(rownames(clin.ln4),colnames(HC.ln4)),
                      "Diagnosis"=c(rep("SLE",ncol(SLE.ln4)),rep("NHV",ncol(HC.ln4))),
                      "pLN"=c(clin.ln4$group,rep("NHV",ncol(HC.ln4))),
                      "Dataset"=rep("ln4",ncol(SLE.ln4)+ncol(HC.ln4)))
           )

## Gene expression
sharedGenes <- Reduce(intersect, lapply(list(SLE.ln1,SLE.ln2,SLE.ln3,SLE.ln4), rownames))

DATA[["geneExp"]]<-cbind(SLE.ln1[sharedGenes,],HC.ln1[sharedGenes,],
                         SLE.ln2[sharedGenes,],HC.ln2[sharedGenes,],
                         SLE.ln3[sharedGenes,],HC.ln3[sharedGenes,],
                         SLE.ln4[sharedGenes,],HC.ln4[sharedGenes,])

## Gene expression normalized by zscore by samples
DATA[["geneExp_norm"]] <- apply(DATA$geneExp, 2, function(x) (x - mean(x)) / sd(x))

## Score matrices
for(m in methods){
  sharedPaths <- Reduce(intersect, lapply(list(SCORES$ln1[[m]],
                                               SCORES$ln2[[m]],
                                               SCORES$ln3[[m]],
                                               SCORES$ln4[[m]]), rownames))
  
  DATA[[m]]<-cbind(SCORES$ln1[[m]][sharedPaths,c(colnames(SLE.ln1),colnames(HC.ln1))],
                   SCORES$ln2[[m]][sharedPaths,c(colnames(SLE.ln2),colnames(HC.ln2))],
                   SCORES$ln3[[m]][sharedPaths,c(colnames(SLE.ln3),colnames(HC.ln3))],
                   SCORES$ln4[[m]][sharedPaths,c(colnames(SLE.ln4),colnames(HC.ln4))] )
}




RESULTS.BE<-lapply(1:length(DATA),function(d){
  
  sb<-NULL
  print(names(DATA)[d])
  ## Silouette pLN
  tmp.clin<-MET[MET$pLN!="NHV",]
  tmp.dat<-DATA[[d]][,tmp.clin$SampleID]
  tmp.dat <- tmp.dat[apply(tmp.dat, 1, function(row) all(is.finite(row))), ]
  
  x<-prcomp(t(tmp.dat))
  D<-as.data.frame(x$x[,1:2])
  D<-cbind(tmp.clin,D)
  
  dist_matrix<-dist(D[,grepl("PC",colnames(D))])
  sil_batch <- silhouette(ifelse(D$pLN=="pLN",2,ifelse(D$pLN=="NoLN",1,0)),
                          dist = dist_matrix)
  sb<-c(sb,summary(sil_batch)$si.summary["Median"])
  
  
  ## Silouette Dataset and Diagnosis
  tmp.clin<-MET
  tmp.dat<-DATA[[d]][,tmp.clin$SampleID]
  tmp.dat <- tmp.dat[apply(tmp.dat, 1, function(row) all(is.finite(row))), ]
  
  x<-prcomp(t(tmp.dat))
  D<-as.data.frame(x$x[,1:2])
  D<-cbind(tmp.clin,D)
  
  dist_matrix<-dist(D[,grepl("PC",colnames(D))])
  
  sil_batch <- silhouette(ifelse(D$pLN=="pLN",2,ifelse(D$pLN=="NoLN",1,0)),
                          dist = dist_matrix)
  sb<-c(sb,summary(sil_batch)$si.summary["Median"])
  
  sil_batch <- silhouette(ifelse(D$Dataset=="ln1",1,ifelse(D$Dataset=="ln2",2,
                                                           ifelse(D$Dataset=="ln3",3,4))), dist = dist_matrix)
  sb<-c(sb,summary(sil_batch)$si.summary["Median"])
  names(sb)<-c("Nephritis","Diagnosis","Dataset"); print(sb)
  
  
  ## Remove outliers
  outliers <- apply(scale(D[,c("PC1","PC2")]), 1, function(row) any(abs(row) > 3))
  print(table(outliers))
  D<-D[!outliers,]
  if(names(DATA)[d]=="norm_WMEAN" | names(DATA)[d]=="norm_WSUM"){
    D<-D[D$PC1<0 & D$PC2>-6000,]
  }
  
  p1<-ggplot(D,aes(x=PC1,y=PC2,fill=Dataset)) + theme_classic() + geom_point(shape=21,size=0.8,color="black",stroke=0.1)+
    theme(legend.text = element_text(size=7),
          axis.text.x=element_text(size=7, color = "black",hjust=1, vjust=0.5,angle = 90), #,angle = 90
          axis.text.y = element_text(size=7, color = "black"),
          axis.title = element_text(size=7, color = "black", face = "bold"), 
          strip.text = element_text(size=7,face = "bold"), 
          strip.background = element_blank(),
          plot.title = element_text(size=9,face="bold"))+
    theme(legend.position = "none")+
    scale_fill_manual(values=c("ln1"="#cd5554","ln2"="#0198b1","ln3"="#8fbc8f","ln4"="orange1"))+labs(title=names(DATA)[d])
  plot(p1)
  
  p2<-ggplot(D,aes(x=PC1,y=PC2,fill=pLN)) + theme_classic() + geom_point(shape=21,size=0.8,color="black",stroke=0.1)+
    theme(legend.text = element_text(size=7),
          axis.text.x=element_text(size=7, color = "black",hjust=1, vjust=0.5,angle = 90), #,angle = 90
          axis.text.y = element_text(size=7, color = "black"),
          axis.title = element_text(size=7, color = "black", face = "bold"), 
          strip.text = element_text(size=7,face = "bold"), 
          strip.background = element_blank(),
          plot.title = element_text(size=9,face="bold"))+
    theme(legend.position = "none")+
    scale_fill_manual(values=c("pLN"="#cd5554","NoLN"="grey","NHV"="#8fbc8f"))+labs(title=names(DATA)[d])
  plot(p2)
  
  res<-list("SB"=sb,"p.dataset"=p1,"p.pln"=p2)
  return(res)
})
names(RESULTS.BE)<-names(DATA)


ggarrange(RESULTS.BE[[1]]$p.dataset,RESULTS.BE[[2]]$p.dataset,
          RESULTS.BE[[3]]$p.dataset,RESULTS.BE[[4]]$p.dataset,
          RESULTS.BE[[5]]$p.dataset,RESULTS.BE[[6]]$p.dataset,
          RESULTS.BE[[7]]$p.dataset,RESULTS.BE[[8]]$p.dataset,
          RESULTS.BE[[9]]$p.dataset,RESULTS.BE[[10]]$p.dataset,
          RESULTS.BE[[11]]$p.dataset,RESULTS.BE[[12]]$p.dataset,
          RESULTS.BE[[13]]$p.dataset,RESULTS.BE[[14]]$p.dataset,
          RESULTS.BE[[15]]$p.dataset,RESULTS.BE[[16]]$p.dataset,
          RESULTS.BE[[17]]$p.dataset,RESULTS.BE[[18]]$p.dataset,
          RESULTS.BE[[19]]$p.dataset,RESULTS.BE[[20]]$p.dataset,
          ncol=5,nrow=4,common.legend = T,legend = "bottom")


ggarrange(RESULTS.BE[[1]]$p.pln,RESULTS.BE[[2]]$p.pln,
          RESULTS.BE[[3]]$p.pln,RESULTS.BE[[4]]$p.pln,
          RESULTS.BE[[5]]$p.pln,RESULTS.BE[[6]]$p.pln,
          RESULTS.BE[[7]]$p.pln,RESULTS.BE[[8]]$p.pln,
          RESULTS.BE[[9]]$p.pln,RESULTS.BE[[10]]$p.pln,
          RESULTS.BE[[11]]$p.pln,RESULTS.BE[[12]]$p.pln,
          RESULTS.BE[[13]]$p.pln,RESULTS.BE[[14]]$p.pln,
          RESULTS.BE[[15]]$p.pln,RESULTS.BE[[16]]$p.pln,
          RESULTS.BE[[17]]$p.pln,RESULTS.BE[[18]]$p.pln,
          RESULTS.BE[[19]]$p.pln,RESULTS.BE[[20]]$p.pln,
          ncol=5,nrow=4,common.legend = T,legend = "bottom")


## Plot Silhouette scores
m<-do.call("rbind",lapply(RESULTS.BE,function(x){
  x$SB
}))
m<-reshape::melt(m)
colnames(m)<-c("method","var","value")

m<-m[m$var!="Diagnosis",]

line_positions <- seq(1.5, length(levels(factor(m$method))) - 0.5, by = 1)

ggplot(m, aes(x = value, 
              y = factor(method, levels = c(rev(methods), "geneExp_norm", "geneExp")), 
              fill = var)) +
  geom_col(stat = "identity", 
           position = position_dodge(width = 0.7), 
           color = "black", 
           width = 0.6) +  # Cambiar ancho de las barras si es necesario
  theme_classic() +  
  labs(x = "Silhouette", y = "", fill = "") +
  geom_vline(xintercept = 0, color = "black") +  
  geom_hline(aes(yintercept = as.numeric(method)), color = "gray", linetype ="dashed", size = 0.5) +
  #geom_hline(yintercept = line_positions, color = "gray2", size = 0.5,alpha=0.4) +
  theme(
    legend.text = element_text(size = 7),
    axis.text.x = element_text(size = 7, color = "black"), # ,angle = 90
    axis.text.y = element_text(size = 7, color = "black", face = "bold"),
    axis.title = element_text(size = 7, color = "black", face = "bold"), 
    strip.text = element_text(size = 7, face = "bold"), 
    strip.background = element_blank(),
    plot.title = element_text(size = 9, face = "bold"),
    legend.position = "bottom"
  ) +
  scale_fill_manual(values = c("Nephritis" = "steelblue4", "Dataset" = "grey"))







