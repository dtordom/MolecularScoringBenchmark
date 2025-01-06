################################################################################
## Personalized molecular profile scoring benchmark
## R version 4.3.1 (2024-04-10)
## 
################################################################################
## Testing Computational costs based on cohort size and pathway database size


# devtools::install_github("jordimartorell/pathMED")
##······································································· Step 0
## Set environment

setwd("D:/Work/BMS/PMPS/PMSPbenchmark")
set.seed(1234)

source(paste0(getwd(),"/code/utils.R"))

check.packages(c("vctrs","stats","stringr","caret","singscore","GSVA","parallel","lsa",
                 "metrica","BiocParallel","pbapply","reshape","ggplot2","ggpubr","psych",
                 "jaccard","NbClust","ConsensusClusterPlus","mclust","decoupleR",
                 "reshape2","qvalue","ggbreak","pheatmap","pathMED","scales"))


load("D:/Work/BMS/PMPS/PMSPbenchmark/RData/SCORES.RData")

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
  

##······································································· Step 1
## Computational time vs Sample size and length of the database ----

db.length<-c(1,3,5,10,15,25,50,75,100,250,500,750,1000,1250,1500)
cohort.size<-c(2,3,5,10,15,25,50,75,100,150,250,300,400,500,600)
  
## Get times for different input (sample size, database length)
Results.cp1<-as.data.frame(do.call("rbind",lapply(1:length(db.length),function(i){

  ## Prepare data
  if(cohort.size[i]>ncol(SLE.prec)){
    data.tmp<-data.frame(cbind(SLE.prec,SLE.prec))
  }else{
    data.tmp<-SLE.prec
  }
  
  genesetDB.rd<-geneset.DB[sample(1:length(geneset.DB),db.length[i])]
  genesetDB.200<-geneset.DB[1:100]
  data.sub.rd<-data.tmp[,sample(1:ncol(data.tmp),cohort.size[i]),drop=FALSE]
  data.sub.50<-SLE.prec[,1:50,drop=FALSE]
  
  
  times<-as.data.frame(do.call("rbind",lapply(methods,function(method){
    
    ncores<-ifelse(method %in% c("M-scores", "MDT", "UDT", "ULM"), 12, 1)
    
    if(method=="M-score"){ ## Get Mscores also for HC
      data.sub.rd.m<-cbind(data.sub.rd,HC.prec)
      labs.sub.rd.m<-c(rep("SLE",ncol(data.sub.rd)),rep("Healthy",ncol(HC.prec)))
      
      data.sub.50.m<-cbind(data.sub.50,HC.prec)
      labs.sub.50.m<-c(rep("SLE",ncol(data.sub.50)),rep("Healthy",ncol(HC.prec)))
    }else{
      data.sub.rd.m<-data.sub.rd
      data.sub.50.m<-data.sub.50
    }
    
    ## Pathways
    startTime <- Sys.time()
    tmp<-getScores(inputData=data.sub.50.m, geneSets=genesetDB.rd, method = method,
                   labels = labs.sub.50.m, cores = ncores, minSize=3, 
                   times=ifelse(grepl("norm",method),100,ifelse(grepl("corr",method),200,2)))
    finalTime<-as.numeric(difftime(Sys.time(),startTime,units = "min")) 
    gc()
    
    ## Patients
    startTime <- Sys.time()
    tmp<-getScores(inputData=data.sub.rd.m,geneSets=genesetDB.200, method = method,
                   labels = labs.sub.rd.m, cores = ncores, minSize=3, 
                   times=ifelse(grepl("norm",method),100,ifelse(grepl("corr",method),200,2)))
    finalTime<-c(finalTime,as.numeric(difftime(Sys.time(),startTime,units = "min")))
    names(finalTime)<-c("time.db","time.size")
    gc()
    return(finalTime)
  })))
  
  times$cohort.size<-cohort.size[i]
  times$length.db<-db.length[i]
  times$method<-methods
  print(times)
  return(times)
})))

save.image("D:/Work/BMS/PMPS/PMSPbenchmark/RData/CompTime1_1c.RData")


##······································································· Step 2
## Computational time comparing length of paths database elements ----

genes<-intersect(unique(unlist(geneset.DB)),rownames(SLE.prec))
sizes<-c(5,8,10,15,20,30,40,50,75,100,150,200,300,400,500,600)

Results.cp2<-as.data.frame(do.call("rbind",lapply(sizes,function(sz){
  
  geneset.rd<-lapply(1:50,function(sz.i){ tmp<-sample(genes,sz)})
  names(geneset.rd)<-paste0("Random.",1:50)
  data.sub.50<-SLE.prec[,1:50]
  
  times<-as.numeric(lapply(methods,function(method){
    
    ncores<-ifelse(method %in% c("M-scores", "MDT", "UDT", "ULM"), 12, 1)
    
    if(method=="M-score"){ ## Get Mscores also for HC
      data.sub.50.m<-cbind(data.sub.50,HC.prec)
      labs.sub.50.m<-c(rep("SLE",ncol(data.sub.50)),rep("Healthy",ncol(HC.prec)))
    }else{
      data.sub.50.m<-data.sub.50
    }
    
    startTime <- Sys.time()
    tmp<-getScores(inputData=data.sub.50.m, geneSets=geneset.rd, method = method,
                   labels = labs.sub.50.m, cores = ncores, minSize=3, 
                   times=ifelse(grepl("norm",method),100,ifelse(grepl("corr",method),200,2)))
    finalTime<-as.numeric(difftime(Sys.time(),startTime,units = "min")) 
    gc()
    return(finalTime)
  }))
  
  res<-data.frame("method"=methods,"pathsLength.time"=times,"pathsLength"=sz)
  return(res)
})))

save.image("D:/Work/BMS/PMPS/PMSPbenchmark/RData/CompTime1_1c.RData")


##······································································· Step 3
## Get Plots ----
# load("D:/Work/BMS/PMPS/PMSPbenchmark/RData/CompTime1_1c.RData")

custom_trans <- trans_new( name = "custom", 
  transform = function(x) ifelse(x <= 100, x / 200, 0.5 + (x - 100) / 1000),
  inverse = function(x) ifelse(x <= 0.5, x * 200, 100 + (x - 0.5) * 1000))

p1<-ggplot(Results.cp1,aes(x=cohort.size,y=time.size,color=method))+
  geom_point(size=1.75) +
  geom_line(aes(group=method),linewidth=1) +
  labs(x="Cohort size",y="Time (min)",title = "")+ theme_bw()+
  theme(legend.text = element_text(size=8),
        axis.text.x=element_text(size=7, color = "black"),
        axis.text.y = element_text(size=7, color = "black"),
        axis.title = element_text(size=8, color = "black", face = "bold"), 
        strip.text = element_text(size=7), 
        plot.title = element_text(size=10,hjust=0.5,face="bold"))+
  scale_color_manual(values=colors)+
  scale_x_continuous(trans = custom_trans,
    breaks = c(0, 50, 100, 200, 400, 600),
    labels = c("0", "50", "100", "200", "400", "600") )
plot(p1)


custom_trans <- trans_new(name = "custom", 
  transform = function(x) ifelse(x <= 200, x / 400, 0.5 + (x - 200) / 2600),
  inverse = function(x) ifelse(x <= 0.5, x * 400, 200 + (x - 0.5) * 2600))

p2<-ggplot(Results.cp1,aes(x=length.db,y=time.db,color=method))+
  geom_point(size=1.75) +
  geom_line(aes(group=method),linewidth=1) +
  labs(x="Number of genesets",y="Time (min)",title = "")+ theme_bw()+
  theme(legend.text = element_text(size=8),
        axis.text.x=element_text(size=7, color = "black"),
        axis.text.y = element_text(size=7, color = "black"),
        axis.title = element_text(size=8, color = "black", face = "bold"), 
        strip.text = element_text(size=7), 
        plot.title = element_text(size=10,hjust=0.5,face="bold"))+
  scale_color_manual(values=colors)+
  scale_x_continuous(trans = custom_trans,
                     breaks = c(0, 100, 200, 500, 1000, 1500),
                     labels = c("0", "100", "200", "500", "1000", "1500") )
plot(p2)


custom_trans <- trans_new( name = "custom", 
                           transform = function(x) ifelse(x <= 100, x / 200, 0.5 + (x - 100) / 1000),
                           inverse = function(x) ifelse(x <= 0.5, x * 200, 100 + (x - 0.5) * 1000))

p3<-ggplot(Results.cp2,aes(x=pathsLength,y=pathsLength.time,color=method))+
  geom_point(size=1.75) +
  geom_line(aes(group=method),linewidth=1) +
  labs(x="Mean length of genesets",y="Time (min)",title = "")+ theme_bw()+
  theme(legend.text = element_text(size=8),
        axis.text.x=element_text(size=7, color = "black"),
        axis.text.y = element_text(size=7, color = "black"),
        axis.title = element_text(size=8, color = "black", face = "bold"), 
        strip.text = element_text(size=7), 
        plot.title = element_text(size=10,hjust=0.5,face="bold"))+
  scale_color_manual(values=colors)+
  scale_x_continuous(trans = custom_trans,
                     breaks = c(0, 50, 100, 200, 400, 600),
                     labels = c("0", "50", "100", "200", "400", "600") )
plot(p3)


ggarrange(p1,p2,p3,ncol=3,nrow=1,common.legend = T,legend = "bottom",
          labels = c("A)", "B)", "C)"),
          font.label = list(size = 12, color = "black", face = "bold"))








