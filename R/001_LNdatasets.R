#'··············································································
#' Personalized molecular profile scoring benchmark
#' R version 4.3.1 (2024-04-10)
#' 
#'··············································································
#' Load datasets with LN patients and NHV

#' GSE65391: 30 NoLN, 89 pLN, 72 NHV
#' GSE99967: 13 NoLN, 23 pLN, 17 NHV
#' GSE72326: 36 NoLN, 14 pLN, 20 NHV
#' Petri Cohort (GSE45291): 135 NoLN, 12 pLN, 20 NHV

##······································································· Step 0
## Set environment ----


setwd("D:/Work/BMS/PMPS/PMSPbenchmark")
set.seed(1234)
source(paste0(getwd(),"/code/utils.R"))

pkgs<-c("parallel","matrixStats","biomaRt","NOISeq","stringr","dplyr","tidyr",
        "doParallel","caret","pbapply","BiocParallel","tibble","GEOquery","stringi")
check.packages(pkgs); rm(pkgs)


##······································································· Step 1
## Pascual et. al (GSE65391) ----

## gene expression for GSE65391 can be also downloaded from GEO using getGEO function
data<-read.table(file="D:/Work/Update_MyPROSLE/Datasets/GSE65391_expressionData.txt",
                 sep="\t",header=T,row.names = 1)

gset<-getGEO(GEO="GSE65391",GSEMatrix = TRUE)
clin<-phenoData(gset[[1]])
clin<-pData(clin)
data<-data[,rownames(clin)];

## Preprocessing (Log normalization and filtering non variable genes)
data<-norm.log(data)
nonVar.genes<-summ.var(data = data)
data<-data[!nonVar.genes$nzv,] ## remove genes with near-zero variance

clin<-clin[,c("subject:ch1","visit:ch1","disease state:ch1","nephritis_class:ch1",
              "biopsy_history:ch1","days_since_kidney_biopsy:ch1")]
colnames(clin)<-c("Patient","Visit","Diagnosis","NeprClass","NeprBX","TimeFromBX")


## Separate Healthy and SLE samples
HC.ln1<-data[,clin$Diagnosis=="Healthy"]
SLE<-data[,clin$Diagnosis=="SLE"]
clin<-clin[clin$Diagnosis=="SLE",]

## Select Nephritis and NonLN patients
clin<-clin[clin$TimeFromBX!="Data Not Available",]
clin$TimeFromBX<-as.numeric(clin$TimeFromBX)

clin<-clin[clin$NeprBX %in% c("Proliferative","NoLN"),]
table(clin$NeprBX)

clin<-clin[abs(clin$TimeFromBX)<=365 ,]
table(clin$NeprBX)

clin<-clin[clin$NeprBX=="NoLN" | (clin$NeprBX=="Proliferative" & grepl("Prolif",clin$NeprClass)),]
table(clin$NeprBX)

SLE.ln1<-SLE[,rownames(clin)]

## Remove genes with near-zero variance
nonVar.genes<-summ.var(data = HC.ln1); table(!nonVar.genes$nzv)
HC.ln1<-HC.ln1[!nonVar.genes$nzv,]
SLE.ln1<-SLE.ln1[!nonVar.genes$nzv,]

apply(SLE.ln1,1,sd)[order(apply(SLE.ln1,1,sd),decreasing=FALSE)][1:10]
apply(HC.ln1,1,sd)[order(apply(HC.ln1,1,sd),decreasing=FALSE)][1:10]

## PCA
x<-prcomp(t(SLE.ln1))
D<-as.data.frame(x$x[,1:2])
D<-cbind(clin,D)

ggplot(D,aes(x=PC1,y=PC2,fill=NeprBX)) + theme_classic() + geom_point(shape=21,size=3)+
  scale_fill_manual(values=c("#cd5554","#0198b1","#8fbc8f","gold"))

clin.ln1<-clin


rm(clin,data,SLE,nonVar.genes,x,D,gset)



##······································································· Step 2
## GSE99967 ----

gset <- getGEO("GSE99967", GSEMatrix = TRUE) ## Download. data from NCBI GEO
gset <- gset[[1]]

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
data <- temp

metadata <- phenoData(gset) ## Download clinical data from NCBI GEO
metadata <- pData(metadata)

## Select relevant variables for the analysis
metadata <- metadata[, c(
  "title", "disease state:ch1",
  "nephritis state:ch1", "biopsy class:ch1"
)]
colnames(metadata) <- c("title", "diagnosis", "LN", "biopsy")

## Separe Healthy controls and SLE samples
sel <- ifelse(metadata$diagnosis == "control", F, T)
HC.ln2 <- data[, rownames(metadata)[!sel]]
clin.ln2 <- metadata[sel, ]

clin.ln2$LN <- ifelse(clin.ln2$LN == "active LN (ALN)", "YES", "NO")
clin.ln2<-clin.ln2[clin.ln2$biopsy=="NA" | grepl("III",clin.ln2$biopsy) | grepl("IV",clin.ln2$biopsy),]

SLE.ln2 <- data[, rownames(clin.ln2)]

## Remove genes with near-zero variance
nonVar.genes<-summ.var(data = HC.ln2); table(!nonVar.genes$nzv)
HC.ln2<-HC.ln2[!nonVar.genes$nzv,]
SLE.ln2<-SLE.ln2[!nonVar.genes$nzv,]

apply(SLE.ln2,1,sd)[order(apply(SLE.ln2,1,sd),decreasing=FALSE)][1:10]
apply(HC.ln2,1,sd)[order(apply(HC.ln2,1,sd),decreasing=FALSE)][1:10]


## PCA
x<-prcomp(t(SLE.ln2))
D<-as.data.frame(x$x[,1:2])
D<-cbind(clin.ln2,D)

ggplot(D,aes(x=PC1,y=PC2,fill=LN)) + theme_classic() + geom_point(shape=21,size=3)+
  scale_fill_manual(values=c("#cd5554","#0198b1","#8fbc8f","gold"))


rm(data,D, genome, gset, metadata, nonVar.genes, temp, finalGenes, sel,x)



##······································································· Step 3
## GSE72326 ----

gset <- getGEO("GSE72326", GSEMatrix = TRUE) ## Download data from NCBI GEO
gset <- gset[[1]]

data <- exprs(gset) ## Get gene-expression matrix
data <- norm.log(data) ## Log2 transformation

nonVar.genes <- summ.var(data = data) ## Remove genes with near to zero variance
data <- data[!nonVar.genes$nzv, ]


## Annotate Probes to GeneSymbol
## Download GPL10558 db from NCBI GEO
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL21970
genome <- getGEO("GPL10558")@dataTable@table[, c("Symbol", "ID")]
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
data <- temp


metadata <- phenoData(gset) ## Download clinical data from NCBI GEO
metadata <- pData(metadata)

## Separe Healthy controls and SLE samples
HC.ln3 <- data[, rownames(metadata)[metadata$`group:ch1` == "Healthy control of SLE"]]

metadata <- metadata[ifelse(metadata$"group:ch1" == "SLE", T, F), ]

## Select relevant variables for the analysis
metadata <- metadata[, c(
  "title", "collection date:ch1", "group:ch1",
  "renal:ch1", "class kb:ch1", "ever renal:ch1", "condition:ch1", "sledai:ch1"
)]
colnames(metadata) <- c("title", "date", "diagnosis", "LN", "biopsy", "everRenal", "condition", "sledai")

patientInfo <- do.call("rbind", str_split(metadata$title, "V"))

metadata <- data.frame(
  "PatientID" = patientInfo[, 1],
  "Visit" = as.numeric(patientInfo[, 2]), metadata
)

metadata <- metadata[order(metadata$PatientID, metadata$Visit, decreasing = F), ]

## Date to numeric measurement
metadata$date <- stri_replace_all_regex(
  metadata$date,
  pattern = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"),
  replacement = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"),
  vectorize = FALSE
)

metadata$TimeFromVisits <- NA
patients <- unique(metadata$PatientID)
for (i in 1:length(patients)) {
  tmp <- metadata[metadata$PatientID == patients[i], ]
  tmp <- tmp[order(tmp$Visit, decreasing = F), ]
  ref <- tmp[1, "date"]
  
  for (j in 1:nrow(tmp)) {
    x <- as.numeric(difftime(strptime(tmp[j, "date"], format = "%d-%m-%Y"),
                             strptime(ref, format = "%d-%m-%Y"),
                             units = "days"
    ))
    metadata[metadata$PatientID == patients[i] &
               metadata$Visit == tmp$Visit[j], "TimeFromVisits"] <- x
  }
}

metadata$biopsyClass <- ifelse(is.na(metadata$biopsy), NA,
                               ifelse(str_detect(metadata$biopsy, "III") | str_detect(metadata$biopsy, "IV"), "Prolif", "Other")
)

SLE.ln3 <- data[, rownames(metadata)]
table(metadata$biopsy)

## Select only one sample by patient (only have Biopsy information for the first sample of each patient)
pats <- unique(metadata$PatientID)
sel <- NULL
for (i in 1:length(pats)) {
  tmp <- metadata[metadata$PatientID == pats[i], ]
  #if ("Y" %in% tmp$LN) {
  tmp <- tmp[order(tmp$biopsy), ]
  sel <- c(sel, rownames(tmp)[1])
  # } else {
  #   sel <- c(sel, c(rownames(tmp)[sample(1:nrow(tmp), 1)]))
  # }
}

metadata <- metadata[sel, ]

metadata<-metadata[is.na(metadata$biopsyClass) | metadata$biopsyClass=="Prolif",]

metadata$pLN <- ifelse(is.na(metadata$biopsyClass), "NO",
                       ifelse(metadata$biopsyClass == "Prolif", "YES", "NO")
)

SLE.ln3 <- data[, rownames(metadata)]

## Remove genes with near-zero variance
nonVar.genes<-summ.var(data = HC.ln3); table(!nonVar.genes$nzv)
HC.ln3<-HC.ln3[!nonVar.genes$nzv,]
SLE.ln3<-SLE.ln3[!nonVar.genes$nzv,]

apply(SLE.ln3,1,sd)[order(apply(SLE.ln3,1,sd),decreasing=FALSE)][1:10]
apply(HC.ln3,1,sd)[order(apply(HC.ln3,1,sd),decreasing=FALSE)][1:10]


x<-prcomp(t(SLE.ln3))
D<-as.data.frame(x$x[,1:2])
D<-cbind(metadata,D)

ggplot(D,aes(x=PC1,y=PC2,fill=pLN)) + theme_classic() + geom_point(shape=21,size=3)+
  scale_fill_manual(values=c("#cd5554","#0198b1","#8fbc8f","gold"))

clin.ln3<-metadata
table(clin.ln3$pLN)

rm(D,data,genome,gset,metadata,nonVar.genes,patientInfo,temp,tmp,x,i,j,
   patients,pats,ref,sel,finalGenes)




##······································································· Step 4
## Petri Cohort ----
## All clinical data can be requested to the original authors (GSE45291)


## Onset renal manifestation and 
ren<-read.csv("D:/Work/BMS/PMPS/PMSPbenchmark/Datasets/History1_Biogen_20170731.csv",sep=",")
ren<-ren[,c("Cohort.ID","Renal.SLE","PROTEINU","PROTN_MO","PROTN_YR","HEMATURI","HEMAT_MO","HEMAT_YR",
            "NEPHROTC","NEPHR_MO","NEPHR_YR","INSUFF","INSUFF_MO","INSUFF_YR","FAILURE","FAIL_MO","FAIL_YR",
            "REBIOPSY","REBX_MO1","REBX_YR1","HISTOL1",
            "REBX_MO2","REBX_YR2","HISTOL2",
            "REBX_MO3","REBX_YR3","HISTOL3")]

## No LN
noren<-ren[ren$Renal.SLE==0 & ren$INSUFF==0 & ren$FAILURE==0,
           c("Cohort.ID","Renal.SLE")]

## LN
ren<- ren %>% filter((HISTOL1 !="" & !is.na(HISTOL1)) |
                       (HISTOL2 !="" & !is.na(HISTOL2)) |
                       (HISTOL3 !="" & !is.na(HISTOL3)))

## Formating dates
ren$BX_YR1<-NA
ren$BX_YR2<-NA
ren$BX_YR3<-NA
for(i in 1:nrow(ren)){
  if(!is.na(ren$REBX_YR1[i])){
    if(is.na(ren$REBX_MO1[i])){ren$REBX_MO1[i]<-"6"}
    ren$BX_YR1[i]<-paste0(ren$REBX_YR1[i],"-",ren$REBX_MO1[i],"-01")
  }
  if(!is.na(ren$REBX_YR2[i])){
    if(is.na(ren$REBX_MO2[i])){ren$REBX_MO2[i]<-"6"}
    ren$BX_YR2[i]<-paste0(ren$REBX_YR2[i],"-",ren$REBX_MO2[i],"-01")
  }
  if(!is.na(ren$REBX_YR3[i])){
    if(is.na(ren$REBX_MO3[i])){ren$REBX_MO3[i]<-"6"}
    ren$BX_YR3[i]<-paste0(ren$REBX_YR3[i],"-",ren$REBX_MO3[i],"-01")
  }
}
ren<-ren %>% dplyr::select(Cohort.ID,BX_YR1,HISTOL1,BX_YR2,HISTOL2,BX_YR3,HISTOL3)

## Clinical - by Sample

clin<-read.csv("D:/Work/BMS/PMPS/PMSPbenchmark/Datasets/metadataPetri.csv",sep="\t")
clin<-clin[,c("Visit.IDN","CohortID","Date")]
colnames(clin)<-c("visitID","Cohort.ID","date")
clin$date<-dates_converted <- format(as.Date(clin$date, format = "%m/%d/%Y"), format = "%Y-%m-%d")


clinnoren<-clin[clin$Cohort.ID %in% noren$Cohort.ID,]
#clinnoren<-clinnoren[!duplicated(clinnoren$Cohort.ID),] ## No LN (only one sample by patient)

## LN
clin<-clin[clin$Cohort.ID %in% ren$Cohort.ID,]

clin<- left_join(clin, ren, by = "Cohort.ID")
clin[clin==""]<-NA

clin <- clin %>%
  mutate(across(contains("HISTOL"), ~ str_replace_all(., c("Focal Glomerulosclerosis" = NA,
                                                           "lupus nephritis" = NA,
                                                           "Calcium phosphate, glomerulosclerosis" = NA,
                                                           "Plasma cell" = NA,
                                                           "fibrillary GN"=NA,
                                                           "Interstitial"=NA,
                                                           "Thrombotic"=NA,
                                                           "Transplant rejection"=NA,
                                                           "II sclerosis"="II",
                                                           "Class V"="V",
                                                           "III,V"="III, V",
                                                           "Class 2, 6"= "II, VI",
                                                           "3"="III",
                                                           "No Lupus"=NA,
                                                           "Class 2" = "II",
                                                           "CLASS VI"="VI",
                                                           "II, FSGS" = "II",
                                                           "class V"="V",
                                                           "VI glomerular hypertrophy"="VI"
  ))))

rm(dates_converted,i)


#' Medir el tiempo para las 3 biopsias
#' si esta dentro de 1year, seleccionar

clin$Selection<-NA
for(i in 1:nrow(clin)){
  res<-c("HISTOL1"=ifelse(is.na(clin[i,"HISTOL1"]),10,as.numeric(difftime(clin[i,"date"], clin[i,"BX_YR1"], units = "days"))/365),
         "HISTOL2"=ifelse(is.na(clin[i,"HISTOL2"]),10,as.numeric(difftime(clin[i,"date"], clin[i,"BX_YR2"], units = "days"))/365),
         "HISTOL3"=ifelse(is.na(clin[i,"HISTOL3"]),10,as.numeric(difftime(clin[i,"date"], clin[i,"BX_YR3"], units = "days"))/365))
  res<-res[order(abs(res),decreasing = F)]
  #print(res)
  
  if(abs(res[1])<=1){
    print(i)
    clin$Selection[i]<-names(res)[1]
  }
}
clin<-clin[!is.na(clin$Selection),c("Cohort.ID","date","BX_YR1","BX_YR2","BX_YR3","HISTOL1",
                                    "HISTOL2","HISTOL3","Selection")]


## Get sample ID
matchID<-read.csv("D:/Work/BMS/PMPS/PMSPbenchmark/Datasets/cellFile2sampleID.csv",
                  sep=";")
matchID<-matchID[!is.na(matchID$Visit_Date),c("Subject_ID","Visit_Date","GZ_Filenames")]
colnames(matchID)<-c("Cohort.ID","date","sampleID")
matchID$date<-dates_converted <- format(as.Date(matchID$date, format = "%d/%m/%Y"), format = "%Y-%m-%d")


## No LN
clinnoren<-left_join(clinnoren, matchID, by = c("Cohort.ID","date"))
#clinnoren<-clinnoren[sample(1:nrow(clinnoren),nrow(clinnoren),replace = F),]
clinnoren<-clinnoren[!duplicated(clinnoren$Cohort.ID),] ## No LN (only one sample by patient)


## Select only patients with pLN
sel<-NULL
for(i in 1:nrow(clin)){
  if(grepl("III",clin[i,clin$Selection[i]]) | grepl("IV",clin[i,clin$Selection[i]])){
    sel<-c(sel,i)
  }
}
clin<-clin[sel,]


## Select one visit for each patient (close to BX): For test
clin$sel2<-ifelse(clin$Selection=="HISTOL1","BX_YR1",
                  ifelse(clin$Selection=="HISTOL2","BX_YR2","BX_YR3"))
CLIN<-NULL
pats<-unique(clin$Cohort.ID)
for(i in 1:length(pats)){
  print(i)
  tmp<-clin[clin$Cohort.ID==pats[i],,drop=FALSE]
  
  TimeFromBX<-NULL
  nephrClass<-NULL
  for(j in 1:nrow(tmp)){
    x<-as.numeric(difftime(tmp[j,"date"], tmp[j,tmp$sel2[j]], units = "days"))/365
    TimeFromBX<-c(TimeFromBX,x)
    nephrClass<-c(nephrClass,tmp[j,tmp$Selection[j]])
  }
  tmp$nephrClass<-nephrClass
  tmp$TimeFromBX<-TimeFromBX
  tmp<-tmp[order(tmp$TimeFromBX,decreasing = FALSE),]
  CLIN<-rbind(CLIN,tmp[1,,drop=FALSE])
  
}
clin<-CLIN[,c("Cohort.ID","date","nephrClass","TimeFromBX")]


## LN
clin<- left_join(clin, matchID, by = c("Cohort.ID","date"))

clinnoren$TimeFromBX<-NA
clinnoren$nephrClass<-NA
clinnoren<-clinnoren[,c("Cohort.ID","date","nephrClass","TimeFromBX","sampleID")]

clin$pLN<-"YES"
clinnoren$pLN<-"NO"

clin.ln4<-rbind(clin,clinnoren)


rm(clin,CLIN,clinnoren,noren,ren,tmp,dates_converted,i,j,nephrClass,pats,res,
   sel,TimeFromBX,x)


## Load GeneExpression
## Gene expression can be also downloaded from GEO (GSE45291)
data<-read.csv(file="D:/Work/Update_MyPROSLE/Datasets/PetriALL.txt",
               sep="\t",row.names = "GeneSymbol")

## Preprocessing (Log normalization and filtering non variable genes)
data<-norm.log(data)
nonVar.genes<-summ.var(data = data)
data<-data[!nonVar.genes$nzv,] ## remove genes with near-zero variance
colnames(data)<-gsub(pattern = "X",replacement = "",x = colnames(data))

## Separate Healthy and SLE samples
clin<-read.csv(file="D:/Work/Update_MyPROSLE/Datasets/Metadata.petri.csv",
               sep=";",row.names = "GZ_Filenames")

HC.ln4<-data[,rownames(clin)[clin$Diagnosis=="Healthy"]]
SLE.ln4<-data[,clin.ln4$sampleID]
rownames(clin.ln4)<-clin.ln4$sampleID


## Remove genes with near-zero variance
nonVar.genes<-summ.var(data = HC.ln4); table(!nonVar.genes$nzv)
HC.ln4<-HC.ln4[!nonVar.genes$nzv,]
SLE.ln4<-SLE.ln4[!nonVar.genes$nzv,]

apply(SLE.ln4,1,sd)[order(apply(SLE.ln4,1,sd),decreasing=FALSE)][1:10]
apply(HC.ln4,1,sd)[order(apply(HC.ln4,1,sd),decreasing=FALSE)][1:10]



x<-prcomp(t(SLE.ln4))
D<-as.data.frame(x$x[,1:2])
D<-cbind(clin.ln4,D)

ggplot(D,aes(x=PC1,y=PC2,fill=pLN)) + theme_classic() + geom_point(shape=21,size=3)+
  scale_fill_manual(values=c("#cd5554","#0198b1","#8fbc8f","gold"))


rm(list=setdiff(ls(),c("clin.ln1","clin.ln2","clin.ln3","clin.ln4",
                       "SLE.ln1","SLE.ln2","SLE.ln3","SLE.ln4",
                       "HC.ln1","HC.ln2","HC.ln3","HC.ln4")))

save.image("D:/Work/BMS/PMPS/PMSPbenchmark/RData/LNDatasets.RData")


## Plots
##----
x<-prcomp(t(SLE.ln1))
D<-as.data.frame(x$x[,1:2])
D<-cbind(clin.ln1,D)
D$pLN<-ifelse(D$NeprBX=="NoLN","NoLN","pLN")
D$pLN<-factor(D$pLN,levels=c("NoLN","pLN"))

p1<-ggplot(D,aes(x=PC1,y=PC2,fill=pLN)) + theme_classic() + geom_point(shape=21,size=1.25,stroke=0.35)+
  stat_ellipse(geom="polygon", aes(color = pLN),
               alpha = 0.0, linetype=10, show.legend = FALSE, level = 0.9) +
  scale_fill_manual(values=c("pLN"="#cd5554","NoLN"="#B0B0B0"))+
  scale_color_manual(values=c("pLN"="#cd5554","NoLN"="#B0B0B0")) +
  theme(axis.line = element_line(linewidth=0.5),
        axis.text = element_text(size=7),
        axis.ticks = element_line(linewidth=0.5),
        axis.title = element_text(size=8))

##----
x<-prcomp(t(SLE.ln2))
D<-as.data.frame(x$x[,1:2])
D<-cbind(clin.ln2,D)
D$pLN<-ifelse(D$LN=="NO","NoLN","pLN")
D$pLN<-factor(D$pLN,levels=c("NoLN","pLN"))

p2<-ggplot(D,aes(x=PC1,y=PC2,fill=pLN)) + theme_classic() + geom_point(shape=21,size=1.25,stroke=0.35)+
  stat_ellipse(geom="polygon", aes(color = pLN),
               alpha = 0.0, linetype=10, show.legend = FALSE, level = 0.9) +
  scale_fill_manual(values=c("pLN"="#cd5554","NoLN"="#B0B0B0"))+
  scale_color_manual(values=c("pLN"="#cd5554","NoLN"="#B0B0B0")) +
  theme(axis.line = element_line(linewidth=0.5),
        axis.text = element_text(size=7),
        axis.ticks = element_line(linewidth=0.5),
        axis.title = element_text(size=8))


##----
x<-prcomp(t(SLE.ln3))
D<-as.data.frame(x$x[,1:2])
D<-cbind(clin.ln3,D)
D$pLN<-ifelse(D$pLN=="NO","NoLN","pLN")
D$pLN<-factor(D$pLN,levels=c("NoLN","pLN"))

p3<-ggplot(D,aes(x=PC1,y=PC2,fill=pLN)) + theme_classic() + geom_point(shape=21,size=1.25,stroke=0.35)+
  stat_ellipse(geom="polygon", aes(color = pLN),
               alpha = 0.0, linetype=10, show.legend = FALSE, level = 0.9) +
  scale_fill_manual(values=c("pLN"="#cd5554","NoLN"="#B0B0B0"))+
  scale_color_manual(values=c("pLN"="#cd5554","NoLN"="#B0B0B0")) +
  theme(axis.line = element_line(linewidth=0.5),
        axis.text = element_text(size=7),
        axis.ticks = element_line(linewidth=0.5),
        axis.title = element_text(size=8))


##----
x<-prcomp(t(SLE.ln4))
D<-as.data.frame(x$x[,1:2])
D<-cbind(clin.ln4,D)
D$pLN<-ifelse(D$pLN=="NO","NoLN","pLN")
D$pLN<-factor(D$pLN,levels=c("NoLN","pLN"))

p4<-ggplot(D,aes(x=PC1,y=PC2,fill=pLN)) + theme_classic() + geom_point(shape=21,size=1.25,stroke=0.35)+
  stat_ellipse(geom="polygon", aes(color = pLN),
               alpha = 0.0, linetype=10, show.legend = FALSE, level = 0.9) +
  scale_fill_manual(values=c("pLN"="#cd5554","NoLN"="#B0B0B0"))+
  scale_color_manual(values=c("pLN"="#cd5554","NoLN"="#B0B0B0")) +
  theme(axis.line = element_line(linewidth=0.5),
        axis.text = element_text(size=7),
        axis.ticks = element_line(linewidth=0.5),
        axis.title = element_text(size=8))


library(ggpubr)

ggarrange(p1,p2,p3,p4,ncol=2,nrow=2,common.legend = TRUE,legend = "bottom")

rm(D,x,p1,p2,p3,p4)

