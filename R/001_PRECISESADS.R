#'··············································································
#' Personalized molecular profile scoring benchmark
#' R version 4.3.1 (2024-04-10)
#' 
#'··············································································
#' Process data from PRECISESADS cohort:

#' SLE and HC samples from 376 and 263 individuals. Demographical, serological,
#' cellular and other clinical data must be also correcly formatted.

#' Clinical and gene expression data from PRECISESADS can be requested to the 
#' original authors: https://acrjournals.onlinelibrary.wiley.com/doi/10.1002/art.41610

##······································································· Step 0
## Set environment ----

setwd("C:/Users/danie/Desktop/WORK/BENCHMARK_25")
set.seed(1234)
source(paste0(getwd(),"/code/utils.R"))

pkgs<-c("parallel","matrixStats","biomaRt","NOISeq","stringr","dplyr","tidyr",
        "doParallel","caret","pbapply","BiocParallel","tibble")
check.packages(pkgs); rm(pkgs)

##······································································· Step 1
## Get Gene Expression data (PRECISESADS) ----
#' 2 sets of patients and 1 of healthy controls, without batch effect)

## First SLE set (and controls) # "/mnt/data/PMPS_Datasets/Count.PRECISESADS.csv"
Count.precisesads<-read.table(paste0(getwd(),"/Datasets/Count.PRECISESADS.csv"),
                              header=T,sep=";",row.names=1,dec=",",stringsAsFactors = F)
Count.precisesads<-type.convert(x = Count.precisesads,as.is=F)

## Filter non-expressed genes
Count.precisesads<-Count.precisesads[rownames(Count.precisesads)[rowCounts(Count.precisesads>=10)>=10],]

## Annotate to Gene symbol
#load(paste0(getwd(),"/RData/ensembl.RData")) # if biomaRt is down
Count.precisesads<-annotateGenes(data=Count.precisesads,
                                 toGenes='external_gene_name',
                                 fromGenes='ensembl_gene_id',
                                 method = "median") #,ensembl = ensembl)

## Tmm and Log2 normalization and filtering non variable gene
data<-tmm(Count.precisesads)
data<-log2(data+1)
nonVar.genes<-summ.var(data = data)
data<-data[!nonVar.genes$nzv,] 
colnames(data)<-gsub(pattern = "X",replacement = "",x = colnames(data))

precisesads1<-data
rm(Count.precisesads,data,nonVar.genes)

## Second SLE set # "/mnt/data/PMPS_Datasets/Count.PRECISESADS2.csv"
Count.precisesads2<-read.table(paste0(getwd(),"/Datasets/Count.PRECISESADS2.csv"),
                               header=T,sep=";",row.names=1,dec=",",stringsAsFactors = F)
Count.precisesads2<-type.convert(x = Count.precisesads2,as.is=F)

## Filter non-expressed genes
Count.precisesads2<-Count.precisesads2[rownames(Count.precisesads2)[rowCounts(Count.precisesads2>=10)>=10],]

## Annotate to Gene symbol
Count.precisesads2<-annotateGenes(data=Count.precisesads2,
                                  toGenes='external_gene_name',
                                  fromGenes='ensembl_gene_id',
                                  method = "median") #,ensembl = ensembl)

## Tmm and Log2 normalization and filtering non variable gene
data<-tmm(Count.precisesads2)
data<-log2(data+1)
nonVar.genes<-summ.var(data = data)
data<-data[!nonVar.genes$nzv,] 
colnames(data)<-gsub(pattern = "X",replacement = "",x = colnames(data))

precisesads2<-data
rm(Count.precisesads2,data,nonVar.genes)

## Select common genes
genes<-as.character(intersect(rownames(precisesads1),rownames(precisesads2)))
precisesads1<-as.data.frame(precisesads1[genes,]) 
precisesads2<-as.data.frame(precisesads2[genes,])


## Select SLE and Healthy samples # "/mnt/data/PMPS_Datasets/Metadata.PRECISESADS.csv"
clin1<-read.csv(file=paste0(getwd(),"/Datasets/Metadata.PRECISESADS.csv"),
                header=T,sep=";",row.names=1,dec=",")
clin1<-clin1[colnames(precisesads1),c(2,8)]

SLE<-precisesads1[,rownames(clin1[ifelse(clin1$Diagnosis=="SLE",T,F),])]
HC.prec<-precisesads1[,rownames(clin1[ifelse(clin1$Diagnosis=="CTRL",T,F),])]  


# "/mnt/data/PMPS_Datasets/Metadata.PRECISESADS2.csv"
clin2<-read.csv(file= paste0(getwd(),"/Datasets/Metadata.PRECISESADS2.csv"),
                header=T,sep=";",row.names=1,dec=",")
clin2<-clin2[colnames(precisesads2),c(1,5)]

SLE2<-precisesads2[,rownames(clin2[ifelse(clin2$Diagnosis=="SLE",T,F),])]
SLE.prec<-cbind(SLE,SLE2)

## Remove flat genes in HC
nonVar.genes<-summ.var(data = HC.prec); table(nonVar.genes$nzv)
SLE.prec<-SLE.prec[!nonVar.genes$nzv,]
HC.prec<-HC.prec[!nonVar.genes$nzv,]


rm(clin1,clin2,SLE,genes,SLE2,precisesads1,precisesads2,nonVar.genes)

##······································································· Step 2
## Clinical manifestations ----

clin.prec<-list()

# "/mnt/data/PMPS_Datasets/PRECISESADSclin.tsv"
symt<-read.csv(paste0(getwd(),"/Datasets/PRECISESADSclin.tsv"),sep="\t",
               row.names = "Sampling_OMIC.number")
rownames(symt)<-gsub(pattern = "N",replacement = "",x = rownames(symt)) 

symt<-symt[intersect(rownames(symt),colnames(SLE.prec)),
           c("Diagnosis_Arm",
             "Symptom_Disease.activity",
             "Diagnosis_Age.at.onset",
             "Sampling_Disease.duration",
             "Comorbidity_Abdominal.pain",
             "Comorbidity_Diarrhea..recurrent.",
             "Comorbidity_Dyslipidemia",
             "Comorbidity_Hypertension",
             "Comorbidity_Obesity..BMI....30.",
             "Comorbidity_Stipsis.Constipation",
             "Comorbidity_Thyroiditis",
             "Gastro_History.of.esophageal.reflux.disease",
             "Heart_Systemic.hypertension",
             "Heart_Valve.lesions",
             "Kidney_Abnormal.Creatinine",
             "Kidney_Abnormal.lipid.profile",
             "Kidney_Abnormal.urine.analysis",
             "Kidney_Biopsy.proven.nephritis",
             "Kidney_Proteinuria",
             "Kidney_Urine.proteins",
             "Lab_Reduced.C3.levels",
             "Lab_Reduced.C4.levels",
             "Lab_Reduced.White.Blood.Cell.Count",
             "Muscle.and.Skeletal_Arthritis",
             "Skin_Cutaneous.Lupus.active",
             "Skin_Cutaneous.Lupus.chronic",
             "Skin_Photosensitivity",
             "Skin_Sicca.syndrome",
             "Skin_Telangectasia",
             "Symptom_Abnormal.inflammatory.indexes",
             "Symptom_Fever",
             "Symptom_Hypergammabulinemia",
             "Vascular_History.of.Raynauds.phenomenon",
             "Vascular_History.of.recurrent.miscarriage.or.pregnancy.complications")]

symt[symt =="Unknown"]<-NA
symt[symt =="Past"]<-NA
symt[symt =="N/A"]<-NA
symt[symt =="Present"]<-"Yes"
symt[symt =="Moderate"]<-"Yes"
symt[symt =="Severe"]<-"Yes"
symt$Symptom_Disease.activity<-as.numeric(symt$Symptom_Disease.activity)
symt$Diagnosis_Age.at.onset<-as.numeric(symt$Diagnosis_Age.at.onset)
symt$Sampling_Disease.duration<-as.numeric(symt$Sampling_Disease.duration)

clin.prec[[1]]<-symt
names(clin.prec)[1]<-"symptoms"
rm(symt)


##······································································· Step 3
## Serology
# "/mnt/data/PMPS_Datasets/prec_cytokines.txt"
cyt<-read.csv(paste0(getwd(),"/Datasets/prec_cytokines.txt"),header=T,
              sep="\t",dec=".",row.names = "OMICID")
cyt<-cyt[intersect(colnames(SLE.prec),rownames(cyt)),c(106:ncol(cyt))]

## Select the most measured cytokines
cyt<-cyt[,colnames(cyt) %in% c("BAFF_ELISA","BLC","CRP_ELISA","FAS_LIGAND",
                               "GDF_15", "IL_1_RA","IL_1_RII","IL_6_ELISA",
                               "IP_10","MCP_2","MCP_4","MIP_1_BETA",
                               "MMP_2_ELISA","MMP_8","TARC","TGF_BETA_ELISA",
                               "TNF_ALPHA_ELISA","TNF_RI")]
colnames(cyt)<-gsub(pattern = "_ELISA",replacement = "",x = colnames(cyt))

clin.prec[[2]]<-cyt
names(clin.prec)[2]<-"cytokines"
rm(cyt)


##······································································· Step 4
## Autoantybodies
# "/mnt/data/PMPS_Datasets/prec_autoantibodies.txt"
aant<-read.csv(paste0(getwd(),"/Datasets/prec_autoantibodies.txt"),header=T,
               sep="\t",dec=",",row.names = "OMICID")
aant<-aant[intersect(colnames(SLE.prec),rownames(aant)),
           str_detect(string = colnames(aant),pattern = "CALL",negate = F)]
colnames(aant)<-gsub(pattern = "_CALL",replacement = "",x = colnames(aant))

## Filter antibodies
aant<-aant[,unlist(lapply(apply(aant,2,table),function(x){
  ifelse(dim(x)<2 | min(x)<10,F,T)}))]

clin.prec[[3]]<-aant
names(clin.prec)[3]<-"autoantibodies"
rm(aant)


##······································································· Step 5
## Cell percentages
# "/mnt/data/PMPS_Datasets/prec_flowcytometry.csv"
cells<-read.csv(paste0(getwd(),"/Datasets/prec_flowcytometry.csv"),header=T,
                sep=";",dec=".",row.names = "OMICID")
cells<-cells[intersect(colnames(SLE.prec),rownames(cells)),]

cells<- cells %>% mutate(rowIDs = rownames(cells)) %>% drop_na(PBMC,PMN) %>% 
  mutate(total = PBMC + PMN) %>% column_to_rownames("rowIDs")
cells <- cells/cells$total

cells<- cells %>% select(-all_of(c("LEUKOCYTES","LYMPHOCYTES","MONOCYTES",
                                   "PBMC","PMN","total")))

clin.prec[[4]]<-cells
names(clin.prec)[4]<-"cells"

## rm all
rm(list=setdiff(ls(),c("clin.prec","HC.prec","SLE.prec")))

# "/mnt/data/PMPS_Datasets/PRECISESADS.RData"
save.image(paste0(getwd(),"/RData/PRECISESADS.RData"))


