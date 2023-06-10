setwd("D:/高端paper work/AI/血管/AI结合循环血小板mRNA预测隐匿性肿瘤/前期/正式图表/投稿使用/Fundamental Research/InTVeld_Pancancer_TSOO-main")
library(doMC)
library(foreach)
library(RUVSeq)
library(ppso)
library(GEOquery)
library(readxl)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(tibble)
library(tableone)
library(pROC)
library(plotrix)
library(ggplot2)
library(circlize)
library(tidyverse)
library(limma)
library(dplyr)
library(GseaVis)
library(BiocManager)
library(dplyr)
library(GOplot)
library(readxl)
library(ggplot2)
library(cowplot)
library(enrichplot)
library(GOplot)
source('bin/thromboSeqTools_PreProcessing_2.R')
source('bin/thromboSeqTools_ANOVA.R')
source('bin/thromboSeqTools_PSO.R')
nCores = 50
###################################################################################################################################Prepare Training, Evaluation and Validation 
#切割训练集和验证集
data <- read.csv("TableS2.csv")
data$group <- data$Group
data$Group[which(data$Group %in% c("Breast cancer","Cholangiocarcinoma","Colorectal cancer",       "Endometrial cancer",
                                   "Esophageal carcinoma","Glioma","Head and neck cancer","Hepatocellular carcinoma",
                                   "Lymphoma","Melanoma","Multiple Myeloma","Non-small-cell lung cancer","Ovarian cancer",
                                   "Pancreatic cancer","Prostate cancer","Renal cell cancer","Sarcoma","Urothelial cancer"))] <- "Cancer"
# here all non-cancer groups are names 'asymptomaticControls' to enable validation purposes
data$Group[which(data$Group %in% c("Asymptomatic controls","Angina pectoris","Bowel disease","Former sarcoma","Hematuria",
                                   "Medically-intractable epilepsy","Multiple sclerosis","nSTEMI",
                                   "Pancreatic diseases","Pulmonary Hypertension"))] <- "asymptomaticControls"
data <- subset(data,data$Group=="Cancer")
nrow(data)
set.seed(121111)
sample <- data$Sample.ID
validation <- sample(nrow(data), round(nrow(data)*0.4))
data.validation <- data[validation, ]
data.validation$Validation.series = 1
data.validation$Training.series = 0
data.validation$Evaluation.series = 0
data1 <- data[-validation, ]
Training <- sample(nrow(data1), round(nrow(data1)*0.667))
data.Training <- data1[Training, ]
data.Training$Validation.series = 0
data.Training$Training.series = 1
data.Training$Evaluation.series = 0
data.Evaluation <- data1[-Training, ]
data.Evaluation$Validation.series = 0
data.Evaluation$Training.series = 0
data.Evaluation$Evaluation.series = 1
data <- rbind(data.validation,data.Training,data.Evaluation)
setwd("D:/高端paper work/AI/血管/AI结合循环血小板mRNA预测隐匿性肿瘤/前期/正式图表/投稿使用/返修/cluster3 4")
data_derivation_clusters <- read.csv("output_train_3cluster_evaluation.csv")
data.Training$X %in% data_derivation_clusters$Unnamed..0
data_derivation_clusters <- read.csv("output_test_3cluster_evaluation.csv")
data.Evaluation$X %in% data_derivation_clusters$Unnamed..0
data_derivation_clusters <- read.csv("output_validation_3cluster_evaluation.csv")
data.validation$X %in% data_derivation_clusters$Unnamed..0

nrow(data)

setwd("D:/高端paper work/AI/血管/AI结合循环血小板mRNA预测隐匿性肿瘤/前期/正式图表/投稿使用/返修/cluster3 4")
data$cohort <- ifelse(data$Training.series==1,"Training",ifelse(data$Evaluation.series==1,"Evaluation","Validation"))
data$Stage <- ifelse(data$Stage=="n.a.",NA,data$Stage)
data$Sex <- ifelse(data$Sex=="n.a.",NA,data$Sex)
data$Age <- ifelse(data$Age=="n.a.",NA,data$Age)
data$Age  <- as.numeric(data$Age)

data$`BRCA` <- ifelse(data$group=="Breast cancer",1,0)
data$`CHOL` <- ifelse(data$group=="Cholangiocarcinoma",1,0)
data$`CRC` <- ifelse(data$group=="Colorectal cancer",1,0)
data$`ENDO` <- ifelse(data$group=="Endometrial cancer",1,0)
data$`ESO` <- ifelse(data$group=="Esophageal carcinoma",1,0)
data$`GLIO` <- ifelse(data$group=="Glioma",1,0)
data$`HNSSC` <- ifelse(data$group=="Head and neck cancer",1,0)
data$`HCC` <- ifelse(data$group=="Hepatocellular carcinoma",1,0)
data$`LYM` <- ifelse(data$group=="Lymphoma",1,0)
data$`MELA` <- ifelse(data$group=="Melanoma",1,0)
data$`MM` <- ifelse(data$group=="Multiple Myeloma",1,0)
data$`NSCLC` <- ifelse(data$group=="Non-small-cell lung cancer",1,0)
data$`OVCA` <- ifelse(data$group=="Ovarian cancer",1,0)
data$`PDAC` <- ifelse(data$group=="Pancreatic cancer",1,0)
data$`PRCA` <- ifelse(data$group=="Prostate cancer",1,0)
data$`RCC` <- ifelse(data$group=="Renal cell cancer",1,0)
data$`SARC` <- ifelse(data$group=="Sarcoma",1,0)
data$`URO` <- ifelse(data$group=="Urothelial cancer",1,0)
data$`I Stage` <- ifelse(is.na(data$Stage)==T,NA,ifelse(data$Stage=="I",1,0))
data$`II Stage` <- ifelse(is.na(data$Stage)==T,NA,ifelse(data$Stage=="II",1,0))
data$`III Stage` <- ifelse(is.na(data$Stage)==T,NA,ifelse(data$Stage=="III",1,0))
data$`IV Stage` <-ifelse(is.na(data$Stage)==T,NA,ifelse(data$Stage=="IV",1,0))

###sTable 1
varsToFactor <- c("Sex","BRCA","CHOL","CRC","ENDO","ESO","GLIO","HNSSC","HCC","LYM","MELA","MM","NSCLC","OVCA","PDAC","PRCA","RCC","SARC","URO","I Stage","II Stage","III Stage","IV Stage")
vars <- c("Sex","Age","BRCA","CHOL","CRC","ENDO","ESO","GLIO","HNSSC","HCC","LYM","MELA","MM","NSCLC","OVCA","PDAC","PRCA","RCC","SARC","URO","I Stage","II Stage","III Stage","IV Stage")
data[varsToFactor] <- lapply(data[varsToFactor],factor)
tableOne <- CreateTableOne(vars = vars,strata = 'cohort', data = data)
tableOne <-print(tableOne,nonnormal = c("Age"))
write.csv(tableOne,"table 1.csv")
###################Prepare_data####################################################################################################Prepare_data
setwd("D:/高端paper work/AI/血管/AI结合循环血小板mRNA预测隐匿性肿瘤/前期/正式图表/投稿使用/Fundamental Research/InTVeld_Pancancer_TSOO-main")
load("TEP_Count_Matrix.RData")
colnames(TEP_Count_Matrix) <- sub(".*?-", "", colnames(TEP_Count_Matrix))
# load and prepare sample info file provided in GEO. Identical to Table S2 of the manuscript.
sampleInfo <- read.csv('TableS2.csv', sep=",", row.names=1)
sampleInfo$rownames <- rownames(sampleInfo)
for (i in 1:3){
  sampleInfo <- sampleInfo[,-7]
}
data <- data[,c("Sample.ID","Training.series","Evaluation.series","Validation.series")]
sampleInfo <- merge(sampleInfo,data,by="Sample.ID",all.x=T)
sampleInfo <- sampleInfo[order(match(sampleInfo$rownames, colnames(TEP_Count_Matrix))), ]
# load gene info and create edgeR object
load('bin/dgeGenesEnsembl75.RData')
library(edgeR)
dge <- DGEList(counts = TEP_Count_Matrix,
               group = sampleInfo$Group,
               genes = genes[which(rownames(genes) %in% rownames(TEP_Count_Matrix)),]
)
dge$samples <- cbind(dge$samples, sampleInfo)
dge$samples$lib.size <- sampleInfo$lib.size

# dichotomize the groups (i.e. Cancer and asymptomatic controls)
dge$samples$group <- factor(dge$samples$group, levels = c(levels(dge$samples$group),"Cancer","asymptomaticControls"))
dge$samples$group[which(dge$samples$group %in% c("Breast cancer","Cholangiocarcinoma","Colorectal cancer",       "Endometrial cancer",
                                                 "Esophageal carcinoma","Glioma","Head and neck cancer","Hepatocellular carcinoma",
                                                 "Lymphoma","Melanoma","Multiple Myeloma","Non-small-cell lung cancer","Ovarian cancer",
                                                 "Pancreatic cancer","Prostate cancer","Renal cell cancer","Sarcoma","Urothelial cancer"))] <- "Cancer"
# here all non-cancer groups are names 'asymptomaticControls' to enable validation purposes
dge$samples$group[which(dge$samples$group %in% c("Asymptomatic controls","Angina pectoris","Bowel disease","Former sarcoma","Hematuria",
                                                 "Medically-intractable epilepsy","Multiple sclerosis","nSTEMI",
                                                 "Pancreatic diseases","Pulmonary Hypertension"))] <- "asymptomaticControls"
dge = dge[,c(colnames(dge)[dge$samples$group=="Cancer"])]
dge$samples$group <- dge$samples$Group
dge$samples <- droplevels(dge$samples)
summary(as.factor(dge$samples$group))
summary(as.factor(dge$samples$Training.series))
summary(as.factor(dge$samples$Evaluation.series))
summary(as.factor(dge$samples$Validation.series))

###############################################excluded the genes that yielded <30 intron-spanning reads in >90% of the dataset for the samples
dge$samples$training.series.only.series <- ifelse(dge$samples$Validation.series==1,0,1)
training.series.only.samples = rownames(dge$samples)[dge$samples$training.series.only.series==1]
dge <- filter.for.platelet.transcriptome.TrainEval.group(dge=dge, 
                                                         minimum.read.counts = 30,
                                                         minimum.prct.cohort = 90,
                                                         training.series.only = TRUE,
                                                         training.series.only.samples = training.series.only.samples,
                                                         groupPlateletomeSelection = TRUE,
                                                         verbose = TRUE)
ncol(dge$counts)
nrow(dge$counts)

###############################################################################################Thromboseq QC for Training and Evaluation Cohort
dge_TrainingEvaluation = dge[,c(colnames(dge)[dge$samples$Training.series==1],
                                colnames(dge)[dge$samples$Evaluation.series==1])]
min.number.reads.per.RNAs.detected = 0
min.number.total.RNAs.detected = 1500
k.variables = 2
variable.to.assess = c("lib.size")
variable.threshold = c(0.8)
ruvg.pvalue.threshold.group = 1e-2
ruvg.pvalue.threshold.strongest.variable = 1e-2
training.series.only = TRUE
training.series.only.samples = training.series.only.samples
leave.sample.out.threshold = 0.5
figureDir = "figureOutputFolder"
number.cores = 2
verbose = FALSE

registerDoMC(cores = number.cores)
detection.loop <- foreach(i = 1 : ncol(dge_TrainingEvaluation$counts)) %dopar% {
  # calculate number of genes with more than zero raw read counts
  # in case all genes are detected, summary would have only two output
  # options and is.na would result in FALSE, hence the if-statement
  sample.raw.counts <- dge_TrainingEvaluation$counts[, colnames(dge_TrainingEvaluation$counts)[i]]
  if (as.character(is.na(summary(as.numeric(sample.raw.counts) > 0)[3]))){
    sample.RNAs.detected <- as.numeric(summary(as.numeric(sample.raw.counts) > 
                                                 min.number.reads.per.RNAs.detected)[2])  
  } else {
    sample.RNAs.detected <- as.numeric(summary(as.numeric(sample.raw.counts) > 
                                                 min.number.reads.per.RNAs.detected)[3])  
  }
  
  # store data in container
  cont <- list()
  cont[["sample"]] <- colnames(dge_TrainingEvaluation$counts)[i] # collect sample ID
  cont[["sample.libsize"]] <- dge_TrainingEvaluation$samples[colnames(dge_TrainingEvaluation$counts)[i], "lib.size"] # collect lib.size
  cont[["sample.RNAs.detected"]] <- sample.RNAs.detected
  cont
}

# summarize data from loop into a data.frame
detection.data.frame <- data.frame(
  sample = unlist(lapply(detection.loop, function(x){x[["sample"]]})),
  lib.size = unlist(lapply(detection.loop, function(x){x[["sample.libsize"]]})),
  RNAs.detected = unlist(lapply(detection.loop, function(x){x[["sample.RNAs.detected"]]}))
)

# select samples that have to be excluded because of too little RNAs detected
samples.excluded.RNAs.detected <- detection.data.frame[detection.data.frame$RNAs.detected < 
                                                         min.number.total.RNAs.detected, ]
# update DGEList, and remove excluded samples
dgeIncludedSamples <- dge_TrainingEvaluation[, !colnames(dge_TrainingEvaluation) %in% samples.excluded.RNAs.detected$sample]
dgeIncludedSamples$samples <- droplevels(dgeIncludedSamples$samples)
ncol(dgeIncludedSamples)

dge_TrainingEvaluation <- perform.RUVg.correction(dge = dge[,colnames(dge)%in%colnames(dgeIncludedSamples)],
                                                  k.variables = 2, 
                                                  variable.to.assess =  c("lib.size"),
                                                  variable.threshold = c(0.8), 
                                                  ruvg.pvalue.threshold.group = 1e-2,
                                                  ruvg.pvalue.threshold.strongest.variable = 1e-2,
                                                  training.series.only = FALSE)
ncol(dge_TrainingEvaluation$counts)
nrow(dge_TrainingEvaluation$counts)
# replace raw counts by RUV-corrected counts, and recalculate lib-size
dge_TrainingEvaluation$counts <- dge_TrainingEvaluation$ruv.counts
dge_TrainingEvaluation$samples$lib.size <- colSums(dge_TrainingEvaluation$counts)
# perform TMM normalisation
dge_TrainingEvaluation <- calcNormFactorsThromboseq(dge_TrainingEvaluation,
                                                    normalize.on.training.series = FALSE, 
                                                    ref.sample.readout = FALSE) # calculate normalization factors
dge_TrainingEvaluation$samples <- droplevels(dge_TrainingEvaluation$samples)
# calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
normalized.counts <- cpm(dge_TrainingEvaluation, log = T, normalized.lib.sizes = T) 
nrow(normalized.counts)
ncol(normalized.counts)
data_1 <- normalized.counts

#直接进行ENSEMBLE号和基因名之间的转换
ENSGid<-rownames(data_1)
GENE_SYMBOL<- bitr(ENSGid,
                   fromType = "ENSEMBL",
                   toType = c("SYMBOL", "ENTREZID"),
                   OrgDb = org.Hs.eg.db)%>%
  distinct()%>%
  add_count(ENSEMBL)%>%
  filter(n==1)
data_2 <- data.frame(ENSGid)
data_2 <- merge(data_2,GENE_SYMBOL,by.x="ENSGid",by.y="ENSEMBL",all.x=T)
data_2$test <- ENSGid
data_2$test1 <- ifelse(data_2$test==data_2$ENSGid,0,1)
sum(data_2$test1)
data_2$SYMBOL <- ifelse(is.na(data_2$SYMBOL)==T,data_2$ENSGid,data_2$SYMBOL)
rownames(data_1) <- data_2$SYMBOL
data_1 <- t(data_1)
data_Train <- data_1[rownames(data_1)%in%c(colnames(dge)[dge$samples$Training.series==1]),]
nrow(data_Train)
data.Evaluation <- data_1[rownames(data_1)%in%c(colnames(dge)[dge$samples$Evaluation.series==1]),]
nrow(data.Evaluation)
View(data_Train)
View(data.Evaluation)
View(data.validation)
write.csv(data_Train,"Train.csv")
write.csv(data.Evaluation,"Evaluation.csv")
# correlate all samples to all other samples in a leave-one-sample-out cross-correlation analysis
# summarize data and filter those with Pearson's correlation below provided threshold
##数据集已经过LOOV 无洗脱
all.samples <- colnames(normalized.counts)
registerDoMC(cores = number.cores)
tmp.loop <- c()
for (i in 1:length(all.samples)){
  test.sample = all.samples[i]
  a <- cor(apply(normalized.counts[, all.samples[all.samples != test.sample]], 1, median), normalized.counts[, test.sample])
  tmp.loop <- c(tmp.loop,a)
}
leave.sample.out.output <- cbind(all.samples, tmp.loop)
leave.sample.out.output <- leave.sample.out.output[as.numeric(leave.sample.out.output[, 2])>0.5, ]
dge_TrainingEvaluation <-  dge_TrainingEvaluation[,colnames(dge_TrainingEvaluation)%in% leave.sample.out.output[, 1]]
ncol(dge_TrainingEvaluation$counts)
nrow(dge_TrainingEvaluation$counts)

###############################################################################################Thromboseq QC for Validation Cohort
dge_Validation = dge[,c(colnames(dge)[dge$samples$Validation.series==1])]
min.number.reads.per.RNAs.detected = 0
min.number.total.RNAs.detected = 1500
k.variables = 2
variable.to.assess = c("lib.size")
variable.threshold = c(0.8)
ruvg.pvalue.threshold.group = 1e-2
ruvg.pvalue.threshold.strongest.variable = 1e-2
training.series.only = TRUE
training.series.only.samples = training.series.only.samples
leave.sample.out.threshold = 0.5
figureDir = "figureOutputFolder"
number.cores = 2
verbose = FALSE

registerDoMC(cores = number.cores)
detection.loop <- foreach(i = 1 : ncol(dge_Validation$counts)) %dopar% {
  # calculate number of genes with more than zero raw read counts
  # in case all genes are detected, summary would have only two output
  # options and is.na would result in FALSE, hence the if-statement
  sample.raw.counts <- dge_Validation$counts[, colnames(dge_Validation$counts)[i]]
  if (as.character(is.na(summary(as.numeric(sample.raw.counts) > 0)[3]))){
    sample.RNAs.detected <- as.numeric(summary(as.numeric(sample.raw.counts) > 
                                                 min.number.reads.per.RNAs.detected)[2])  
  } else {
    sample.RNAs.detected <- as.numeric(summary(as.numeric(sample.raw.counts) > 
                                                 min.number.reads.per.RNAs.detected)[3])  
  }
  
  # store data in container
  cont <- list()
  cont[["sample"]] <- colnames(dge_Validation$counts)[i] # collect sample ID
  cont[["sample.libsize"]] <- dge_Validation$samples[colnames(dge_Validation$counts)[i], "lib.size"] # collect lib.size
  cont[["sample.RNAs.detected"]] <- sample.RNAs.detected
  cont
}

# summarize data from loop into a data.frame
detection.data.frame <- data.frame(
  sample = unlist(lapply(detection.loop, function(x){x[["sample"]]})),
  lib.size = unlist(lapply(detection.loop, function(x){x[["sample.libsize"]]})),
  RNAs.detected = unlist(lapply(detection.loop, function(x){x[["sample.RNAs.detected"]]}))
)

# select samples that have to be excluded because of too little RNAs detected
samples.excluded.RNAs.detected <- detection.data.frame[detection.data.frame$RNAs.detected < 
                                                         min.number.total.RNAs.detected, ]
# update DGEList, and remove excluded samples
dgeIncludedSamples <- dge_Validation[, !colnames(dge_Validation) %in% samples.excluded.RNAs.detected$sample]
dgeIncludedSamples$samples <- droplevels(dgeIncludedSamples$samples)
ncol(dgeIncludedSamples) #All in dge_validation

dge_Validation <- perform.RUVg.correction(dge =dge[,colnames(dge)%in%colnames(dgeIncludedSamples)],
                                          k.variables = 2, 
                                          variable.to.assess =  c("lib.size"),
                                          variable.threshold = c(0.8), 
                                          ruvg.pvalue.threshold.group = 1e-2,
                                          ruvg.pvalue.threshold.strongest.variable = 1e-2,
                                          training.series.only = FALSE)
# replace raw counts by RUV-corrected counts, and recalculate lib-size
dge_Validation$counts <- dge_Validation$ruv.counts
dge_Validation$samples$lib.size <- colSums(dge_Validation$counts)
# perform TMM normalisation
dge_Validation <- calcNormFactorsThromboseq(dge_Validation,
                                            normalize.on.training.series = FALSE, 
                                            ref.sample.readout = FALSE) # calculate normalization factors
dge_Validation$samples <- droplevels(dge_Validation$samples)
# calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
normalized.counts <- cpm(dge_Validation, log = T, normalized.lib.sizes = T) 
nrow(normalized.counts)
ncol(normalized.counts)
data_1 <- normalized.counts

#直接进行ENSEMBLE号和基因名之间的转换
ENSGid<-rownames(data_1)
GENE_SYMBOL<- bitr(ENSGid,
                   fromType = "ENSEMBL",
                   toType = c("SYMBOL", "ENTREZID"),
                   OrgDb = org.Hs.eg.db)%>%
  distinct()%>%
  add_count(ENSEMBL)%>%
  filter(n==1)
data_2 <- data.frame(ENSGid)
data_2 <- merge(data_2,GENE_SYMBOL,by.x="ENSGid",by.y="ENSEMBL",all.x=T)
data_2$test <- ENSGid
data_2$test1 <- ifelse(data_2$test==data_2$ENSGid,0,1)
sum(data_2$test1)
data_2$SYMBOL <- ifelse(is.na(data_2$SYMBOL)==T,data_2$ENSGid,data_2$SYMBOL)
rownames(data_1) <- data_2$SYMBOL
data_1 <- t(data_1)
write.csv(data_1,"Validation.csv")

#数据集已经过LOOV 无洗脱
all.samples <- colnames(normalized.counts)
registerDoMC(cores = number.cores)
tmp.loop <- c()
for (i in 1:length(all.samples)){
  test.sample = all.samples[i]
  a <- cor(apply(normalized.counts[, all.samples[all.samples != test.sample]], 1, median), normalized.counts[, test.sample])
  tmp.loop <- c(tmp.loop,a)
}
leave.sample.out.output <- cbind(all.samples, tmp.loop)
leave.sample.out.output <- leave.sample.out.output[as.numeric(leave.sample.out.output[, 2])>0.5, ]
dge_Validation <-  dge_Validation[,colnames(dge_Validation)%in% leave.sample.out.output[, 1]]

#########################################
#Supplementary Table 1-3
setwd("D:/高端paper work/AI/血管/AI结合循环血小板mRNA预测隐匿性肿瘤/前期/正式图表/投稿使用/Fundamental Research/InTVeld_Pancancer_TSOO-main")
# load and prepare sample info file provided in GEO. Identical to Table S2 of the manuscript.
sampleInfo <- read.csv('TableS2.csv', sep=",", row.names=1)
sampleInfo$rownames <- rownames(sampleInfo)
setwd("D:/高端paper work/AI/血管/AI结合循环血小板mRNA预测隐匿性肿瘤/前期/正式图表/投稿使用/返修/cluster3 4")
data_derivation_clusters <- read.csv("output_train_3cluster.csv")
data_derivation_clusters <- read.csv("output_test_3cluster.csv")
data_derivation_clusters <- read.csv("output_validation_3cluster.csv")
data_derivation_clusters <- data_derivation_clusters[,c("Unnamed..0","X0")]
colnames(data_derivation_clusters) <- c("patients","cluster")
data_derivation_clusters$cluster <- ifelse(data_derivation_clusters$cluster=="0","cluster 1",
                                           ifelse(data_derivation_clusters$cluster=="1","cluster 2",
                                                  ifelse(data_derivation_clusters$cluster=="2","cluster 3",
                                                         ifelse(data_derivation_clusters$cluster=="3","cluster 4","cluster 5"))))
data_derivation <- merge(data_derivation_clusters,sampleInfo,by.x="patients",by.y="rownames",all.x=T)

data <- data_derivation
data$Stage <- ifelse(data$Stage=="n.a.",NA,as.character(data$Stage))
data$Sex <- ifelse(data$Sex=="n.a.",NA,as.character(data$Sex))
data$Age <- ifelse(data$Age=="n.a.",NA,data$Age)
data$Age  <- as.numeric(data$Age)
data$`BRCA` <- ifelse(data$Group=="Breast cancer",1,0)
data$`CHOL` <- ifelse(data$Group=="Cholangiocarcinoma",1,0)
data$`CRC` <- ifelse(data$Group=="Colorectal cancer",1,0)
data$`ENDO` <- ifelse(data$Group=="Endometrial cancer",1,0)
data$`ESO` <- ifelse(data$Group=="Esophageal carcinoma",1,0)
data$`GLIO` <- ifelse(data$Group=="Glioma",1,0)
data$`HNSSC` <- ifelse(data$Group=="Head and neck cancer",1,0)
data$`HCC` <- ifelse(data$Group=="Hepatocellular carcinoma",1,0)
data$`LYM` <- ifelse(data$Group=="Lymphoma",1,0)
data$`MELA` <- ifelse(data$Group=="Melanoma",1,0)
data$`MM` <- ifelse(data$Group=="Multiple Myeloma",1,0)
data$`NSCLC` <- ifelse(data$Group=="Non-small-cell lung cancer",1,0)
data$`OVCA` <- ifelse(data$Group=="Ovarian cancer",1,0)
data$`PDAC` <- ifelse(data$Group=="Pancreatic cancer",1,0)
data$`PRCA` <- ifelse(data$Group=="Prostate cancer",1,0)
data$`RCC` <- ifelse(data$Group=="Renal cell cancer",1,0)
data$`SARC` <- ifelse(data$Group=="Sarcoma",1,0)
data$`URO` <- ifelse(data$Group=="Urothelial cancer",1,0)
data$`I Stage` <- ifelse(is.na(data$Stage)==T,NA,ifelse(data$Stage=="I",1,0))
data$`II Stage` <- ifelse(is.na(data$Stage)==T,NA,ifelse(data$Stage=="II",1,0))
data$`III Stage` <- ifelse(is.na(data$Stage)==T,NA,ifelse(data$Stage=="III",1,0))
data$`IV Stage` <-ifelse(is.na(data$Stage)==T,NA,ifelse(data$Stage=="IV",1,0))

###输出临床信息
data <- data[,c("patients","Group","Stage","Sex","Age")]
write.csv(data,"Derivation_Clinical_Info.csv")
write.csv(data,"Evaluation_Clinical_Info.csv")
write.csv(data,"Validation_Clinical_Info.csv")

###sTable 1
varsToFactor <- c("Sex","BRCA","CHOL","CRC","ENDO","ESO","GLIO","HNSSC","HCC","LYM","MELA","MM","NSCLC","OVCA","PDAC","PRCA","RCC","SARC","URO","I Stage","II Stage","III Stage","IV Stage")
vars <- c("Sex","Age","BRCA","CHOL","CRC","ENDO","ESO","GLIO","HNSSC","HCC","LYM","MELA","MM","NSCLC","OVCA","PDAC","PRCA","RCC","SARC","URO","I Stage","II Stage","III Stage","IV Stage")
data[varsToFactor] <- lapply(data[varsToFactor],factor)
tableOne <- CreateTableOne(vars = vars,strata = 'cluster', data = data)
tableOne <-print(tableOne,nonnormal = c("Age"))
write.csv(tableOne,"sTable 1.csv")
write.csv(tableOne,"sTable 2.csv")
write.csv(tableOne,"sTable 3.csv")

#############Supplementary Table 4-6
setwd("D:/高端paper work/AI/血管/AI结合循环血小板mRNA预测隐匿性肿瘤/前期/正式图表/投稿使用/Fundamental Research/InTVeld_Pancancer_TSOO-main")
# load and prepare sample info file provided in GEO. Identical to Table S2 of the manuscript.
sampleInfo <- read.csv('TableS2.csv', sep=",", row.names=1)
sampleInfo$rownames <- rownames(sampleInfo)
setwd("D:/高端paper work/AI/血管/AI结合循环血小板mRNA预测隐匿性肿瘤/前期/正式图表/投稿使用/返修/cluster3 4")
##
data_derivation_clusters <- read.csv("output_train_3cluster.csv")
##
data_derivation_clusters <- read.csv("output_test_3cluster.csv")
##
data_derivation_clusters <- read.csv("output_validation_3cluster.csv")
##
data_derivation_clusters <- data_derivation_clusters[,c("Unnamed..0","X0")]
colnames(data_derivation_clusters) <- c("patients","cluster")
data_derivation_clusters$cluster <- ifelse(data_derivation_clusters$cluster=="0","cluster 1",
                                           ifelse(data_derivation_clusters$cluster=="1","cluster 2",
                                                  ifelse(data_derivation_clusters$cluster=="2","cluster 3",
                                                         ifelse(data_derivation_clusters$cluster=="3","cluster 4","cluster 5"))))
data_derivation <- merge(data_derivation_clusters,sampleInfo,by.x="patients",by.y="rownames",all.x=T)
data <- data_derivation
data$Stage <- ifelse(data$Stage=="n.a.",NA,as.character(data$Stage))
data$Sex <- ifelse(data$Sex=="n.a.",NA,as.character(data$Sex))
data$Age <- ifelse(data$Age=="n.a.",NA,data$Age)
data$Age  <- as.numeric(data$Age)
data$`I Stage` <- ifelse(is.na(data$Stage)==T,NA,ifelse(data$Stage=="I",1,0))
data$`II Stage` <- ifelse(is.na(data$Stage)==T,NA,ifelse(data$Stage=="II",1,0))
data$`III Stage` <- ifelse(is.na(data$Stage)==T,NA,ifelse(data$Stage=="III",1,0))
data$`IV Stage` <-ifelse(is.na(data$Stage)==T,NA,ifelse(data$Stage=="IV",1,0))
setwd("D:/高端paper work/AI/血管/AI结合循环血小板mRNA预测隐匿性肿瘤/前期/正式图表/投稿使用/返修/cluster3 4/Stable 4")
##
b <- "Derivation_"
##
b <- "Evaluation_"
##
b <- "Validation_"
for (i in 1:18){
  a <- sort(as.character(c(unique(data$Group))))[i]
  data1 <- subset(data,data$Group==a)
  data1 <- subset(data1,is.na(data1$`I Stage`)==F)
  if(nrow(data1>0)){
    if(length(sort(as.character(c(unique(data1$cluster)))))==1){
    }else{
      varsToFactor <- c("I Stage","II Stage","III Stage","IV Stage")
      vars <- c("I Stage","II Stage","III Stage","IV Stage")
      data1[varsToFactor] <- lapply(data1[varsToFactor],factor)
      tableOne <- CreateTableOne(vars = vars,strata = 'cluster', data = data1)
      tableOne <-print(tableOne)
      tableOne <- data.frame(tableOne)
      tableOne$Cancer = a
      write.csv(tableOne,paste(b,a,".csv",sep=""))
    }
  }
}

#############Supplementary Table 7
setwd("D:/高端paper work/AI/血管/AI结合循环血小板mRNA预测隐匿性肿瘤/前期/正式图表/投稿使用/Fundamental Research/InTVeld_Pancancer_TSOO-main")
# load and prepare sample info file provided in GEO. Identical to Table S2 of the manuscript.
sampleInfo <- read.csv('TableS2.csv', sep=",", row.names=1)
sampleInfo$rownames <- rownames(sampleInfo)
data_final <- data.frame()

setwd("D:/高端paper work/AI/血管/AI结合循环血小板mRNA预测隐匿性肿瘤/前期/正式图表/投稿使用/返修/cluster3 4")
data_derivation_clusters <- read.csv("output_train_3cluster.csv")
data_derivation_clusters <- read.csv("output_test_3cluster.csv")
data_derivation_clusters <- read.csv("output_validation_3cluster.csv")

data_derivation_clusters <- data_derivation_clusters[,c("Unnamed..0","X0")]
colnames(data_derivation_clusters) <- c("patients","cluster")
data_derivation_clusters$cluster <- ifelse(data_derivation_clusters$cluster=="0","cluster 1",
                                           ifelse(data_derivation_clusters$cluster=="1","cluster 2",
                                                  ifelse(data_derivation_clusters$cluster=="2","cluster 3",
                                                         ifelse(data_derivation_clusters$cluster=="3","cluster 4","cluster 5"))))
data_derivation <- merge(data_derivation_clusters,sampleInfo,by.x="patients",by.y="rownames",all.x=T)
cancers <- sort(as.character(c(unique(data_derivation$Group))))
cancers <- c(cancers,"I","II","III","IV")
##
data <- data_derivation
data$Stage <- ifelse(data$Stage=="n.a.",NA,as.character(data$Stage))
data$Sex <- ifelse(data$Sex=="n.a.",NA,as.character(data$Sex))
data$Age <- ifelse(data$Age=="n.a.",NA,data$Age)
data$Age  <- as.numeric(data$Age)


data_a <- data.frame()
data_b <- data.frame()
data_d <- data.frame()
##
CandidateVariables <- c("cluster")
##
CandidateVariables <- c("Age","Sex","cluster")
##

for (i in 1:length(cancers)){
  data_xx <- data
  z <- cancers[i]
  if(i <= 18){
    data_xx$outcome <- ifelse(data_xx$Group==z,1,0)
  }else{
    data_xx$outcome <- ifelse(is.na(data_xx$Stage)==T,NA,ifelse(data_xx$Stage==z,1,0))
  }
  
  Outcome <- "outcome"
  data_xx <- subset(data_xx,is.na(data_xx$outcome)==F)
  Formula <- formula(paste(paste(Outcome,"~", collapse=" "), 
                           paste(CandidateVariables, collapse=" + ")))
  model.lasso <- glm(Formula,data=data_xx,family=binomial)
  a <- summary(model.lasso)$coefficients
  a <- data.frame(a)
  a$cancer = z
  b <- exp(model.lasso$coefficients)
  b <- data.frame(b)
  b$cancer = z
  d <- exp(confint(model.lasso))
  d <- data.frame(d)
  d$cancer = z
  data_a <- rbind(data_a,a)
  data_b <- rbind(data_b,b)
  data_d <- rbind(data_d,d)
}
data_1 <- cbind(data_a,data_b,data_d)
data_1 <- data_1[,c("cancer","b","X2.5..","X97.5..","Pr...z..")]
##
data_1 <- subset(data_1,data_1$Pr...z..<0.05)
##
data_1$rownames <- rownames(data_1)
data_2 <- data_1[grepl("cluster", data_1$rownames),]

data_2$adjusted <- 0
data_2$adjusted <- 1
##
data_2$cohort <- "Training"
data_2$cohort <- "Evaluation"
data_2$cohort <- "Validation"
data_final <- rbind(data_final,data_2)
write.csv(data_final,"sTable 71.csv")

#Figure 1 Association with Clusters and Cancer Type and Stage in Derivation cohort(弦图)
#肿瘤类型弦图 以肿瘤为上标
setwd("D:/高端paper work/AI/血管/AI结合循环血小板mRNA预测隐匿性肿瘤/前期/正式图表/投稿使用/Fundamental Research/InTVeld_Pancancer_TSOO-main")
# load and prepare sample info file provided in GEO. Identical to Table S2 of the manuscript.
sampleInfo <- read.csv('TableS2.csv', sep=",", row.names=1)
sampleInfo$rownames <- rownames(sampleInfo)
setwd("D:/高端paper work/AI/血管/AI结合循环血小板mRNA预测隐匿性肿瘤/前期/正式图表/投稿使用/返修/cluster3 4")
##
data_derivation_clusters <- read.csv("output_train_3cluster.csv")
##
data_derivation_clusters <- read.csv("output_test_3cluster.csv")
##
data_derivation_clusters <- read.csv("output_validation_3cluster.csv")

data_derivation_clusters <- data_derivation_clusters[,c("Unnamed..0","X0")]
colnames(data_derivation_clusters) <- c("patients","cluster")
data_derivation_clusters$cluster <- ifelse(data_derivation_clusters$cluster=="0","cluster 1",
                                           ifelse(data_derivation_clusters$cluster=="1","cluster 2",
                                                  ifelse(data_derivation_clusters$cluster=="2","cluster 3",
                                                         ifelse(data_derivation_clusters$cluster=="3","cluster 4","cluster 5"))))
data_derivation <- merge(data_derivation_clusters,sampleInfo,by.x="patients",by.y="rownames",all.x=T)

cancers = sort(as.character(c(unique(data_derivation$Group))))
clusters = sort(as.character(c(unique(data_derivation$cluster))))
data_final <- data.frame(clusters)
data_final$`BRCA` <- rep(0,nrow(data_final))
data_final$`CHOL` <- rep(0,nrow(data_final))
data_final$`CRC` <- rep(0,nrow(data_final))
data_final$`ENDO` <- rep(0,nrow(data_final))
data_final$`ESO` <- rep(0,nrow(data_final))
data_final$`GLIO` <- rep(0,nrow(data_final))
data_final$`HNSSC` <- rep(0,nrow(data_final))
data_final$`HCC` <- rep(0,nrow(data_final))
data_final$`LYM` <- rep(0,nrow(data_final))
data_final$`MELA` <- rep(0,nrow(data_final))
data_final$`MM` <- rep(0,nrow(data_final))
data_final$`NSCLC` <- rep(0,nrow(data_final))
data_final$`OVCA` <- rep(0,nrow(data_final))
data_final$`PDAC` <- rep(0,nrow(data_final))
data_final$`PRCA` <- rep(0,nrow(data_final))
data_final$`RCC` <- rep(0,nrow(data_final))
data_final$`SARC` <- rep(0,nrow(data_final))
data_final$`URO` <- rep(0,nrow(data_final))
for (i in 1:length(cancers)){
  a <- cancers[i]
  data_1 <- subset(data_derivation,data_derivation$Group==a)
  for (z in 1:length(clusters)){
    b <- clusters[z]
    data_2 <- subset(data_1,data_1$cluster==b)
    c <- round(nrow(data_2)/nrow(data_1),8)*100
    data_final[z,i+1] <- c
  }
}
for (i in 1:length(cancers)){
  data_final[3,i+1] <- 100-sum(data_final[1:2,i+1])
}
rownames(data_final) <- data_final$clusters
data_final <- data_final[,-1]
data_final = as.matrix(data_final)
chordDiagram(data_final)

#肿瘤类型弦图 以分期为上标
setwd("D:/高端paper work/AI/血管/AI结合循环血小板mRNA预测隐匿性肿瘤/前期/正式图表/投稿使用/Fundamental Research/InTVeld_Pancancer_TSOO-main")
# load and prepare sample info file provided in GEO. Identical to Table S2 of the manuscript.
sampleInfo <- read.csv('TableS2.csv', sep=",", row.names=1)
sampleInfo$rownames <- rownames(sampleInfo)
setwd("D:/高端paper work/AI/血管/AI结合循环血小板mRNA预测隐匿性肿瘤/前期/正式图表/投稿使用/返修/cluster3 4")
##
data_derivation_clusters <- read.csv("output_train_3cluster.csv")
##
data_derivation_clusters <- read.csv("output_test_3cluster.csv")
##
data_derivation_clusters <- read.csv("output_validation_3cluster.csv")

data_derivation_clusters <- data_derivation_clusters[,c("Unnamed..0","X0")]
colnames(data_derivation_clusters) <- c("patients","cluster")
data_derivation_clusters$cluster <- ifelse(data_derivation_clusters$cluster=="0","cluster 1",
                                           ifelse(data_derivation_clusters$cluster=="1","cluster 2",
                                                  ifelse(data_derivation_clusters$cluster=="2","cluster 3",
                                                         ifelse(data_derivation_clusters$cluster=="3","cluster 4","cluster 5"))))
data_derivation <- merge(data_derivation_clusters,sampleInfo,by.x="patients",by.y="rownames",all.x=T)


Stage = sort(as.character(c(unique(subset(data_derivation$Stage,data_derivation$Stage!="n.a.")))))
clusters = sort(as.character(c(unique(data_derivation$cluster))))
data_final <- data.frame(clusters)
data_final$`I Stage` <- rep(0,nrow(data_final))
data_final$`II Stage` <- rep(0,nrow(data_final))
data_final$`III Stage` <- rep(0,nrow(data_final))
data_final$`IV Stage` <- rep(0,nrow(data_final))

for (i in 1:length(Stage)){
  a <- Stage[i]
  data_1 <- subset(data_derivation,data_derivation$Stage==a)
  for (z in 1:length(clusters)){
    b <- clusters[z]
    data_2 <- subset(data_1,data_1$cluster==b)
    c <- round(nrow(data_2)/nrow(data_1),8)*100
    data_final[z,i+1] <- c
  }
}
for (i in 1:length(Stage)){
  data_final[3,i+1] <- 100-sum(data_final[1:2,i+1])
}
rownames(data_final) <- data_final$clusters
data_final <- data_final[,-1]
data_final = as.matrix(data_final)
chordDiagram(data_final)

###########################################################################################################heatmap
data_derivation <- read.csv("output_train_3cluster.csv")

data_train <- data_derivation
rownames(data_train)<-data_train$Unnamed..0
data_train$cluster <- ifelse(data_train$X0=="0","cluster 1",
                             ifelse(data_train$X0=="1","cluster 2",
                                    ifelse(data_train$X0=="2","cluster 3",
                                           ifelse(data_train$X0=="3","cluster 4","cluster 5"))))

gene_list<-c('PCNX4', 'FKBP5', 'ITGB3BP', 'RECQL', 'FAM107B', 'ITGA4', 'KSR1', 'AP1B1', 'NELL2', 'HVCN1', 
             'TUBA1C', 'ZNF346', 'HTT', 'CRYM', 'ALAS2', 'ACADVL', 'EIF4G1', 'PRKAB2', 'ALDH16A1', 'ZNF542P',
             'RAB5B', 'ATP5MC1', 'WFDC1', 'STRADB', 'UPP1', 'EEF1B2', 'TEX9', 'TBC1D14', 'MAOB', 'ARL2')

data_train<-data_train[,c(gene_list,"cluster")]

cluster1<-subset(data_train,data_train$cluster=="cluster 1")
cluster1 <- cluster1[,-31]
cluster1_d<-colMeans(cluster1)
cluster2<-subset(data_train,data_train$cluster=="cluster 2")
cluster2 <- cluster2[,-31]
cluster2_d<-colMeans(cluster2)
cluster3<-subset(data_train,data_train$cluster=="cluster 3")
cluster3 <- cluster3[,-31]
cluster3_d<-colMeans(cluster3)


data_evaluation <- read.csv("output_test_3cluster_evaluation.csv")
rownames(data_evaluation)<-data_evaluation$Unnamed..0
data_evaluation$cluster <- ifelse(data_evaluation$X0=="0","cluster 1",
                                  ifelse(data_evaluation$X0=="1","cluster 2",
                                         ifelse(data_evaluation$X0=="2","cluster 3",
                                                ifelse(data_evaluation$X0=="3","cluster 4","cluster 5"))))


gene_list<-c('PCNX4', 'FKBP5', 'ITGB3BP', 'RECQL', 'FAM107B', 'ITGA4', 'KSR1', 'AP1B1', 'NELL2', 'HVCN1', 
             'TUBA1C', 'ZNF346', 'HTT', 'CRYM', 'ALAS2', 'ACADVL', 'EIF4G1', 'PRKAB2', 'ALDH16A1', 'ZNF542P',
             'RAB5B', 'ATP5MC1', 'WFDC1', 'STRADB', 'UPP1', 'EEF1B2', 'TEX9', 'TBC1D14', 'MAOB', 'ARL2')

data_evaluation<-data_evaluation[,c(gene_list,"cluster")]

cluster1<-subset(data_evaluation,data_evaluation$cluster=="cluster 1")
cluster1 <- cluster1[,-31]
cluster1_e<-colMeans(cluster1)
cluster2<-subset(data_evaluation,data_evaluation$cluster=="cluster 2")
cluster2 <- cluster2[,-31]
cluster2_e<-colMeans(cluster2)
cluster3<-subset(data_evaluation,data_evaluation$cluster=="cluster 3")
cluster3 <- cluster3[,-31]
cluster3_e<-colMeans(cluster3)


data_validation <- read.csv("output_validation_3cluster_evaluation.csv")

data_test<- data_validation
rownames(data_test)<-data_test$Unnamed..0
data_test$cluster <- ifelse(data_test$X0=="0","cluster 1",
                            ifelse(data_test$X0=="1","cluster 2",
                                   ifelse(data_test$X0=="2","cluster 3",
                                          ifelse(data_test$X0=="3","cluster 4","cluster 5"))))


gene_list<-c('PCNX4', 'FKBP5', 'ITGB3BP', 'RECQL', 'FAM107B', 'ITGA4', 'KSR1', 'AP1B1', 'NELL2', 'HVCN1', 
             'TUBA1C', 'ZNF346', 'HTT', 'CRYM', 'ALAS2', 'ACADVL', 'EIF4G1', 'PRKAB2', 'ALDH16A1', 'ZNF542P',
             'RAB5B', 'ATP5MC1', 'WFDC1', 'STRADB', 'UPP1', 'EEF1B2', 'TEX9', 'TBC1D14', 'MAOB', 'ARL2')

data_test<-data_test[,c(gene_list,"cluster")]

cluster1<-subset(data_test,data_test$cluster=="cluster 1")
cluster1 <- cluster1[,-31]
cluster1_v<-colMeans(cluster1)
cluster2<-subset(data_test,data_test$cluster=="cluster 2")
cluster2 <- cluster2[,-31]
cluster2_v<-colMeans(cluster2)
cluster3<-subset(data_test,data_test$cluster=="cluster 3")
cluster3 <- cluster3[,-31]
cluster3_v<-colMeans(cluster3)


df<-rbind(cluster1_d,cluster2_d,cluster3_d,cluster1_e,cluster2_e,cluster3_e,cluster1_v,cluster2_v,cluster3_v)
df<-log2(df+1)
df<-as.data.frame(df)
df<-t(df)

Groups=c(rep("Derivation",3),rep("Evaluation",3),rep("Validation",3))

annotation_c<-data.frame(Groups)

rownames(annotation_c)<-colnames(df)


library(pheatmap)
library(RColorBrewer)
pheatmap::pheatmap(df, 
                   
                   border=F,
                   cluster_cols = F,
                   fontsize_row = 8,
                   height = 5,
                   width=8,
                   annotation_col = annotation_c,
                   show_colnames = T,
                   angle_col = "45",
                   show_rownames = T,
                   annotation_colors=list(Groups=c(Derivation="#91CAB6",Evaluation="light yellow",Validation="#A8DDE0")),
                   color = colorRampPalette(brewer.pal(n = 7, name ="GnBu" ))(100))

#火山图
#Figure 4 Functional Annotation of Each clusters in Derivation and Validation cohort (GSEA)
#火山图
##############################
##limma差异分析
data_derivation <- read.csv("output_train_3cluster.csv")
data_derivation <- read.csv("output_test_3cluster.csv")
data_derivation <- read.csv("output_validation_3cluster.csv")
library(limma)
dat<-data_derivation
dat<-dat[,-1:-2]
for (i in 1:5){
  dat<-dat[,-5441]
}
ncol(dat)
dat<-t(dat)
group_list<-ifelse(data_derivation$X0==0,"cluster1","other")
design <- model.matrix(~0+factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(dat)
head(design)
exp="cluster1"
ctr="other"
contrast.matrix <- makeContrasts(contrasts=paste0(exp,'-',ctr),
                                 levels = design)
fit1 <- lmFit(dat,design)                 #拟合模型
fit2 <- contrasts.fit(fit1, contrast.matrix) #统计检验
efit <- eBayes(fit2)                         #修正
summary(decideTests(efit,lfc=1, p.value=0.05)) #统计查看差异结果
tempOutput <- topTable(efit, coef=paste0(exp,'-',ctr), n=Inf)
degs_derivation_cluster1<- na.omit(tempOutput) 

group_list<-ifelse(data_derivation$X0==1,"cluster2","other")
design <- model.matrix(~0+factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(dat)
head(design)
exp="cluster2"
ctr="other"
contrast.matrix <- makeContrasts(contrasts=paste0(exp,'-',ctr),
                                 levels = design)
fit1 <- lmFit(dat,design)                 #拟合模型
fit2 <- contrasts.fit(fit1, contrast.matrix) #统计检验
efit <- eBayes(fit2)                         #修正
summary(decideTests(efit,lfc=1, p.value=0.05)) #统计查看差异结果
tempOutput <- topTable(efit, coef=paste0(exp,'-',ctr), n=Inf)
degs_derivation_cluster2<- na.omit(tempOutput) 

group_list<-ifelse(data_derivation$X0==2,"cluster3","other")
design <- model.matrix(~0+factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(dat)
head(design)
exp="cluster3"
ctr="other"
contrast.matrix <- makeContrasts(contrasts=paste0(exp,'-',ctr),
                                 levels = design)
fit1 <- lmFit(dat,design)                 #拟合模型
fit2 <- contrasts.fit(fit1, contrast.matrix) #统计检验
efit <- eBayes(fit2)                         #修正
summary(decideTests(efit,lfc=1, p.value=0.05)) #统计查看差异结果
tempOutput <- topTable(efit, coef=paste0(exp,'-',ctr), n=Inf)
degs_derivation_cluster3<- na.omit(tempOutput) 


degs_derivation_cluster1$group<-ifelse(degs_derivation_cluster1$logFC>0,"up","down")
degs_derivation_cluster1$cluster=1
degs_derivation_cluster1$symbol <- rownames(degs_derivation_cluster1)

degs_derivation_cluster2$group<-ifelse(degs_derivation_cluster2$logFC>0,"up","down")
degs_derivation_cluster2$cluster=2
degs_derivation_cluster2$symbol <- rownames(degs_derivation_cluster2)

degs_derivation_cluster3$group<-ifelse(degs_derivation_cluster3$logFC>0,"up","down")
degs_derivation_cluster3$cluster=3
degs_derivation_cluster3$symbol <- rownames(degs_derivation_cluster3)



df<-rbind(degs_derivation_cluster1,degs_derivation_cluster2,degs_derivation_cluster3)

# for循环挑选每个cluster的top前5 gene symbol
tm <- function(data){
  for (i in c(1:3)) {
    assign(paste("topgene",i,sep = ""),
           filter(data,cluster==i) %>% distinct(symbol,.keep_all = T) %>% top_n(5,abs(logFC)))
  }
  topgene <- rbind(topgene1,topgene2,topgene3)
  return(topgene)
}
top <- tm(df)

#相关包的载入
library(ggplot2)
library(ggrepel)
# 先画背景柱，根据数据log2FC的max值,min值来确定
#根据数据中log2FC区间确定背景柱长度：
col1<-data.frame(x=c(1,2,3),
                 y=c(max(degs_derivation_cluster1$logFC),max(degs_derivation_cluster2$logFC),max(degs_derivation_cluster3$logFC)))
col2<-data.frame(x=c(1,2,3),
                 y=c(min(degs_derivation_cluster1$logFC),min(degs_derivation_cluster2$logFC),min(degs_derivation_cluster3$logFC) ))
# 绘制背景柱
p1 <- ggplot()+
  geom_col(data = col1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = col2,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)
p1

#把散点火山图叠加到背景柱上：
p2 <- ggplot()+
  geom_col(data = col1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = col2,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data = df,
              aes(x =cluster , y = logFC, color = group),
              size = 1,
              width =0.4)+
  scale_color_manual(name=NULL,
                     values = c("#4393C3","#FC4E2A"))+
  labs(x="",y="log2(FoldChange)")

p2

# 添加X轴的分组色块标签：
dfcol<-data.frame(x=c(1:3),
                  y=0,
                  label=c(1:3))
# 添加分组色块标签
dfcol$group <- c("cluster 1","cluster 2","cluster 3")
# 加载包
library(RColorBrewer)
library(MetBrewer)
# 自定义分组色块的颜色
tile_color <- met.brewer("Thomas",3)
# 在图中镶嵌色块
p3 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=0.2,
                     color = "black",
                     fill = tile_color,
                     alpha = 1,
                     show.legend = F)+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=group),
            size =3,
            color ="white")
p3


# 添加基因标签
p4<-p3+geom_text_repel(
  data=top,
  aes(x=cluster,y=logFC,label=symbol),
  size=4,
  color="red",
  force = 1.2,
  arrow = arrow(length = unit(0.008, "npc"),
                type = "open", ends = "last"))
p4

# 去除背景，美化图片
p5 <- p4+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 10)
  )
p5



#差异基因按foldchange排序
degs <- degs_derivation_cluster1
degs <- degs_derivation_cluster2
degs <- degs_derivation_cluster3

genelist = degs$logFC
names(genelist) = rownames(degs)
genelist = sort(genelist,decreasing = T)
head(genelist)
library(clusterProfiler)
head(hallmark)
egmt <- GSEA(genelist, TERM2GENE = hallmark, 
             minGSSize = 10,maxGSSize = 10000,
             pvalueCutoff = 1,pAdjustMethod = 'BH')
gesa_res <- egmt@result
gesa_res <- gesa_res[gesa_res$pvalue < 0.05,]
gesa_res <- gesa_res[order(gesa_res$NES,decreasing = T),]
head(gesa_res)
genelist <- data.frame(genelist)
genelist$gene <- rownames(genelist)

genelist$gene <- as.character(genelist$gene)
data_1 <- genelist[genelist$gene%in% c("KLK3","FOLH1","NPY"),]
colnames(gesa_res)
library(enrichplot)
p1<-gseaNb(object = egmt,
           geneSetID = rownames(gesa_res)[1:3],lineSize=1,
           curveCol=c("blueviolet","brown1","pink2"))
p1

#Table 2 Important omics indics associated with cacner type and stage
setwd("D:/高端paper work/AI/血管/AI结合循环血小板mRNA预测隐匿性肿瘤/前期/正式图表/投稿使用/返修/cluster3 4")
data_feature_type <- read.csv("feature_type.csv") 
data_feature_type <- read.csv("feature_Stage.csv") 
data_feature_type <- read.csv("1+2.csv")
colnames(data_feature_type)
genelist <- data_feature_type$Breast.cancer
genelist <- data_feature_type$Cholangiocarcinoma
genelist <- data_feature_type$Colorectal.cancer
genelist <- data_feature_type$Endometrial.cancer
genelist <- data_feature_type$Esophageal.carcinoma ##
genelist <- data_feature_type$Glioma
genelist <- data_feature_type$Hepatocellular.carcinoma
genelist <- data_feature_type$Head.and.neck.cancer
genelist <- data_feature_type$Lymphoma ##
genelist <- data_feature_type$Melanoma ##
genelist <- data_feature_type$Multiple.Myeloma
genelist <- data_feature_type$Non.small.cell.lung.cancer
genelist <- data_feature_type$Ovarian.cancer
genelist <- data_feature_type$Pancreatic.cancer ##
genelist <- data_feature_type$Prostate.cancer
genelist <- data_feature_type$Renal.cell.cancer
genelist <- data_feature_type$Sarcoma
genelist <- data_feature_type$Urothelial.cancer
genelist <- data_feature_type$I ##
genelist <- data_feature_type$II ##
genelist <- data_feature_type$III ##
genelist <- data_feature_type$IV ##
genelist <- data_feature_type$I.II
#Top 3 genes
genelist[1:3]
genelist <- genelist[1:25]
#Functional Annoation
DEG_trans<- bitr(genelist,
                 fromType = "SYMBOL",
                 toType = c("ENSEMBL", "ENTREZID"),
                 OrgDb = org.Hs.eg.db)%>%
  filter(!duplicated(SYMBOL)) 

head(DEG_trans)
gene<-(DEG_trans$ENTREZID)

#KEGG富集分析
kk <- enrichKEGG(gene = gene,
                 organism = "hsa", 
                 pvalueCutoff =0.1, 
                 qvalueCutoff =0.1)
kk@result$Description[1:3] #输出KEGG最显著的

