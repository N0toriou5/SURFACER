# ML predictor
# https://machinelearningmastery.com/machine-learning-in-r-step-by-step/
setwd("F:/projects/Breast/")

load("results/000_spmra.rda")
tcga<-spmra
rm(spmra)
load("results/000_metabric_spmra.rda")
metabric<-spmra
rm(spmra)
common<-intersect(rownames(tcga),rownames(metabric))
tcga<-tcga[common,]
metabric<-metabric[common,]
load("surfacer_2020.rda")
library(caret)
library(kernlab)
library(randomForest)
library(arm)
#library(RSNNS)
#validation_index <- createDataPartition(dataset$Species, p=0.80, list=FALSE)
load("results/000_patient_mapping.rda")
my_patients$`Surface Genes`[my_patients$`Surface Genes`=="plum2"]<-"Basal-enriched"
my_patients$`Surface Genes`[my_patients$`Surface Genes`=="tomato"]<-"Lum1"
my_patients$`Surface Genes`[my_patients$`Surface Genes`=="cyan"]<-"Lum2"
my_patients$`Surface Genes`[my_patients$`Surface Genes`=="seagreen3"]<-"Lum3"
my_patients$`Surface Genes`[my_patients$`Surface Genes`=="darkkhaki"]<-"Mixed"

dataset<-as.data.frame(t(tcga))
identical(rownames(my_patients),rownames(dataset)) ### TRUE
dataset$subtype<-factor(my_patients$`Surface Genes`,levels<-c("Basal-enriched","Lum1","Lum2","Lum3","Mixed"))
# Run algorithms using 10-fold cross validation
control <- trainControl(method="cv", number=10)
metric <- "Accuracy"

#### Test several algorithms

# a) non-linear algorithms
set.seed(3)
fit.svm <- train(subtype~., data=dataset, method="svmRadial", metric=metric, trControl=control)
summary(fit.svm)

# estimate skill of LDA on the validation dataset
validation<-t(metabric)
ML_predictions <- predict(fit.svm, validation)
str(ML_predictions)
length(ML_predictions) #1949
met_srftyp<-data.frame(patients=rownames(validation),surftypes=as.character(ML_predictions))
rownames(met_srftyp)<-met_srftyp[,1]

library(survival)
library(survminer)
# try the Survival thing
load("data/010_METABRIC_survival.rda")
#surftypes<-met_srftyp
surftypes<-setNames(met_srftyp[,2],rownames(met_srftyp))
common<-intersect(names(surftypes),rownames(survival))
subtypes<-surftypes[common]
survival<-survival[common,]
times<-survival[,1]
# times[times>8000]<-8000
# survival[,1]<-times
# survival[survival[,1]>8000]
sfit<-survfit(survival~subtypes)
palette<-c("plum2","tomato","cyan","seagreen3","darkkhaki")
png("plots/007_survival_METABRIC_subtypes_SVM.png",w=3000,h=4000,res=350)
ggsurvplot(sfit, data=survival,palette=palette,
           conf.int=FALSE, pval=TRUE, risk.table=TRUE, 
           title="Kaplan-Meier Curve for METABRIC survival, Surface Genes")
dev.off()
save(subtypes,file="data/000_SVM_subtypes.rda")
