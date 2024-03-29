# # nearest centroid classifier
# library(lolR)
# data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
# X <- data$X; Y <- data$Y
# model <- lol.classify.nearestCentroid(X, Y)
# 
# ### repeat the example
# # assign metabric patients to SURFACER phenotypes
# setwd("F:/projects/Breast/")
# load("results/000_patient_mapping.rda")
# #load("results/000_BRCA-expmat.rda")
# load("results/000_spmra.rda")
# subtypes<-setNames(my_patients$`Surface Genes`,rownames(my_patients))
# subs<-intersect(colnames(spmra),names(subtypes))
# spmra<-spmra[,subs]
# subtypes[subtypes=="tomato"]<-1#"Lum1"
# subtypes[subtypes=="darkkhaki"]<-5#"Mixed"
# subtypes[subtypes=="seagreen3"]<-3#"Lum3"
# subtypes[subtypes=="cyan"]<-2#"Lum2"
# subtypes[subtypes=="plum2"]<-4#"Basal-enriched"
# subtypes<-as.numeric(subtypes)
# table(subtypes)
# datamatrix<-t(spmra)
# model <- lol.classify.nearestCentroid(datamatrix, subtypes)
# yh<-predict(model,datamatrix) # if you substitute the model tih metabric you should assign it to a cluster
# table(yh)
# 
# load("results/000_metabric_spmra.rda")
# datamatrix<-t(spmra)
# yh<-predict(model,t(spmra))
# table(yh)
# library(survival)
# library(survminer)
# # try the Survival thing
# load("data/010_METABRIC_survival.rda")
# names(yh)<-rownames(datamatrix)
# common<-intersect(names(yh),rownames(survival))
# subtypes<-yh[common]
# survival<-survival[common,]
# times<-survival[,1]
# # times[times>8000]<-8000
# # survival[,1]<-times
# # survival[survival[,1]>8000]
# sfit<-survfit(survival~subtypes)
# palette<-c("plum2","tomato","cyan","seagreen3","darkkhaki")
# png("plots/010_survival_METABRIC_subtypes.png",w=3000,h=4000,res=350)
# ggsurvplot(sfit, data=survival,palette=palette,
#            conf.int=FALSE, pval=TRUE, risk.table=TRUE, 
#            title="Kaplan-Meier Curve for METABRIC survival, Surface Genes")
# dev.off()
# 
# ### try with pamr library
# library(pamr)
# x<-spmra
# mydata<-list(x=spmra,y=factor(subtypes), geneid=as.character(1:nrow(x)))
# mytrain<- pamr.train(mydata)
# saved<-mytrain$centroids
# colnames(saved)<-colnames(centroids)
# centroids<-saved

# assign metabric patients to SURFACER phenotypes
setwd("F:/projects/Breast/")
load("results/000_patient_mapping.rda")
load("results/000_BRCA-expmat.rda")
load("results/000_spmra.rda")
# create subtypes
subtypes<-setNames(my_patients$`Surface Genes`,rownames(my_patients))
subs<-intersect(colnames(spmra),names(subtypes))
spmra<-spmra[,subs]
subtypes[subtypes=="tomato"]<-"Lum1"
subtypes[subtypes=="darkkhaki"]<-"Mixed"
subtypes[subtypes=="seagreen3"]<-"Lum3"
subtypes[subtypes=="cyan"]<-"Lum2"
subtypes[subtypes=="plum2"]<-"Basal-enriched"

table(subtypes)
# Basal-enriched           Lum1           Lum2           Lum3          Mixed 
#           139            261            160            204            145 
# create activity matrix
library(pamr)
x<-spmra
mydata<-list(x=spmra,y=factor(subtypes), geneid=as.character(1:nrow(x)))
mytrain<- pamr.train(mydata,n.threshold = 100)
centroids<-mytrain$centroids
# create a MRA matrix


# assign subtypes at metabric on TCGA models
load("results/000_metabric_spmra.rda")
assignPAM50<-function(expmat=expmat,centroids=centroids){
  genes<-intersect(rownames(centroids),rownames(expmat))
  submat<-expmat[genes,]
  centromat<-centroids[genes,]
  cormat<-cor(submat,centromat,method="spearman")
  subtypes<-apply(cormat,1,function(x){
    subtype<-names(which.max(x))
    return(subtype)
  })
  return(subtypes)
}
# Pure assignment
surftypes<-assignPAM50(expmat=spmra,centroids=centroids)
# Save as a subtype object
save(surftypes,file="data/metabric_BRCA-surftypes.rda")
library(survival)
library(survminer)
# try the Survival thing
load("data/010_METABRIC_survival.rda")
common<-intersect(names(surftypes),rownames(survival))
subtypes<-surftypes[common]
survival<-survival[common,]
times<-survival[,1]
# times[times>8000]<-8000
# survival[,1]<-times
# survival[survival[,1]>8000]
sfit<-survfit(survival~subtypes)
palette<-c("plum2","tomato","cyan","seagreen3","darkkhaki")
png("plots/007_survival_METABRIC_subtypes.png",w=3000,h=4000,res=350)
ggsurvplot(sfit, data=survival,palette=palette,
           conf.int=FALSE, pval=TRUE, risk.table=TRUE, 
           title="Kaplan-Meier Curve for METABRIC survival, Surface Genes")
dev.off()
