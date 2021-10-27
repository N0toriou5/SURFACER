# METABRIC test set
setwd("F:/projects/Breast/")
library(stringr)
library(matrixStats)
library(corto)
library(ggforce)
library(edgeR)
# surfaceR approach on TCGA ovarian cancer
load("F:/Projects/Pan-Cancer/data/surfacer_2020.rda")
load("results/000_Breast_surfNet.rda")
# Run Breast Surface Markers Classifier

# load METABRIC cohort
load("data/metabric_BRCA-expmat.rda")
dim(expmat) # 18011  1980
expmat<-na.omit(expmat)
dim(expmat) # 18004  1980
save(expmat,file="data/metabric_BRCA-expmat.rda")
load("data/metabric_BRCA-subtypes.rda")
table(subtypes)
dim(expmat) # 18004  1980
tum<-names(which(subtypes!="Normal"))
brca<-expmat[,tum]
dim(brca) # 18004  1949
# clin <- read.delim("F:/Projects/Orazio/TEPA/data/brca_metabric_clinical_data.tsv",sep="\t")
brca<-brca[rowVars(brca)>0.01,]
dim(brca) # 18004   1949
spmra<-mra(brca,regulon=surfNet,nthreads=7,verbose=TRUE) # single-patient MRA
save(spmra,file="results/000_metabric_spmra.rda")
