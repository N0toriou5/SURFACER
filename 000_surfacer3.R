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
cluster1<-names(subtypes[subtypes=="Basal-enriched"])
mat<-spmra[,cluster1]
signature_basal_like<-apply(mat,1,median)
cluster2<-names(subtypes[subtypes=="Lum1"])
mat<-spmra[,cluster2]
signature_lum1<-apply(mat,1,median)
cluster3<-names(subtypes[subtypes=="Lum2"])
mat<-spmra[,cluster3]
signature_lum2<-apply(mat,1,median)
cluster4<-names(subtypes[subtypes=="Lum3"])
mat<-spmra[,cluster4]
signature_lum3<-apply(mat,1,median)
cluster5<-names(subtypes[subtypes=="Mixed"])
mat<-spmra[,cluster5]
signature_mixed<-apply(mat,1,median)

# create a MRA matrix
surfacer_breast<-cbind(signature_basal_like,signature_lum1,signature_lum2,signature_lum3,signature_mixed)
range(surfacer_breast)
range(spmra)

# assign subtypes at metabric on TCGA models
load("results/000_metabric_spmra.rda")
centroids<-as.data.frame(surfacer_breast)
assignPAM50<-function(expmat=expmat,centroids=centroids){
  genes<-intersect(rownames(centroids),rownames(expmat))
  submat<-expmat[genes,]
  centromat<-centroids[genes,]
  cormat<-cor(submat,centromat,method="pearson")
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
png("plots/010_survival_METABRIC_subtypes.png",w=3000,h=4000,res=350)
ggsurvplot(sfit, data=survival,palette=palette,
           conf.int=FALSE, pval=TRUE, risk.table=TRUE, 
           title="Kaplan-Meier Curve for METABRIC survival, Surface Genes")
dev.off()

### explore if diffentially expressed genes are the same
