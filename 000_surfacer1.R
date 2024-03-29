setwd("F:/projects/Breast/")
library(stringr)
library(matrixStats)
library(corto)
library(ggforce)
library(edgeR)
library(DESeq2)
library(xlsx)
#source("F:/Archive/vst.R")
# surfaceR approach on TCGA BRCA cancer
load("F:/Projects/Pan-Cancer/data/surfacer_2020.rda")

# Load normal breast tissue counts (we used count matrices from https://doi.org/10.1038/sdata.2018.61)
# Normal tissue counts
gtex_counts <- read.delim("F:/Projects/Pan-Cancer/data/breast-rsem-count-gtex.txt.gz",sep="\t",header=TRUE,as.is=TRUE)
graw <- sapply(gtex_counts[,3:ncol(gtex_counts)],as.integer)
samplenames<-gsub("\\.","-",colnames(graw))
colnames(graw) <- samplenames
rownames(graw) <- gtex_counts$Hugo_Symbol

# Tumor tissue counts
tcga_counts <- read.delim("F:/Projects/Pan-Cancer/data/brca-rsem-count-tcga-t.txt.gz",sep="\t",header=TRUE,as.is=TRUE)
traw <- sapply(tcga_counts[,3:ncol(tcga_counts)],as.integer)
samplenames<-gsub("\\.","-",colnames(traw))
samplenames<-substr(samplenames,1,16)
colnames(traw) <- samplenames
rownames(traw) <- tcga_counts$Hugo_Symbol
traw<-traw[,unique(colnames(traw))]

# Intersect common genes
common<-intersect(rownames(graw),rownames(traw))
traw<-traw[common,]
ncol(traw) #968
graw<-graw[common,]
ncol(graw) #212
rawcounts<-cbind(traw,graw)

# What is normal, what is cancer
tcgacodes<-substr(colnames(traw),14,16)
table(tcgacodes) # 01A 954, 01B 14
tumor<-traw[,tcgacodes=="01A"]
ncol(tumor) #954
load("data/tcga_BRCA-subtypes.rda") # load subtypes and keep only tumor samples according to PAM50
common<-intersect(colnames(tumor),names(subtypes))
subtypes<-subtypes[common]
tumor<-tumor[,common]
ncol(tumor) #909
gtexcodes <- substring(colnames(graw),1,1)
table(gtexcodes)
gnorm <- graw[,gtexcodes=="G"]
ncol(gnorm) #212
rawcounts <- cbind(gnorm,tumor)
group<-factor(c(rep("normal",ncol(gnorm)),rep("tumor",ncol(tumor))))
names(group)<-colnames(rawcounts)
group<-relevel(group,ref="normal")

# VST
#expmat<-vst(rawcounts)
## Use DESeq2 function
expmat <- vst(rawcounts, blind = TRUE, nsub = 1000, fitType = "parametric")
save(expmat,file = "results/000_BRCA-expmat.rda")
save(rawcounts,gnorm,tumor,file="results/000_BRCA-rawcounts.rda")

if(!file.exists("data/GTEx_BRST-regulon.rda")){
# create the surface activity network for normal breast tissue
regulon<-corto(expmat[,colnames(gnorm)],centroids=surfacer,nbootstraps = 1000,p=1e-10,nthreads=7,
               verbose = TRUE)
save(regulon,file = "data/GTEx_BRST-regulon.rda")
} else {load("data/GTEx_BRST-regulon.rda")}

# Perform classic DE + Network
## go for edgeR
edger<-DGEList(rawcounts,group=group,genes=rownames(rawcounts))
# Design Matrix
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
# Filtering to remove low counts
keep <- filterByExpr(edger, design)
table(keep) #FALSE 3272 TRUE 16970
edger<-edger[keep,,keep.lib.sizes=FALSE]
# Normalization for composition bias
edger<-calcNormFactors(edger)


# Dispersion estimation
edger<-estimateDisp(edger,design,robust=TRUE)

# Quasi-Likelihood dispersion estimate
fit<-glmQLFit(edger,design,robust=TRUE)

# Differential Expression (positive logFC will be upregulated in tumor)
contrast<-makeContrasts(tumor-normal,levels=design)
res<-glmQLFTest(fit,contrast=contrast)

pv<-0.05
is.de <- decideTestsDGE(res,p.value=pv)
summary(is.de)
results<-res$table
results$adjPValue<-p.adjust(results$PValue,method="BH")
sum(results$adjPValue<=pv&results$logFC> 1) #3921
sum(results$adjPValue<=pv&results$logFC< -1) #2759
save(results, file = "data/000_tcga_BRCA-edgeres.rda")
write.xlsx(results,file="results/000_Supp_Table1.xlsx",sheetName = "BRCA",append=F)

# MRA Breast vs. Normal
newcounts<-expmat[rowVars(expmat)>0.01,]
mr<-mra(newcounts[,colnames(tumor)],newcounts[,colnames(gnorm)],regulon=regulon,minsize=15,nperm=1000,nthreads = 7,verbose=TRUE)
save(mr,file="results/000_tcga_BRCA-mra.rda")

# single-tissue surfaceR
common<-intersect(rownames(results),surfacer)
surfmk<-results[common,]
sigup<-surfmk[surfmk$logFC>1&surfmk$adjPValue<=0.01,]
sigup<-sigup[order(sigup$logFC,decreasing = T),]
mraup<-names(which(mr$nes>2&mr$pvalue<=0.01))
surf_mark<-intersect(rownames(sigup),mraup)
sigdn<-surfmk[surfmk$logFC< -1&surfmk$adjPValue<=0.01,]
mradn<-names(which(mr$nes< -2&mr$pvalue<=0.01))
downmark<-intersect(rownames(sigdn),mradn)

# Plot
source("F:/Archive/textplot3.R")
#marks<-c("MUC1","ERBB2","PGR","ESR1","SLC39A6","CEACAM6","CLCN3","CHRNB2","KCNK1","LRRC8E")
#marks<-sort(marks)
marks<-surf_mark[1:10]
# plot archetype MD plot
results<-results[order(results$logFC,decreasing = T),]
x<-setNames(results$logCPM,rownames(results))
y<-setNames(results$logFC,rownames(results))
up<-rownames(results)[results$adjPValue<=0.05 & results$logFC>1]
dn<-rownames(results)[results$adjPValue<=0.05 & results$logFC< -1]
png("plots/000_BRCA-top_targets.png",w=4000,h=3000,res=600)
plot(x,y,pch=".",xlab="Average log CPM",ylab="logFC in tumor samples",main = "BRCA vs. Normal Tissue",cex=1,col="darkgrey")
grid()
points(x[dn],y[dn],pch=20,col="lightblue")
points(x[up],y[up],pch=20,col="lightpink")
points(x[surf_mark],y[surf_mark],pch=20,col="red")
legend("bottomright",pch=20,legend=c(paste0("Upregulated Genes: ",length(up)),
                                     paste0("Critical Surface Markers: ",length(surf_mark))),col=c("lightpink","red"),pt.cex=2)
#textplot3(x[marks],y[marks],words = marks,cex=1,line.col = "grey")
#gp<-mraplot(mrs,mrs=marks,title="Surface Markers Activity")
#plot(gp)
dev.off()

# Run Breast Surface Markers Classifier
brca<-expmat[,colnames(tumor)]
dim(brca) #20242   909
brca<-brca[rowVars(brca)>0.01,]
dim(brca) # 18873   909
spmra<-mra(brca,regulon=regulon,nthreads=7,verbose=TRUE) # single-patient MRA
save(spmra,file="results/000_spmra.rda")

# Plot SP-MRA
genes<-intersect(rownames(spmra),surfacer)
# Create the heatmap
## Prepare table
toshow<-spmra[genes,]
# Pairwise correlation between samples (columns)
cols.cor <- cor(toshow, use = "pairwise.complete.obs", method = "pearson")
clustering_distance_cols = as.dist(1 - cols.cor)
res.hc <- hclust(clustering_distance_cols, method = "ward.D2" )
library(ggplot2)
library(factoextra)
png("plots/000_S2.png",w=2000,h=1500,res=450)
fviz_nbclust(t(toshow), kmeans, diss=clustering_distance_cols,method = "wss")+
geom_vline(xintercept = 5, linetype = 2)+
  labs(subtitle = "Elbow method")
dev.off()
# fviz_nbclust(t(toshow), kmeans, diss=clustering_distance_cols,method = "silhouette")
# fviz_nbclust(t(toshow), kmeans, diss=clustering_distance_cols,method = "gap_stat")

# Plot the heatmap
colside<-cutree(res.hc,k=5)
# Color function
colfun<-colorRampPalette(c("navy","navy","blue","blue","white","red","red","red3","red3"))
toplot<-toshow
toplot[toplot>10]<-10
toplot[toplot< -10]<-(-10)
source("F:/Archive/heatmaps.R")
# png("plots/000_heatmap.png",w=4000,h=3000,res=300)
# heatmap.3(toplot,KeyValueName = "NES",hcc=res.hc,col=colfun,ColSideColors = colside)
# dev.off()


# now add further information in the given heatmap
my_gene_col <- colside
my_gene_col[my_gene_col==1]<-"tomato"
my_gene_col[my_gene_col==2]<-"darkkhaki"
my_gene_col[my_gene_col==3]<-"seagreen3"
my_gene_col[my_gene_col==4]<-"cyan"
my_gene_col[my_gene_col==5]<-"plum2"

# https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
# Now add apm 50 informations
#load("data/tcga_BRCA-subtypes.rda")
patients<-intersect(colnames(toplot),names(subtypes))
subtypes<-subtypes[patients]
my_gene_col<-my_gene_col[patients]
my_patients<-as.data.frame(my_gene_col)
# PLos comp bio Bilal et al., 2013
subtypes[subtypes=="Basal"]<-"red"
subtypes[subtypes=="Her2"]<-"purple"
subtypes[subtypes=="LumA"]<-"darkblue"
subtypes[subtypes=="LumB"]<-"lightblue"
subtypes[subtypes=="Tumor, Normal-Like"]<-"green"
my_patients$subtypes<-subtypes
colnames(my_patients)<-c("Surface Genes","PAM50")
# load surface markers classifier
genelist<-read.delim("F:/Datasets/Datasets/Liste/surfacer.txt",as.is = TRUE)
rownames(genelist)<-genelist$Gene_Name
genes<-intersect(rownames(toshow),rownames(genelist))
genelist<-genelist[genes,]
table(genelist$Functional_Main_Class)
genelist$Functional_Main_Class[genelist$Functional_Main_Class=="Enzymes"]<-"gray14"
genelist$Functional_Main_Class[genelist$Functional_Main_Class=="Receptors"]<-"red4"
genelist$Functional_Main_Class[genelist$Functional_Main_Class=="Transporters"]<-"lightseagreen"
genelist$Functional_Main_Class[genelist$Functional_Main_Class=="Structural/Adhesion_Molecules"]<-"yellow2"
rowside<-setNames(genelist$Functional_Main_Class,rownames(genelist))
rows_labs<-rownames(toplot)
GoI<-c("ERBB2","FOLR1","CD274","SLC9A1","ABCA12")
match(GoI,rows_labs) #537  602  263 1800    3
rows_labs2<-str_replace(rows_labs,".*[A-Z0-9]", " ")
rows_labs2[c(537,602,263,1800,3)]<-GoI
png("plots/000_heatmap_clustering_nolabs.png",w=4000,h=50000,res=300)
heatmap.3(toplot,KeyValueName = "NES",hcc=res.hc,col=colfun,ColSideColors = my_patients,RowSideColors = rowside,labRow = rows_labs2, labCol = FALSE)
dev.off()

# Survival analysis based on clusters
# Examples on https://bioconnector.github.io/workshops/r-survival.html
library(survival)
library(survminer)
#load("data/tcga_BRCA-clinical.rda")
load("data/tcga_BRCA-survival.rda")
#load("F:/Projects/Survivals/TCGA/BRCA/survival.rda")
patients<-intersect(names(colside),names(survival))
survival<-survival[patients]
colside<-colside[patients]
Subtype<-colside
Subtype[Subtype==1]<-"Lum1"
Subtype[Subtype==2]<-"Mixed"
Subtype[Subtype==3]<-"Lum3"
Subtype[Subtype==4]<-"Lum2"
Subtype[Subtype==5]<-"Basal-enriched"
# times<-survival[,1]
# times[times>4000]<-4000
# survival[,1]<-times
# survival[survival[,1]>4000]
sfit<-survfit(survival~Subtype)
summary(sfit)
#plot(sfit)
palette<-c("plum2","tomato","cyan","seagreen3","darkkhaki")
png("plots/000_survival_surfaceR_subtypes.png",w=3000,h=4000,res=350)
ggsurvplot(sfit, data=survival,palette=palette,
           conf.int=FALSE, pval=TRUE, risk.table=TRUE, 
           title="Kaplan-Meier Curve for Breast Cancer Survival")
dev.off()
save(my_patients,colside,file="results/000_patient_mapping.rda")
# perform the same curve for PAM50 aptients
table(my_patients$PAM50)
mypatients2<-my_patients[my_patients$PAM50!="green",]
colside<-setNames(mypatients2$PAM50,rownames(mypatients2))
patients<-intersect(names(colside),names(survival))
survival<-survival[patients]
colside<-colside[patients]
Subtype<-colside
sfit<-survfit(survival~Subtype)
summary(sfit)
# survival PAM50 subtypes
png("plots/000_survival_PAM50_subtypes.png",w=3000,h=4000,res=400)
ggsurvplot(sfit, data=survival,
           conf.int=FALSE, pval=TRUE, risk.table=TRUE, 
           title="Kaplan-Meier Curve for Breast Cancer Survival", 
           risk.table.height=.15, palette = c("darkblue","cyan","purple","red"))
dev.off()

