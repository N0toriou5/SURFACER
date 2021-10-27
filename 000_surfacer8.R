setwd("F:/projects/Breast/")
#rppa<-read.delim("data/BRCA.rppa.txt")
library(xlsx)
rppa<-read.xlsx2("data/BRCA_rppa_curated.xlsx",sheetIndex = 1,row.names=TRUE)
#rownames(rppa)<-rppa[,1]
rppa2 <- sapply(rppa[,1:ncol(rppa)],as.numeric)
rownames(rppa2)<-rownames(rppa)
rppa<-rppa2
load("results/000_spmra.rda")
colnames(rppa)<-gsub("\\.","-",colnames(rppa))
colnames(rppa)<-substr(colnames(rppa),1,16)
common<-intersect(colnames(spmra),colnames(rppa))
genes<-intersect(rownames(spmra),rownames(rppa))
spmra<-spmra[genes,common]
rppa<-rppa[genes,common]

# scatterplot
x<-setNames(spmra["ANXA1",],colnames(spmra))
y<-setNames(rppa["ANXA1",],colnames(rppa))
library(corto)
png("plots/8_activity_vs_RPPA.png",w=1500,h=1500,res=300)
scatter(y,x,xlab="Inferred Activity",ylab="RPPA", main="ANXA1")
dev.off()
for (name in genes){
  x<-setNames(spmra[name,],colnames(spmra))
  y<-setNames(rppa[name,],colnames(rppa))
  png(paste0("plots/8_",name,"_activity_vs_RPPA.png"),w=1500,h=1500,res=300)
  scatter(y,x,xlab="Inferred Activity",ylab="RPPA", main=name)
  dev.off()
}

load("results/000_BRCA-expmat.rda")
common<-intersect(colnames(expmat),colnames(rppa))
genes<-intersect(rownames(expmat),rownames(rppa))
expmat<-expmat[genes,common]
rppa<-rppa[genes,common]
x<-setNames(expmat["ERBB2",],colnames(expmat))
y<-setNames(rppa["ERBB2",],colnames(rppa))
plot(x,y)

# SLC family receptors in BRCA
setwd("F:/projects/Breast/")
library(stringr)
library(matrixStats)
library(corto)
library(ggforce)
library(edgeR)
library(DESeq2)
library(xlsx)
#source("F:/Archive/vst.R")
# surfaceR approach on TCGA ovarian cancer
load("F:/Projects/Pan-Cancer/data/surfacer_2020.rda")
#load("data/000_tcga_BRCA-edgeres.rda")
#load("results/000_tcga_BRCA-mra.rda")
#load("000_subtypes_DE.rda")
load("F:/Projects/Breast/000_PAM50subtypes_DE.rda")
load("F:/Projects/Breast/results/001_tcga_LumA-mra.rda")
# single-tissue surfaceR
common<-intersect(rownames(results),surfacer)
surfmk<-results[common,]
sigup<-surfmk[surfmk$logFC.LumA.normal>1&surfmk$adjPValue<=1e-50,]
sigup<-sigup[order(sigup$logFC.LumA.normal,decreasing = T),]
mraup<-names(which(mr$nes>15&mr$pvalue<=1e-50))
surf_mark<-intersect(rownames(sigup),mraup)
sigdn<-surfmk[surfmk$logFC.LumA.normal< -1&surfmk$adjPValue<=1e-50,]
mradn<-names(which(mr$nes< -15&mr$pvalue<=1e-50))
downmark<-intersect(rownames(sigdn),mradn)

slcs<-grep("SLC",c(surf_mark,downmark),value=TRUE)
slcsup<-grep("SLC",surf_mark,value=TRUE)
slcsdn<-grep("SLC",downmark,value=TRUE)
#### DEGs plots
# Volcano Plots
library(EnhancedVolcano)
library(airway)
library(magrittr)
png("plots/008_slcs.png", w=2500,h=2500, res=300)
EnhancedVolcano(results, subtitle = "",
                lab = rownames(results),
                selectLab = slcs,
                x = 'logFC.LumA.normal',
                y = 'adjPValue',
                xlim = c(-10, 10),
                ylim = c(0,300),
                title = 'LumA vs. CTRL',
                pCutoff = 0.01, #0.05 cutoff
                FCcutoff = 1, # 2-fold change
                labFace = "bold",
                labSize = 4,
                col = c('black', 'ligthpink', 'lightblue', 'salmon'),
                colAlpha = 4/5,
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.3,colConnectors = 'gray18',
                caption = paste0('Upregulated SLCs= ',length(slcsup),"\n","Downregulated SLCs= ",length(slcsdn)))+#coord_flip()+.
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

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