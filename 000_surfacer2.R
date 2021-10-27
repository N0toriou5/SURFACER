setwd("F:/projects/Breast/")
library(stringr)
library(matrixStats)
library(corto)
library(ggforce)
library(edgeR)
library(xlsx)
# surfaceR approach on TCGA ovarian cancer
load("F:/Projects/Pan-Cancer/data/surfacer_2020.rda")
load("results/000_BRCA-expmat.rda")
load("results/000_BRCA-rawcounts.rda")
load("results/000_patient_mapping.rda")
# Perform DE Analysis for each SURFACER subtype
# create subtypes
subtypes<-setNames(my_patients$`Surface Genes`,rownames(my_patients))
subs<-intersect(colnames(rawcounts),names(subtypes))
subcounts<-rawcounts[,subs]
subtypes<-subtypes[subs]
rawcounts <- cbind(gnorm,subcounts)
annotation<-factor(c(rep("normal",ncol(gnorm)),subtypes))
names(annotation)<-colnames(rawcounts)
# Differential Expression Analysis
annotation<-relevel(annotation,ref="normal")
# check if different samples are correctly assigned
table(annotation)
# normal      cyan darkkhaki     plum2 seagreen3    tomato 
#   212       160       145       139       204       261 
ncol(rawcounts[,annotation=="normal"]) # 212
ncol(rawcounts[,annotation=="plum2"]) # 139
# Perform classic DE + Network
## go for edgeR
edger<-DGEList(rawcounts,group=annotation,genes=rownames(rawcounts))
# Design Matrix
design <- model.matrix(~0+annotation)
colnames(design) <- levels(annotation)
# Filtering to remove low counts
keep <- filterByExpr(edger, design)
table(keep) #FALSE 2933 TRUE 17309
edger<-edger[keep,,keep.lib.sizes=FALSE]
# Normalization for composition bias
edger<-calcNormFactors(edger)
# Dispersion estimation
edger<-estimateDisp(edger,design,robust=TRUE)
# Quasi-Likelihood dispersion estimate
#fit<-glmQLFit(edger,design,robust=TRUE)
fit<-glmQLFit(edger,design,robust=TRUE)
colnames(fit)
## create contrasts
contrasts <- makeContrasts("plum2-normal","cyan-normal","tomato-normal","seagreen3-normal","darkkhaki-normal",levels = design)
res<-glmQLFTest(fit,contrast=contrasts)
colnames(res)
topTags(res)
pv<-0.05
is.de <- decideTestsDGE(res,p.value=pv)
summary(is.de)
results<-res$table
results$adjPValue<-p.adjust(results$PValue,method="BH")
save(results,file="000_subtypes_DE.rda")
write.xlsx(results,file="results/000_Supp_Table1.xlsx",sheetName = "PAM50",append=T)
############ Extract all surface genes
common<-intersect(rownames(results),surfacer)
results<-results[common,]
## Add the network activity step
# MRA Breast vs. Normal
# Load the network
load("data/GTEx_BRST-regulon.rda")
expmat<-expmat[rowVars(expmat)>0.1,]
tissues<-c("plum2","cyan","tomato","seagreen3","darkkhaki")
for (tissue in tissues){
  mr<-mra(expmat[,annotation==tissue],expmat[,annotation=="normal"],regulon=regulon,minsize=15,nperm=1000,nthreads = 7,verbose=TRUE)
  save(mr,file=paste0("results/001_tcga_",tissue,"-mra.rda"))
}
dim(results) #2279 9
# evaluate and plot the results
results<-results[results$adjPValue<0.05,]
dim(results) # 2274 9
surfres<-results[,1:5]
allsiggenes<-c()
for (tissue in tissues){
  load(paste0("results/001_tcga_",tissue,"-mra.rda"))
  masters<-names(mr$nes[abs(mr$nes)>2&mr$pvalue<=0.01])
  allsiggenes<-unique(c(allsiggenes,masters))
}
dim(surfres) # 2274
keep<-intersect(rownames(surfres),allsiggenes)
surfres<-surfres[keep,]
dim(surfres) #1704
save(surfres,file="results/000_surfres.rda")
library(pheatmap)
library(RColorBrewer)
breaksList = seq(-5, 5, by = 0.5)
# load surface markers classifier
genelist<-read.delim("F:/Datasets/Datasets/Liste/surfacer.txt",as.is = TRUE)
surfgenes<-genelist[genelist$Gene_Name%in%common,]
rownames(surfgenes)<-surfgenes$Gene_Name
surfgenes[surfgenes$Functional_Main_Class=="Structural/Adhesion_Molecules",]<-"Structural"
classes<-c("Enzymes","Receptors","Structural","Transporters")
# for (class in classes){
#   elements<-rownames(surfgenes[surfgenes$Functional_Main_Class==class,])
#   el_res<-as.matrix(surfres[elements,])
#   el_res<-na.omit(el_res)
#   # trans_res<-trans_res[trans_res[,9]<0.05,]
#   # #trans_res<-trans_res[abs(rowMeans(trans_res[,1:5]))>1,]
#   # trans_res<-trans_res[rowVars(trans_res[,1:5])>1,]
#   el_res<-el_res[rowVars(el_res)>1,]
#   source("F:/Archive/heatmaps.R")
#   #trans_res<-trans_res[,1:5]
#   colnames(el_res)<-c("Basal-enriched","Lum2","Lum1","Lum3","Mixed")
#   el_res[el_res>5]<-5
#   el_res[el_res< -5]<-(-5)
#   png(paste0("plots/002_heatmap_",class,".png"),w=3500,h=7500,res=500)
#   #heatmap.3(el_res)
#   pheatmap(el_res,
#            color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
#            breaks = breaksList,clustering_distance_cols="correlation") # Sets the breaks of the color scale as in breaksList)
#   
#   dev.off()
# }

# Global
# dim(surfres) #1704
# plotres<-surfres
# colnames(plotres)<-c("Basal-enriched","Lum2","Lum1","Lum3","Mixed")
# plotres[plotres>5]<-5
# plotres[plotres< -5]<-(-5)
# paletteLength <- 50
# myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# myBreaks <- c(seq(min(plotres), 0, length.out=ceiling(paletteLength/2) + 1), 
#               seq(max(plotres)/paletteLength, max(plotres), length.out=floor(paletteLength/2)))
# png("plots/000_surfacer_expression.png",w=1500,h=3000,res=500)
# pheatmap(plotres,cluster_cols = T,cluster_rows = T, color=myColor, breaks=myBreaks,border_color = "grey30",clustering_distance_cols="correlation",
#          show_rownames = FALSE)
# dev.off()


# pheatmap(plotres,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
#          breaks = breaksList,clustering_distance_cols="correlation",show_rownames=F)
# dev.off()
#heatmap.3(plotres)

dim(surfgenes) # 2279
dim(surfres) #1704    5
# now go for the Venn part
library(VennDiagram)
library(ggplot2)
a<-rownames(surfres[abs(surfres$logFC.plum2.normal)>1,]) #Basal-like
b<-rownames(surfres[abs(surfres$logFC.cyan.normal)>1,]) #LumA-2
c<-rownames(surfres[abs(surfres$logFC.tomato.normal)>1,]) # LumA-1
d<-rownames(surfres[abs(surfres$logFC.seagreen3.normal)>1,]) #lumA-3
e<-rownames(surfres[abs(surfres$logFC.darkkhaki.normal)>1,])# Mixed 
gp<-venn.diagram(x=list(a,b,c,d,e),
                 category.names = c("Basal-enriched" , "Lum2" , "Lum1","Lum3","Mixed"),filename="plots/002_Venn_surf.png",
                 output = TRUE ,
                 imagetype="png" ,
                 height = 1000 , 
                 width = 1700 , 
                 resolution = 450,
                 compression = "lzw",
                 lwd = 1,
                 col=c("plum2", 'cyan', 'tomato','seagreen3','darkkhaki'),
                 fill = c(alpha('plum2',0.3), alpha('cyan',0.3), alpha('tomato',0.3),
                          alpha('seagreen3',0.3),alpha('darkkhaki',0.3)),
                 cex = 0.5,
                 fontfamily = "sans",
                 cat.cex = 0.3,
                 cat.default.pos = "outer",
                 cat.fontfamily = "sans",
                 cat.col = c("plum2", 'cyan', 'tomato','seagreen3','darkkhaki'))
print(gp)

allsubs<-list(a,b,c,d,e)
names(allsubs)<-c(letters[1:5])
ol = calculate.overlap(x = allsubs)
ol_size=sapply(ol, length)
ol_size #see numbers of interest
overlap <- calculate.overlap(allsubs)
lum2genes<-overlap$a2
library(pheatmap)
colnames(surfres)<-c("Basal-enriched" , "Lum2" , "Lum1","Lum3","Mixed")
plotres<-surfres[lum2genes,]
png("plots/000_lum2_heatmap.png",w=1500,h=2000,res=450)
pheatmap(plotres,color=myColor)
dev.off()
mixedgenes<-overlap$a5
plotres<-surfres[mixedgenes,]
png("plots/000_mixed_heatmap.png",w=1500,h=2000,res=450)
pheatmap(plotres,color=myColor)
dev.off()
basalenrich<-overlap$a1
plotres<-surfres[basalenrich,]
png("plots/000_basalenrich_heatmap.png",w=1500,h=6750,res=450)
pheatmap(plotres,color=myColor)
dev.off()
lum3genes<-overlap$a4
plotres<-surfres[lum3genes,]
png("plots/000_lum3_heatmap.png",w=1500,h=3000,res=450)
pheatmap(plotres,color=myColor)
dev.off()

lum1genes<-overlap$a3
plotres<-surfres[lum1genes,]
png("plots/000_lum1_heatmap.png",w=1500,h=4000,res=450)
pheatmap(plotres,color=myColor)
dev.off()

# THE SURFACER GENE SIGNATURE
genesignature<-c(lum2genes,mixedgenes,basalenrich,lum3genes,lum1genes)
save(genesignature,file = "results/000_surfsignature.rda")
#####################################

setwd("F:/projects/Breast/")
load("results/000_patient_mapping.rda")
load("results/000_surfsignature.rda")
load("results/000_surfres.rda")
source("F:/Archive/heatmaps.R")
# Global
dim(surfres) #1704
plotres<-surfres
colnames(plotres)<-c("Basal-enriched","Lum2","Lum1","Lum3","Mixed")
plotres<-plotres[genesignature,]
plotres[plotres>2.5]<-2.5
plotres[plotres< -2.5]<-(-2.5)
# load surface markers classifier
genelist<-read.delim("F:/Datasets/Datasets/Liste/surfacer.txt",as.is = TRUE)
rownames(genelist)<-genelist$Gene_Name
genes<-intersect(rownames(plotres),rownames(genelist))
genelist<-genelist[genes,]
table(genelist$Functional_Main_Class)
genelist$Functional_Main_Class[genelist$Functional_Main_Class=="Enzymes"]<-"red"
genelist$Functional_Main_Class[genelist$Functional_Main_Class=="Receptors"]<-"yellow"
genelist$Functional_Main_Class[genelist$Functional_Main_Class=="Transporters"]<-"blue"
genelist$Functional_Main_Class[genelist$Functional_Main_Class=="Structural/Adhesion_Molecules"]<-"darkgrey"
rowside<-setNames(genelist$Functional_Main_Class,rownames(genelist))
rowside<-as.data.frame(rowside)
rowside$cluster<-c(rep("cyan",length(lum2genes)),rep("darkkhaki",length(mixedgenes)),rep("plum2",length(basalenrich)),
                  rep("seagreen3",length(lum3genes)),rep("tomato",length(lum1genes)))
colnames(rowside)[1]<-"class"
png("plots/000_mapofgois.png",w=2000,h=4500,res=350)
heatmap.3(plotres,KeyValueName = "Log2FC",RowSideColors = rowside)
dev.off()
