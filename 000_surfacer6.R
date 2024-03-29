setwd("F:/projects/Breast/")
library(stringr)
library(matrixStats)
library(corto)
library(ggforce)
library(edgeR)
# surfaceR approach on TCGA ovarian cancer
load("F:/Projects/Pan-Cancer/data/surfacer_2020.rda")
load("results/000_BRCA-expmat.rda")
load("results/000_BRCA-rawcounts.rda")
load("results/000_patient_mapping.rda")
# Perform DE Analysis for each PAM50 subtype
# load subtypes
load("data/tcga_BRCA-subtypes.rda")
subtypes[subtypes=="Tumor, Normal-Like"]<-"NormLike"
subtypes<-subtypes[rownames(my_patients)]
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
# normal    Basal     Her2     LumA     LumB NormLike 
#   212      148       53      496      199       13
ncol(rawcounts[,annotation=="Basal"]) # 148
ncol(rawcounts[,annotation=="NormLike"]) # 13
# since only 13 patients are classifed as N-like, they were excluded from further analysis
rawcounts<-rawcounts[,annotation!="NormLike"]
dim(rawcounts)
subannot<-subset(annotation, annotation!="NormLike",drop=TRUE)
subannot <- factor(subannot)
subannot
# Perform classic DE + Network
## go for edgeR
edger<-DGEList(rawcounts,group=subannot,genes=rownames(rawcounts))
# Design Matrix
design <- model.matrix(~0+subannot)
colnames(design) <- levels(subannot)
# Filtering to remove low counts
keep <- filterByExpr(edger, design)
table(keep) #FALSE 2343 TRUE 17899
edger<-edger[keep,,keep.lib.sizes=FALSE]
# Normalization for composition bias
edger<-calcNormFactors(edger)
# Dispersion estimation
edger<-estimateDisp(edger,design,robust=TRUE)

# Quasi-Likelihood dispersion estimate
#fit<-glmQLFit(edger,design,robust=TRUE)
fit<-glmQLFit(edger,design)
colnames(fit)
# "normal" "Basal"  "Her2"   "LumA"   "LumB" 
## create contrasts
contrasts <- makeContrasts("Basal-normal","LumA-normal","LumB-normal","Her2-normal",levels = design)
res<-glmQLFTest(fit,contrast=contrasts)
colnames(res)
topTags(res)
pv<-0.05
is.de <- decideTestsDGE(res,p.value=pv)
summary(is.de)
results<-res$table
results$adjPValue<-p.adjust(results$PValue,method="BH")
save(results,file="000_PAM50subtypes_DE.rda")
write.xlsx(results,file="results/000_Supp_Table1.xlsx",sheetName = "SURFACER",append=T)
############ Extract all surface genes
common<-intersect(rownames(results),surfacer)
results<-results[common,]
## Add the network activity step
# MRA Breast vs. Normal
# Load the network
load("data/GTEx_BRST-regulon.rda")
expmat<-expmat[rowVars(expmat)>0.1,]
tissues<-c("Basal","LumA","LumB","Her2")
for (tissue in tissues){
  mr<-mra(expmat[,annotation==tissue],expmat[,annotation=="normal"],regulon=regulon,minsize=15,nperm=1000,nthreads = 7,verbose=TRUE)
  save(mr,file=paste0("results/001_tcga_",tissue,"-mra.rda"))
}
surfres<-results[,1:4]
allsiggenes<-c()
for (tissue in tissues){
  load(paste0("results/001_tcga_",tissue,"-mra.rda"))
  masters<-names(mr$nes[abs(mr$nes)>2])
  allsiggenes<-unique(c(allsiggenes,masters))
}
dim(surfres) # 2368
keep<-intersect(rownames(surfres),allsiggenes)
surfres<-surfres[keep,]
dim(surfres) #1723

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
#   colnames(el_res)<-c("Basal","LumA","LumB","Her2+")
#   el_res[el_res>5]<-5
#   el_res[el_res< -5]<-(-5)
#   png(paste0("plots/000_heatmap_",class,".png"),w=3500,h=7500,res=500)
#   #heatmap.3(el_res)
#   pheatmap(el_res,
#            color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
#            breaks = breaksList,clustering_distance_cols="correlation") # Sets the breaks of the color scale as in breaksList)
#   
#   dev.off()
# }
dim(surfgenes)
dim(surfres) #1738    4
# now go for the Venn part
library(VennDiagram)
library(ggplot2)
colnames(surfres)<-c("Basal" , "LumA" , "LumB","Her2+")
a<-rownames(surfres[abs(surfres$Basal)>1,]) #Basal
b<-rownames(surfres[abs(surfres$`Her2+`)>1,]) #Her 
c<-rownames(surfres[abs(surfres$LumA)>1,]) # lumA
d<-rownames(surfres[abs(surfres$LumB)>1,]) #lumB
#e<-rownames(surfres[abs(surfres$`N-like`)>1,])# normlike 
gp<-venn.diagram(x=list(a,b,c,d),
                 category.names = c("Basal" , "Her2+" , "LumA","LumB"),filename="plots/000_Venn_PAM50.png",
                 output = TRUE ,
                 imagetype="png" ,
                 height = 1000 , 
                 width = 1500 , 
                 resolution = 450,
                 compression = "lzw",
                 lwd = 1,
                 col=c("red", 'purple', 'darkblue','cyan'),
                 fill = c(alpha('red',0.3), alpha('purple',0.3), alpha('darkblue',0.3),
                          alpha('cyan',0.3)),
                 cex = 0.5,
                 fontfamily = "sans",
                 cat.cex = 0.3,
                 cat.default.pos = "outer",
                 cat.fontfamily = "sans",
                 cat.col = c("red", 'purple', 'darkblue','cyan'))
print(gp)
allsubs<-list(a,b,c,d)
names(allsubs)<-c(letters[1:4])
ol = calculate.overlap(x = allsubs)
ol_size=sapply(ol, length)
ol_size #see numbers of interest
overlap <- calculate.overlap(allsubs)
basalgenes<-overlap$a9
her2genes<-overlap$a14
lumagenes<-overlap$a1
lumbgenes<-overlap$a3
#nlikegenes<-overlap$a5
# THE SURFACER GENE SIGNATURE
genesignature<-c(basalgenes,her2genes,lumagenes,lumbgenes)
save(genesignature,file = "results/000_PAM50signature.rda")
#####################################

source("F:/Archive/heatmaps.R")
# Global
dim(surfres) #1723
plotres<-surfres
plotres<-plotres[genesignature,]
plotres[plotres>2.5]<-2.5
plotres[plotres< -2.5]<-(-2.5)
# load surface markers classifier
genelist<-read.delim("F:/Datasets/Datasets/Liste/surfacer.txt",as.is = TRUE)
rownames(genelist)<-genelist$Gene_Name
genes<-intersect(rownames(plotres),rownames(genelist))
genelist<-genelist[genes,]
table(genelist$Functional_Main_Class)
genelist$Functional_Main_Class[genelist$Functional_Main_Class=="Enzymes"]<-"gray14"
genelist$Functional_Main_Class[genelist$Functional_Main_Class=="Receptors"]<-"red4"
genelist$Functional_Main_Class[genelist$Functional_Main_Class=="Transporters"]<-"lightseagreen"
genelist$Functional_Main_Class[genelist$Functional_Main_Class=="Structural/Adhesion_Molecules"]<-"yellow2"
rowside<-setNames(genelist$Functional_Main_Class,rownames(genelist))
rowside<-as.data.frame(rowside)
rowside$cluster<-c(rep("red",length(basalgenes)),rep("purple",length(her2genes)),rep("darkblue",length(lumagenes)),
                   rep("cyan",length(lumbgenes)))
colnames(rowside)[1]<-"class"
png("plots/000_mapofgois_PAM50.png",w=2000,h=6000,res=350)
heatmap.3(plotres,KeyValueName = "Log2FC",RowSideColors = rowside,Rowv = F)
dev.off()

# add GENE-centered survival info with cox-hazard
load("data/tcga_BRCA-clinical.rda")
load("data/tcga_BRCA-survival.rda")
load("results/000_BRCA-expmat.rda")
expmat<-expmat[,colnames(rawcounts)]
patients<-intersect(colnames(rawcounts),rownames(survival))
survival<-survival[patients,]
expmat<-expmat[genesignature,patients]
surv_target <- data.frame(Symbol = character(),
                          Cox_coef = numeric(),
                          Cox_p = numeric(),
                          Survival_p = numeric(),
                          Pred_sign = character(),
                          stringsAsFactors = F)

pb <- txtProgressBar(min = 0, max = nrow(expmat), initial = 0, style = 3)
if(!file.exists("results/Survival_analysis_target.rda")){
  
  for (i in genesignature) {
    
    library(survival)
    
    mygene <- i
    
    #### Fit Cox model
    res.cox <- summary(coxph(survival ~ expmat[mygene,]))
    
    #### Survfit
    oritrack <- expmat[as.character(mygene),]
    
    #### Define types for survfit
    track<-oritrack
    track[]<-"Other"
    track[oritrack>median(oritrack)]<-paste0(mygene,"up")
    track[oritrack<=median(oritrack)]<-paste0(mygene,"dn")
    track<-as.factor(track)
    track<-relevel(track,ref=paste0(mygene,"dn"))
    
    #### Survdiff for all vs. all p-value
    sdiff<-survdiff(survival~track)
    p <- 1-pchisq(sdiff$chisq, df=length(sdiff$n) - 1)
    
    
    #### Fede function to determine predictor sign
    events_difference<-sdiff$obs-sdiff$exp
    if(events_difference[1]>events_difference[length(events_difference)]){
      sign<-"neg"
    } else {
      sign<-"pos"
    }
    
    #### Group results
    surv_target[i,"Symbol"] <- mygene
    surv_target[i,"Cox_coef"] <- res.cox$coefficients[1,1]
    surv_target[i,"Cox_p"] <- res.cox$coefficients[1,5] 
    surv_target[i,"Survival_p"] <- p
    surv_target[i,"Pred_sign"] <- sign
    
    setTxtProgressBar(pb, value = i)
  }
  save(surv_target, file = "results/Survival_analysis_target.rda")
  
}else{load("results/Survival_analysis_target.rda")}

# colside prognosis and hazard-cox ratio
# univariate Cox's regression model coefficients (pink represents positives coefficients
# (bad prognosis), while light blue are negatives coefficients (good prognosis)) 
# and its corresponding FDR values (black boxes represent
# FDR value for Cox's coefficients < 0.05)
rowside$HR<-ifelse(surv_target$Pred_sign=="neg","dodgerblue4","deeppink2")
rowside$FDR<-ifelse(surv_target$Cox_p<=0.05,"black","white")
rowside<-rowside[,c(4,3,1,2)]

###### newmap
png("plots/000_all_PAM50.png",w=4500,h=7500,res=600)
heatmap.3(plotres,KeyValueName = "Log2FC",RowSideColors = rowside,Rowv = F,Colv = F)
dev.off()




# divide for any cluster
a<-plotres[basalgenes,]
# side<-setNames(rowside$class,rownames(rowside))
# side<-as.data.frame(side)
# side<-side[basalgenes]
side<-rowside[basalgenes,]
png("plots/000_basal_PAM50.png",w=5500,h=7000,res=600)
heatmap.3(a,KeyValueName = "Log2FC",RowSideColors = side,Rowv = T,Colv = F)
dev.off()
b<-plotres[lumagenes,]
side<-rowside[lumagenes,]
png("plots/000_lumA_PAM50.png",w=5500,h=3500,res=600)
heatmap.3(b,KeyValueName = "Log2FC",RowSideColors = side,Rowv = T,Colv = F)
dev.off()
c<-plotres[lumbgenes,]
side<-rowside[lumbgenes,]
png("plots/000_lumB_PAM50.png",w=5500,h=3500,res=600)
heatmap.3(c,KeyValueName = "Log2FC",RowSideColors = side,Rowv = T,Colv = F)
dev.off()
d<-plotres[her2genes,]
side<-rowside[her2genes,]
png("plots/000_her2_PAM50.png",w=5500,h=4500,res=600)
heatmap.3(d,KeyValueName = "Log2FC",RowSideColors = side,Rowv = T,Colv = F)
dev.off()
# e<-plotres[nlikegenes,]
# side<-rowside[nlikegenes,]
# png("plots/000_nlike_PAM50.png",w=3500,h=3500,res=600)
# heatmap.3(e,KeyValueName = "Log2FC",RowSideColors = side,Rowv = T,Colv = F)
# dev.off()

####### ontology (does not have a great meaning at all...)
library(enrichR)
dbs<-listEnrichrDbs()
dbs<-c("WikiPathways_2019_Human","KEGG_2019_Human","MSigDB_Hallmark_2020","GO_Biological_Process_2018")
# Cluster 1
go<-enrichr(basalgenes,dbs)
go$WikiPathways_2019_Human
go<-enrichr(her2genes,dbs)
go$WikiPathways_2019_Human
