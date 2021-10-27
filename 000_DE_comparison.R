##### Compare with METABRIC subtypes for DE genes
setwd("F:/projects/Breast/")
load("results/000_surfres.rda")
load("results/000_BRCA-expmat.rda")
normexp<-expmat[,grep("GTEX",colnames(expmat))]
load("data/metabric_BRCA-expmat.rda")
load("data/metabric_BRCA-surftypes.rda")
common<-intersect(names(surftypes),colnames(expmat))
expmat<-expmat[,common]
surftypes<-surftypes[common]
common<-intersect(rownames(normexp),rownames(expmat))
expmat<-expmat[common,]
normexp<-normexp[common,]
colnames(expmat)
names(surftypes)
dim(expmat)
dim(normexp)
surftypes<-as.vector(surftypes)
surftypes[surftypes=="Basal-enriched"]<-"Basal"
tumormat<-cbind(expmat,normexp)
annotation<-factor(c(surftypes,rep("normal",ncol(normexp))))
names(annotation)<-colnames(tumormat)


# Establish array design
trts=factor(annotation)
design <- model.matrix(~0+trts)     
rownames(design) <- colnames(tumormat)
colnames(design) <- c("Basal","Lum1","Lum2","Lum3","Mixed","normal")

# diff analysis block 
library(limma)
fit <- lmFit(tumormat,design)
contrast.matrix <- makeContrasts(Basal-normal,Lum1-normal,Lum2-normal,Lum3-normal,Mixed-normal,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#Basal
png("plots/001_basal_comp.png",w=1500,h=1500, res=300)
Basal <- topTable(fit2, coef=1, adjust="BH", number = Inf)
surf<-intersect(rownames(surfres),rownames(Basal))
Basal<-Basal[surf,]
surfres<-surfres[surf,]
x<-setNames(surfres$logFC.plum2.normal,rownames(surfres))
y<-setNames(Basal$logFC,rownames(Basal))
plot(x,y,pch=20,col="grey",xlab="Log2FC (TCGA)",ylab="Log2FC (METABRIC)",main="Basal-enriched",
     xlim=1.1*c(min(x),max(x)))
pcc<-cor.test(x,y)
mtext(paste0("R=",signif(pcc$estimate,2),"   p=",signif(pcc$p.value,3)))
lm1<-lm(y~x)
abline(lm1,lwd=2)
dev.off()

# Lum1
png("plots/001_Lum1_comp.png",w=1500,h=1500, res=300)
Lum1 <- topTable(fit2, coef=2, adjust="BH", number = Inf)
surf<-intersect(rownames(surfres),rownames(Lum1))
Lum1<-Lum1[surf,]
surfres<-surfres[surf,]
x<-setNames(surfres$logFC.tomato.normal,rownames(surfres))
y<-setNames(Lum1$logFC,rownames(Lum1))
plot(x,y,pch=20,col="grey",xlab="Log2FC (TCGA)",ylab="Log2FC (METABRIC)",main="Lum1",
     xlim=1.1*c(min(x),max(x)))
pcc<-cor.test(x,y,method="pearson")
mtext(paste0("R=",signif(pcc$estimate,2),"   p=",signif(pcc$p.value,3)))
lm1<-lm(y~x)
abline(lm1,lwd=2)
dev.off()

#Lum2
png("plots/001_Lum2_comp.png",w=1500,h=1500, res=300)
Lum2 <- topTable(fit2, coef=3, adjust="BH", number = Inf)
surf<-intersect(rownames(surfres),rownames(Lum2))
Lum2<-Lum2[surf,]
surfres<-surfres[surf,]
x<-setNames(surfres$logFC.cyan.normal,rownames(surfres))
y<-setNames(Lum2$logFC,rownames(Lum2))
plot(x,y,pch=20,col="grey",xlab="Log2FC (TCGA)",ylab="Log2FC (METABRIC)",main="Lum2",
     xlim=1.1*c(min(x),max(x)))
pcc<-cor.test(x,y,method="pearson")
mtext(paste0("R=",signif(pcc$estimate,2),"   p=",signif(pcc$p.value,3)))
lm1<-lm(y~x)
abline(lm1,lwd=2)
dev.off()

#Lum3
png("plots/001_Lum3_comp.png",w=1500,h=1500, res=300)
Lum3 <- topTable(fit2, coef=4, adjust="BH", number = Inf)
surf<-intersect(rownames(surfres),rownames(Lum3))
Lum3<-Lum3[surf,]
surfres<-surfres[surf,]
x<-setNames(surfres$logFC.seagreen3.normal,rownames(surfres))
y<-setNames(Lum3$logFC,rownames(Lum3))
plot(x,y,pch=20,col="grey",xlab="Log2FC (TCGA)",ylab="Log2FC (METABRIC)",main="Lum3",
     xlim=1.1*c(min(x),max(x)))
pcc<-cor.test(x,y,method="pearson")
mtext(paste0("R=",signif(pcc$estimate,2),"   p=",signif(pcc$p.value,3)))
lm1<-lm(y~x)
abline(lm1,lwd=2)
dev.off()

#Mixed
png("plots/001_Mixed_comp.png",w=1500,h=1500, res=300)
Mixed <- topTable(fit2, coef=5, adjust="BH", number = Inf)
surf<-intersect(rownames(surfres),rownames(Mixed))
Mixed<-Mixed[surf,]
surfres<-surfres[surf,]
x<-setNames(surfres$logFC.darkkhaki.normal,rownames(surfres))
y<-setNames(Mixed$logFC,rownames(Mixed))
plot(x,y,pch=20,col="grey",xlab="Log2FC (TCGA)",ylab="Log2FC (METABRIC)",main="Mixed",
     xlim=1.1*c(min(x),max(x)))
pcc<-cor.test(x,y,method="pearson")
mtext(paste0("R=",signif(pcc$estimate,2),"   p=",signif(pcc$p.value,3)))
lm1<-lm(y~x)
abline(lm1,lwd=2)
dev.off()