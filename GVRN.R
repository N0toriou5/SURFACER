# Genomic Variant Associated Risk Network GVARN
library(GenomicRanges)
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(ChIPpeakAnno)
genomeanno <- toGRanges(TxDb.Hsapiens.UCSC.hg38.knownGene, feature = "gene")

# You need a couple packages, a database (db) one with the information, and an annotation to annotate it 
library(org.Hs.eg.db)
library(annotate)

#Example:
#avs<-GRanges(c("chr12:99455122-99455123","chr1:2523811-2523812")) #create GRanges object of SNPs
#genes<-GVARN(snp) #"FAM71C" "LOC100129534" "PEX10" "MORN1"
# coord is position
# avs is avs name
# fixed arguments are borders
GVARN <- function(avs,lborder=(-2.5e5),rborder=2.5e5){
  out<-matrix(ncol=2)
  for (i in 1:nrow(avs)){
  x <- GRanges(avs[i,1])
  name <- avs[i,2]
  regions <- annotatePeakInBatch(x, AnnotationData = genomeanno,
                                 output = "overlapping",
                                 FeatureLocForDistance = "TSS",select = "all",
                                 ignore.strand=TRUE,
                                 bindingRegion = c(lborder, rborder))
  if(is.null(regions$feature) == TRUE) next
  symbols <- getSYMBOL(regions$feature,data='org.Hs.eg')
  newline <- cbind(symbols,name)
  out <- rbind(out,newline)
  }
  out <- out[-1,]
  rownames(out)<-out[,2]
  return(out)
}						       
#save(ext,file="demo/ext_eample.rda")
