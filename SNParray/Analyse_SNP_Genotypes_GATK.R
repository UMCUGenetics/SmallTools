library("VariantAnnotation")
library("RColorBrewer")
library("ggplot2")
library("gplots")
# -------------------------------------------------------------
vcfdir <- "~/data/24SNPs/VCFbased/vcfs"

# -------------------------------------------------------------
myColors <- c("blue","green","red","white")
names(myColors) <- c("0/0","0/1","1/1",".")

myValues <- c(-1,0,1,NA)
names(myValues) <- c("0/0","0/1","1/1",".")

matrixcols2 <-  colorRampPalette(c("#377eb8","#4daf4a","#e41a1c"))
# -------------------------------------------------------------

project <- "OpenArrayTest"
panel <- "taqman"
vcffiles <- list.files(vcfdir, full.names=T, pattern=".vcf$")

# -------------------------------------------------------------

calls <- lapply(vcffiles, readVcf, "hg19")
samplenames <- unlist(lapply(calls, function(x) samples(header(x)) ))

print(length(vcffiles))
print(length(samplenames))

template <- readVcf("taqman_design.vcf", "hg19")
nrpositions <- nrow(geno(template)$GT)
snpnames <- names(ranges(template))
rm(template)

template <- data.frame(read.table("taqman_design.vcf", sep='\t', header=F))
colnames(template) <- c("Chr","Pos","rsID","Ref","Alt","Qual","Info","Format","Test")
template$Merge <- paste0(as.character(template$Chr),":",template$Pos)
rownames(template) <- snpnames


# -------------------------------------------------------------

genot <- data.frame(  matrix(vector(), length(vcffiles), nrpositions, dimnames=list(samplenames, snpnames)) )


for (s in calls) {
  samplename <- samples(header(s))
  nddf <- data.frame(genotype=unlist(geno(s)$GT))
  colnames(nddf) <- samplename
  
  snps <- names(ranges(s))
  if (grepl("_" ,snps[1])) {
    snpdf <- data.frame(Merge=gsub("_.*","",snps)) 
    snpdf <- merge(template, snpdf, by="Merge")
    rownames(snpdf) <- snpdf$Merge
    snps<-as.character(snpdf[gsub("_.*","",snps),]$rsID)
  } 
  rownames(nddf) <- snps
  genot[samplename,snps] <- nddf[,1]
}
colnames(genot) <- template[colnames(genot),]$Merge

# -------------------------------------------------------------

num_gentyp <- apply(genot, 2, function(x) myValues[x])
rownames(num_gentyp) <- samplenames
colnames(num_gentyp) <- snpnames



write.table(genot, file=paste(project,panel,"Genotypes.txt",sep="_"), sep='\t', quote=F)
write.table(num_gentyp, file=paste(project,panel,"NumGenotypes.txt",sep="_"), sep='\t', quote=F)

rm(calls)
rm(genot)
# -------------------------------------------------------------

# TODO
# Determine SAMPLE GROUPING
# Sample in header 
# Directory for vcfs
# groupcols <- unlist(lapply(samplenames, function(x) sampledetails$Col[sampledetails$ID==x]))

# -------------------------------------------------------------
# PLOT

make_heatmap <- function(data, clusteringflag, groupcols) {
  rowflag=FALSE
  colflag=FALSE
  if (clusteringflag=="row") {
    rowflag=TRUE
  }  
  if (clusteringflag=="col") {
    colflag=TRUE
  }
  if (clusteringflag=="both") {
    rowflag=TRUE
    colflag=TRUE
  }


  heatmappy <- heatmap.2(data, trace="none", scale="none", margins=c(10, 10), col=matrixcols2, cexCol=1.2, cexRow=1.4, dendrogram=clusteringflag, Rowv=rowflag, Colv=colflag, symkey=F,                           sepwidth=c(0.02,0.02), sepcolor="lightgray", colsep=1:ncol(data), rowsep=1:nrow(data), keysize=0.5)
  if (! is.na(groupcols)) {
  heatmappy <- heatmap.2(data, trace="none", scale="none", margins=c(10, 10), col=matrixcols2, cexCol=1.2, cexRow=1.4, dendrogram=clusteringflag, Rowv=rowflag, Colv=colflag, symkey=F,  RowSideColors=groupcols, sepwidth=c(0.02,0.02), sepcolor="lightgray", colsep=1:ncol(data), rowsep=1:nrow(data), keysize=0.5)
  }
  return(heatmappy)
}

pdf(file=paste(project, panel,"genotyping_GATK.pdf", sep="_"), width=30, height=15, pointsize=15)
  # No clustering
  make_heatmap(num_gentyp, "none", NA)

  # Clustering
  make_heatmap(num_gentyp, "row", NA)

  # WITH GROUPING
  #make_heatmap(num_gentyp, "row", groupcols)
  #legend("topright", legend=unique(sampledetails$Individual), col=unique(sampledetails$Col), pch=20, ncol=length(unique(sampledetails$Individual)), bty="n", inset=c(0.01,-0.01))

dev.off()

# -------------------------------------------------------------
rm(num_gentyp)


