library("VariantAnnotation")
library("RColorBrewer")
library("ggplot2")
library("gplots")

# -------------------------------------------------------------
myColors <- c("blue","green","red","white")
names(myColors) <- c("0/0","0/1","1/1",".")

myValues <- c(-1,0,1,NA)
names(myValues) <- c("0/0","0/1","1/1",".")

matrixcols2 <-  colorRampPalette(c("#377eb8","#4daf4a","#e41a1c"))
# -------------------------------------------------------------

project <- "OpenArrayTest"
panel <- "taqman"
vcffiles <- list.files("vcfs", full.names=T, pattern=".vcf$")

# -------------------------------------------------------------

calls <- lapply(vcffiles, readVcf, "hg19")
samplenames <- unlist(lapply(calls, function(x) samples(header(x)) ))

print(length(vcffiles))
print(length(samplenames))

print(vcffiles)
print(samplenames)

nrpositions <- nrow(geno(calls[[1]])$GT)
snpnames <- names(ranges(calls[[1]]))

genot <- data.frame(  matrix(vector(), length(vcffiles), nrpositions, dimnames=list(samplenames, snpnames))  )

i<-0
for (s in calls) {
  i<-i+1
  genot[i,] <- unlist(geno(s)$GT)
}
num_gentyp <- apply(genot, 2, function(x) myValues[x])
rownames(num_gentyp) <- samplenames
colnames(num_gentyp) <- snpnames


write.table(genot, file=paste(project,panel,"Genotypes.txt",sep="_"), sep='\t', quote=F)
write.table(num_gentyp, file=paste(project,panel,"NumGenotypes.txt",sep="_"), sep='\t', quote=F)

rm(calls)
rm(genot)
# -------------------------------------------------------------

# Determine SAMPLE GROUPING
# TODO
# groupcols <- unlist(lapply(samplenames, function(x) sampledetails$Col[sampledetails$ID==x]))

# -------------------------------------------------------------
# PLOT

make_heatmap <- function(data, clusteringflag) {
  rowflag=FALSE
  colflag=FALSE
  if (clusteringflag=="row") {
    rowflag=TRUE
  }  
  if (clusteringflag=="col") {
    colflag=TRUE
  }


  heatmappy <- heatmap.2(data, trace="none", scale="none", margins=c(10, 10), col=matrixcols2, cexCol=1.2, cexRow=1.4, dendrogram=clusteringflag, Rowv=rowflag, Colv=colflag, symkey=F, sepwidth=c(0.02,0.02), sepcolor="lightgray", colsep=1:ncol(data), rowsep=1:nrow(data), keysize=0.5)
  #heatmap.2(data, trace="none", scale="none", margins=c(10, 10), col=matrixcols2, cexCol=1.2, cexRow=1.4, dendrogram="none", Rowv=F, Colv=F, symkey=F,  RowSideColors=groupcols, sepwidth=c(0.02,0.02), sepcolor="lightgray", colsep=1:ncol(inform_genot), rowsep=1:nrow(inform_genot), keysize=0.5)
  #legend("topright", legend=unique(sampledetails$Individual), col=unique(sampledetails$Col), pch=20, ncol=length(unique(sampledetails$Individual)), bty="n", inset=c(0.01,-0.01))
  return(heatmappy)
}

pdf(file=paste(project, panel,"genotyping_GATK.pdf", sep="_"), width=30, height=15, pointsize=15)
  # No clustering
  make_heatmap(num_gentyp, "none")

  # Clustering
  make_heatmap(num_gentyp, "row")

  # WITH GROUPING
  #heatmap.2(genot, trace="none", scale="none", margins=c(10, 10), col=matrixcols2, cexCol=1.2, cexRow=1.4, dendrogram="row", Colv=F, symkey=F, RowSideColors=groupcols, sepwidth=c(0.02,0.02), sepcolor="lightgray", colsep=1:ncol(inform_genot), rowsep=1:nrow(inform_genot), keysize=0.5)
  #legend("topright", legend=unique(sampledetails$Individual), pch=20, ncol=length(unique(sampledetails$Individual)), bty="n")

dev.off()

# -------------------------------------------------------------
rm(num_gentyp)


