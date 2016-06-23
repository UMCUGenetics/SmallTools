library("VariantAnnotation")
library("RColorBrewer")
library("ggplot2")
library("gplots")
# -------------------------------------------------------------
# SET TO THE RELEVANT DIRECTORY FOR ALL VCFs


# EDIT ME
vcfdir <- "~/data/24SNPs/VCFbased/vcfs"
templatevcf <- "32SNPtaqman_design.vcf"


today <- format(Sys.time(), "%d_%b_%Y")
# -------------------------------------------------------------
# COLORS FOR THE DIFFERENT GENOTYPES
myColors <- c("blue","green","red","white")
names(myColors) <- c("0/0","0/1","1/1",".")

myValues <- c(-1,0,1,NA)
names(myValues) <- c("0/0","0/1","1/1",".")

matrixcols2 <-  colorRampPalette(c("#377eb8","#4daf4a","#e41a1c"))
matrixcols3 <-  colorRampPalette(c("#377eb8","#4daf4a","#e41a1c"))
# -------------------------------------------------------------
# CAN BE SET UP FOR DIFFERENT PROJECTS
project <- "OpenArrayTest"
panel <- "taqman_32SNP"
vcffiles <- list.files(vcfdir, full.names=T, pattern=".vcf$")

# -------------------------------------------------------------
# READ VCF FILES
calls <- lapply(vcffiles, readVcf, "hg19")
# EXTRACT ALL SAMPLE NAMES
samplenames <- unlist(lapply(calls, function(x) samples(header(x)) ))

# EXTRACT SNP data from VCF
template <- readVcf(templatevcf, "hg19")
nrpositions <- nrow(geno(template)$GT)
snpnames <- names(ranges(template))
rm(template)

template <- data.frame(read.table(templatevcf, sep='\t', header=F))
colnames(template) <- c("Chr","Pos","rsID","Ref","Alt","Qual","Info","Format","Test")
template$Merge <- paste0(as.character(template$Chr),":",template$Pos)
rownames(template) <- snpnames


# -------------------------------------------------------------
# PREPARE GENOTYPE DATAFRAME OF CORRECT DIMENSIONS
genot <- data.frame(  matrix(vector(), length(vcffiles), nrpositions, dimnames=list(samplenames, snpnames)) )

# FILL THE GENOTYPE DATAFRAME
for (s in calls) {
  samplename <- samples(header(s))
  nddf <- data.frame(genotype=unlist(geno(s)$GT))
  colnames(nddf) <- samplename

	# ENSURE GENOTYPES ARE MERGED CORRECTLY
  snps <- names(ranges(s))
  if (grepl("_" ,snps[1])) {
    snpdf <- data.frame(Merge=gsub("_.*","",snps))
    snpdf <- merge(template, snpdf, by="Merge")
    rownames(snpdf) <- snpdf$Merge
		# PREVENTS UNORDERED DATA TO MESS THINGS UP
    snps<-as.character(snpdf[gsub("_.*","",snps),]$rsID)
  }
  rownames(nddf) <- snps
	# STORE IN COMPLETE DATAFRAME
  genot[samplename,snps] <- nddf[,1]
}

rm(calls)

# -------------------------------------------------------------

# CONVERT GENOTYPES TO NUMERIC VALUES FOR HEATMAP
num_gentyp <- apply(genot, 2, function(x) myValues[x])
rownames(num_gentyp) <- samplenames
colnames(num_gentyp) <- template[snpnames,]$Merge

# WRITE INTERMEDIATE RESULTS TO FILES FOR TRACKING/EXPORT PURPOSES
write.table(genot, file=paste(project,panel,"Genotypes.txt",sep="_"), sep='\t', quote=F)
write.table(num_gentyp, file=paste(project,panel,"NumGenotypes.txt",sep="_"), sep='\t', quote=F)

# -------------------------------------------------------------
# REMOVE NOISY SAMPLES
# FIXME ADD LOW QC TO OUTPUT
NAcount <- rowSums(is.na(num_gentyp))
remove <- names(NAcount[NAcount>=15])
clean_num_gentyp <- num_gentyp[!rownames(num_gentyp) %in% remove, ]

TAD.dist    <- dist(clean_num_gentyp, method="manhattan", diag=FALSE, upper=FALSE)
TAD.cluster <- hclust(TAD.dist, method="average", members=NULL)
grouping <- cutree(TAD.cluster, h=5)

temper <- data.frame(cbind(clean_num_gentyp, grouping))
temper$Class <- "identified"

pdf(file=paste0(today,"_Genotyping_newpanel_perGroup.pdf"), width=30, height=15, pointsize=15)
for (i in c(1:max(grouping))) {
  rows <- rownames(clean_num_gentyp)[grouping==i]
  plotname <- i
  if ( length(rows) <= 1) {
    temper[rows,]$Class <- "single"
    rows <- rownames(clean_num_gentyp)[grouping%in%c(i-1,i,i+1)]
    plotname <- paste0(i-1,"-",i+1)
  }
  dat <- clean_num_gentyp[rows,]

  heatmap.2(as.matrix(dat), trace="none", scale="none", main=plotname, margins=c(10, 25), col=matrixcols3, cexCol=1.2, cexRow=1.2, dendrogram="row", Colv=F, symkey=T, sepwidth=c(0.02,0.02), sepcolor="lightgray", colsep=1:ncol(dat), rowsep=1:nrow(dat), keysize=0.5)

}
dev.off()


write.table(temper, paste0(today,"_Grouped_Calls_",panel,".txt"), sep='\t', row.names=T, quote=F)
rm(genot)
# -------------------------------------------------------------

# -------------------------------------------------------------
# TODO
# ADD GROUPING CHECK BASED ON PROVIDED PAIRS
# STUB
# groupcols <- unlist(lapply(samplenames, function(x) sampledetails$Col[sampledetails$ID==x]))
# -------------------------------------------------------------

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


  heatmappy <- heatmap.2(data, trace="none", scale="none", margins=c(10, 10), col=matrixcols2, cexCol=1.2, cexRow=.3, dendrogram=clusteringflag, Rowv=rowflag, Colv=colflag, symkey=F,                           sepwidth=c(0.002,0.002), sepcolor="lightgray", colsep=1:ncol(data), rowsep=1:nrow(data), keysize=0.5)
  #if (! is.na(groupcols)) {
  #heatmappy <- heatmap.2(data, trace="none", scale="none", margins=c(10, 10), col=matrixcols2, cexCol=1.2, cexRow=1.4, dendrogram=clusteringflag, Rowv=rowflag, Colv=colflag, symkey=F,  RowSideColors=groupcols, sepwidth=c(0.02,0.02), sepcolor="lightgray", colsep=1:ncol(data), rowsep=1:nrow(data), keysize=0.5)
  #}
  return(heatmappy)
}

pdf(file=paste(project, panel,"genotyping_GATK.pdf", sep="_"), width=15, height=8, pointsize=10)

  dd <- dist(num_gentyp, method="euclidean")
  hc <- hclust(dd, method="ward.D2")
  hcd <- as.dendrogram(hc)

  # Dendrogram
  par(mar=c(2,2,2,8))
  nodePar <- list(lab.cex=0.4, pch=c(NA, NA), cex=0.7, col="blue")
  plot(hcd,  xlab="Height", horiz=TRUE, nodePar=nodePar)

  # Clustering
  make_heatmap(num_gentyp, "row", NA)

  #No clustering
  make_heatmap(num_gentyp, "none", NA)

  # WITH GROUPING
  #make_heatmap(num_gentyp, "row", groupcols)
  #legend("topright", legend=unique(sampledetails$Individual), col=unique(sampledetails$Col), pch=20, ncol=length(unique(sampledetails$Individual)), bty="n", inset=c(0.01,-0.01))

dev.off()

# -------------------------------------------------------------
rm(num_gentyp)
