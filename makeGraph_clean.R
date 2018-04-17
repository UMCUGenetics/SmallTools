

maxLevelToPlot <- 3
ploidy <- 2

rat_files <- list.files("/home/hub_cuppen/jdeligt/projects/CNVs_freec/new/", full.names=T, pattern="*_ratio.txt")


make_plot <- function(filename) {
	print(filename)

	ratio<-data.frame(read.table(filename, sep='\t', header=T))
	ratio <- subset(ratio, Ratio!=-1)
	ratio$MedianRatio[ratio$Ratio > maxLevelToPlot] <- maxLevelToPlot

	png(filename=gsub(".txt",".png",filename), width=1180, height=1180, units="px", pointsize=20, bg="white", res=NA)
	plot(1:10)
	op <- par(mfrow = c(5,5))

	for (i in c(1:22)) {
		tt <- which(ratio$Chromosome==i)

		if (length(tt)>0) {
			minx <- max(0,min(ratio$Start[tt], na.rm=T))
			maxx <- max(ratio$Start[tt], na.rm=T)
			#print(c(i, minx, maxx))

			# BACKBONE
			plot(ratio$Start[tt], ratio$MedianRatio[tt]*ploidy, xlim=c(minx, maxx), ylim=c(0,maxLevelToPlot*ploidy), xlab=paste0("position, chr",i), ylab="normalized copy number profile", pch=".", col='darkgray')
			
			# GAINS
			tt <- which(ratio$Chromosome==i  & ratio$CopyNumber>ploidy )
			points(ratio$Start[tt],ratio$MedianRatio[tt]*ploidy,pch = ".",col='blue')
			
			tt <- which(ratio$Chromosome==i  & ratio$MedianRatio==maxLevelToPlot & ratio$CopyNumber>ploidy)	
			points(ratio$Start[tt],ratio$MedianRatio[tt]*ploidy,pch=".",col='darkblue',cex=4)
			 
			# LOSSES
			tt <- which(ratio$Chromosome==i  & ratio$CopyNumber<ploidy & ratio$CopyNumber!=-1)
			points(ratio$Start[tt],ratio$MedianRatio[tt]*ploidy,pch=".",col='red')
			tt <- which(ratio$Chromosome==i)
			 
			#UNCOMMENT HERE TO SEE THE PREDICTED COPY NUMBER LEVEL:
			#points(ratio$Start[tt],ratio$CopyNumber[tt], pch = ".", col = colors()[24],cex=4)
		}	
	}

	dev.off()
}

for (file in rat_files) {
	make_plot(file)	
}

