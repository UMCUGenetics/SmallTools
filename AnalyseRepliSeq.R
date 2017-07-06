library(ggplot2)
library(RColorBrewer)
library(vioplot)
# ----------------------------------------------------------------------

sample_order <- c("STE007623G1b","STE007623S1","STE007623S2","STE007623S3","STE007623S4","STE007623G2")
myColors <- brewer.pal(6,"RdYlBu")
names(myColors) <- sample_order


stats <- data.frame(read.table("persamplestats.txt", header=T, sep='\t'))
summarystats <- data.frame(sums=colSums(stats))
summarystats$corfactor <- 4000000/summarystats$sums

# ----------------------------------------------------------------------

# Load Freely Avaliable Median replication values (mostely cell lines)
med_ratios <-  data.frame(read.table("all_RepliSeq_median.bed", header=F, sep='\t'))
colnames(med_ratios) <- c("Chr","Start","Stop","Ratio")
# med_ratios$Ratio <- med_ratios$Ratio/16

# ----------------------------------------------------------------------

# Load readcount data
counts <- data.frame(read.table("STE007623_50kb_slide1kb_repliseq_depth.txt", header=F, sep='\t'))
colnames(counts) <- c("Chr","Start","Stop","Reads","meanCoverage","Sample")
counts <- counts[,c("Chr","Start","Stop","Reads","Sample")]

# ----------------------------------------------------------------------
# SET AND DEDUCE ORDERING
counts$Sample <- factor(counts$Sample, levels=sample_order)
data_order <- as.character(counts$Sample[1:6])


# ----------------------------------------------------------------------
# Normalize to 4 million reads
counts$Normalized <- counts$Reads * summarystats[data_order, 2]

# Normalise for Sample ReadDepth
normtot <- tapply(counts$Normalized, counts$Sample, FUN=median)
summarystats$normtot <- as.numeric(unname(normtot))
counts$NormFraction <- counts$Normalized / summarystats[data_order, 3]

# Plot chr21
ggplot(subset(counts, Chr==21), aes(Stop, NormFraction)) + geom_line(aes(colour=Sample), alpha=1) + facet_wrap(~Chr, scales="free_x") + scale_colour_manual(values=myColors) + ylim(c(0,20)) + geom_point(data=subset(counts[!saveregions,],Chr==21), aes(x=Stop),y=10)


#ggplot(subset(counts, Chr==21), aes(x=Sample, y=NormFraction, colour=Sample)) + geom_boxplot()
#ggplot(subset(counts, Chr==21), aes(x=Sample, y=Normalized, colour=Sample)) + geom_boxplot()
#ggplot(subset(counts, Chr==21), aes(x=Sample, y=Reads, colour=Sample)) + geom_boxplot()

#ggplot(subset(acounts, Chr==21), aes(x=Sample, y=NormFraction, colour=Sample)) + geom_boxplot()
#ggplot(subset(acounts, Chr==21), aes(x=Sample, y=Normalized, colour=Sample)) + geom_boxplot()
#ggplot(subset(acounts, Chr==21), aes(x=Sample, y=Reads, colour=Sample)) + geom_boxplot()


# ----------------------------------------------------------------------

# Calculate ratios
ratios <- subset(counts, Sample=="STE007623G1b")[,1:3]
 
early <- subset(counts, Sample=="STE007623G1b")$NormFraction + subset(counts, Sample=="STE007623S1")$NormFraction
late <- subset(counts, Sample=="STE007623S4")$NormFraction + subset(counts, Sample=="STE007623G2")$NormFraction
ratios$Ratio <- early/late
rm(early)
rm(late)

a <- ggplot(subset(ratios, Chr%in%c(1,2,3)), aes(Stop, Ratio)) + geom_line() + facet_wrap(~Chr, scales="free_x", ncol=1)  
a + geom_line(data=subset(med_ratios, Chr%in%c(1,2,3)), col="red")


ggplot(subset(ratios, Chr%in%c(1,2,3)), aes(Stop, Ratio)) + geom_point() + facet_wrap(~Chr, scales="free_x", ncol=1) + ylim(c(0,80)) + geom_line(data=subset(med_ratios, Chr%in%c(1,2,3)), col="red")
#ggplot(subset(med_ratios, Chr%in%c(1,2,3)), aes(Stop, Ratio)) + geom_line(col="red") + facet_wrap(~Chr, scales="free_x", ncol=1)

# ----------------------------------------------------------------------

write.table(ratios, file="STE007623_ReplicationTiming.bed", sep='\t', quote=F)

# ----------------------------------------------------------------------
# Make plots for specific regions/genes

cnvs <- data.frame(read.table("CNVs_Intestine.txt", header=T, sep='\t'))
genes <- data.frame(read.table("FragileGenes.txt", header=T, sep='\t'))
wind_size = 4000000


pdf(file="FragileSites_RepliTiming_V2.pdf", width=15, height=5, pointsize=10)
for (i in c(1:nrow(genes))) {
  wind_start <- genes$Start[i] - wind_size
  wind_stop  <- genes$Stop[i]  + wind_size
  
  fhit_counts <- subset(counts,     Chr==as.character(genes$Chr[i]) & Start>=wind_start & Stop<=wind_stop)
  fhit_ratios <- subset(ratios,     Chr==as.character(genes$Chr[i]) & Start>=wind_start & Stop<=wind_stop) 
  fhit_medrat <- subset(med_ratios, Chr==as.character(genes$Chr[i]) & Start>=wind_start & Stop<=wind_stop)
  
  a <- ggplot(fhit_ratios, aes(Stop, Ratio)) + geom_line(alpha=0.2) + stat_smooth(method="loess", span=0.1) + ylim(c(0,8)) + geom_line(data=fhit_medrat, col="red") + geom_rect(xmin=genes$Start[i], xmax=genes$Stop[i], ymin=-0.1, ymax=0.1) + ggtitle(genes$Name[i])
  
  b <- ggplot(fhit_counts, aes(Stop, NormFraction, colour=Sample)) + scale_colour_manual(values=myColors) + stat_smooth(aes(fill=Sample), method="loess", span=0.2) + scale_fill_manual(values=myColors) + theme(panel.background = element_rect(fill = "darkgrey")) + scale_y_continuous(limits = c(0, NA)) + geom_rect(xmin=genes$Start[i], xmax=genes$Stop[i] , ymin=0, ymax=0.5e-7)  + ggtitle(genes$Name[i])
  
  # ADD CNVS
  fhit_cnvs   <- subset(cnvs, Chr==as.character(genes$Chr[i]))
  if (nrow(fhit_cnvs) >= 1) {
    fhit_cnvs$Ratio=0.1
    fhit_cnvs$NormFraction=0.0
    
    a <- a + geom_rect(data=fhit_cnvs, fill='red', color='white', aes(xmin=Start, xmax=Stop, ymin=Ratio-0.1,    ymax=Ratio))
    b <- b + geom_rect(data=fhit_cnvs, fill='red', color='white', aes(xmin=Start, xmax=Stop, ymin=NormFraction, ymax=NormFraction+0.5e-7))
  }
  
  print(a)
  print(b)
}
dev.off()


# ----------------------------------------------------------------------
# ONCE ONLY 
# Find Noise Regions
# Select regions with value above 15
# Use index list to clean counts and redo analysis

#banlist <- subset(ratios, Ratio>=15)[,1:3]
#bannames <- do.call("paste", banlist)

# Select regions with high coverage on at least 3 timepoints
high_cov_regions <- do.call("paste", subset(counts, NormFraction>=10)[,1:3])
region_count <- table(high_cov_regions)
highcovnames <- names(region_count[region_count>=3])


low_cov_regions <- do.call("paste", subset(counts, NormFraction<=0.10)[,1:3])
region_count <- table(low_cov_regions)
lowcovnames <- names(region_count[region_count>=3])


bannames <- unique(c(lowcovnames,highcovnames))

allnames <- do.call("paste", counts[,1:3])

saveregions <- !allnames %in% bannames


rm(banlist, bannames, allnames, high_cov_regions, region_count, highcovnames)
# ----------------------------------------------------------------------

################################################################################################
################################              CLEAN             ################################
################################################################################################


# ----------------------------------------------------------------------
# Remove Noise regions
acounts <- counts[saveregions,]

# ----------------------------------------------------------------------
# SET AND DEDUCE ORDERING
acounts$Sample <- factor(acounts$Sample, levels=sample_order)
data_order <- as.character(acounts$Sample[1:6])

# ----------------------------------------------------------------------
# Normalize to 4 million reads
totcounts <- aggregate(acounts$Reads, by=list(acounts$Sample), FUN=sum)
colnames(totcounts) <- c("Sample","Reads")
rownames(totcounts) <- totcounts$Sample
totcounts$Corfactor <- 4000000/totcounts$Reads

acounts$Normalized <- acounts$Reads * totcounts[data_order, 3]

# Normalise for Sample ReadDepth
cleannormtot <- tapply(acounts$Normalized, acounts$Sample, FUN=median)
totcounts$cleannormtot <- as.numeric(unname(cleannormtot))
acounts$NormFraction <- acounts$Normalized / totcounts[data_order, 4]

# Plot chr21
ggplot(subset(acounts, Chr==21 & Start>=20000000), aes(Stop, NormFraction)) + geom_line(aes(colour=Sample), alpha=1) + facet_wrap(~Chr, scales="free_x") + scale_colour_manual(values=myColors)

# ----------------------------------------------------------------------

# Calculate cleanratios
cleanratios <- subset(acounts, Sample=="STE007623G1b")[,1:3]

early <- subset(acounts, Sample=="STE007623G1b")$NormFraction + subset(acounts, Sample=="STE007623S1")$NormFraction
late <- subset(acounts, Sample=="STE007623S4")$NormFraction + subset(acounts, Sample=="STE007623G2")$NormFraction
cleanratios$Ratio <- early/late
rm(early)
rm(late)

# med_ratios$Ratio <- med_ratios$Ratio/16

cleanratios$OldRatio <- cleanratios$Ratio
cleanratios$Ratio <- cleanratios$OldRatio*13
cleanratios$Stop <- cleanratios$Start+25e3
#cleanratios$Ratio[cleanratios$Ratio>=100] <- 100

pp <- ggplot(subset(cleanratios, Chr==21), aes(Stop, Ratio, group=Chr)) + stat_smooth(se=F, method="loess",n=70000, span=0.01) + facet_wrap(~Chr, scales="free_x", ncol=1) + ylim(c(0,100)) 
pp + geom_line(data=subset(med_ratios, Chr==21), color="red")

#pp + geom_line(data=subset(ratios, Chr%in%c("1","2","3")), color="black")
#geom_line(alpha=0.3) 
#ggplot(subset(med_ratios, Chr%in%c(1,2,3)), aes(Stop, Ratio)) + geom_line(col="red") + facet_wrap(~Chr, scales="free_x", ncol=1)

# ----------------------------------------------------------------------
write.table(cleanratios, file="STE007623_ReplicationTiming_CleanRegions.bed", sep='\t', quote=F)


library("pastecs")
library("gtools")


region_class <- data.frame()
pdf(file="RepliTiming_Peaks_Consensus.pdf", width=15, height=5, pointsize=10)

chromosomes <- mixedsort(as.character(levels((cleanratios$Chr))))
for (chr in chromosomes) {
  this_chr <- subset(cleanratios, Chr==chr)
  that_chr <- subset(med_ratios, Chr==chr)
  print(chr)
  
  #pred <- loess(Ratio ~ Stop, span=0.01, this_chr)
  #fit <- data.frame(pred$x)
  #fit$Ratio <- pred$fitted
  fit <- that_chr
  
  ts_y<-ts(fit$Ratio)
  tp<-turnpoints(ts_y)
  peaks <- data.frame(Stop=fit$Stop[tp$tppos], Ratio=fit$Ratio[tp$tppos])
  
  peaks$slope <- c(diff(peaks$Ratio)/diff(peaks$Stop),0)
  peaks$sig <- "none"
  peaks$sig[peaks$slope>=0.00003] <- "valley"
  peaks$sig[peaks$slope<=-0.00003] <- "peak"
  peaks$sig <- factor(peaks$sig, levels=c("peak","none","valley"))
  
  
  pq <- ggplot(this_chr, aes(Stop, Ratio)) + geom_line(alpha=0.2) +  ylim(c(-10,110)) + geom_line(data=that_chr, color="red") + geom_line(data=fit, color="black") + geom_point(data=peaks, aes(col=sig), pch=20, size=5) + ggtitle(chr)
  

  stapos <- this_chr$Stop[1]
  endpos <- NA
  class <- "undefined"
  
  nonecounter <- 0
  nonestart <- NA
  
  for (i in c(1:nrow(peaks))) {
    
      # Valley -> Peak    or     Peak -> Valley
      if ((class=="left" && peaks$sig[i]=="peak") || (class=="right" && peaks$sig[i]=="valley"))  {
        endpos <- peaks$Stop[i]
        new_df <- data.frame(Chr=chr, Start=stapos, Stop=endpos, Class=class, Ratio=0)
        region_class <- rbind(region_class, new_df)
        
        stapos <- peaks$Stop[i]
        if (class=="left") { class <- "right" } else { class <- "left"}
        nonecounter <- 0
      }
      
      # Peak -> Peak     or     Valley -> Valley
      else if ((class=="left" && peaks$sig[i]=="valley") || (class=="right" && peaks$sig[i]=="peak"))  {
        # Expand region
        endpos <- peaks$Stop[i]
      }
      
      # None -> None
      else if ((class=="undefined" && peaks$sig[i]=="none"))  {
        # Expand region
        endpos <- peaks$Stop[i]
        nonecounter <- 0
      }
    
      # None -> Peak/Valley
      else if ((class=="undefined" && peaks$sig[i]!="none"))  {
        endpos <- peaks$Stop[i]
        new_df <- data.frame(Chr=chr, Start=stapos, Stop=endpos, Class=class, Ratio=0)
        region_class <- rbind(region_class, new_df)
        
        stapos <- peaks$Stop[i]
        if (peaks$sig[i]=="peak")   { class <- "right" } 
        if (peaks$sig[i]=="valley") { class <- "left"}
      }
     
      # Peak/Valley -> None
      else {
        # First none in region
        if (nonecounter==0) {
          nonstart <- peaks$Stop[i]
          
          # Expand region
          endpos <- peaks$Stop[i]
        }
        
        # Add none counter
        nonecounter <- nonecounter + 1
        
        # If more than 1 none
        if (nonecounter >= 2) {
          new_df <- data.frame(Chr=chr, Start=stapos, Stop=nonstart, Class=class, Ratio=0)
          region_class <- rbind(region_class, new_df)
          
          class <- "undefined"
          stapos <- nonstart
          nonecounter <- 0
        }
      }
      
  }
  
  # add last region us undefined
  new_df <- data.frame(Chr=chr, Start=region_class$Stop[nrow(region_class)], Stop=this_chr$Stop[nrow(this_chr)], Class="undefined", Ratio=0)
  
  region_class <- rbind(region_class, new_df)
  region_class$Class <- factor(region_class$Class, levels=c("left","undefined","right"))
  
  # Change telomeric/centromeric regions
  region_class$Class[(region_class$Stop-region_class$Start)>=10000000] <- "undefined"
  
  
  this_region <- subset(region_class, Chr==chr)
  
  print(pq + geom_rect(data=this_region, aes(fill=Class, xmin=Start, xmax=Stop, ymin=Ratio-1, ymax=Ratio+1), alpha=.5, col='black'))
  #+ geom_vline(xintercept=peaks$Stop[1:16])
  
}


dev.off()

write.table(region_class, file="ReplicationTiming_ClassifiedRegions_Consensus.bed", sep='\t', quote=F, row.names=F)


