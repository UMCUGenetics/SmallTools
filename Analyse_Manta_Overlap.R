require(optparse)
require(VariantAnnotation)
require(ggplot2)

library("optparse")
library("VariantAnnotation")
library("ggplot2")

#-------------------------------------------------------------------------------------------------------------------------#
options <- list(
		make_option(c("-v", "--verbose"),	action="store_true",	default=TRUE,		help="Print extra output [default]"),
		make_option(c("-q", "--quietly"),	action="store_false",	dest="verbose",	help="Print little output"),

		make_option(c("-s", "--sample"),	 	type="character",										help="Sample variants",		metavar="vcf"),
		make_option(c("-c", "--control"),		type="character",										help="Control variants",	metavar="vcf"),

		make_option(c("-r", "--reference"),	type="character",	default="hg19"		help="Reference genome build", metavar="ref"),
		make_option("--overlap",						type="double",		default=0.85,			help="Maximum fraction to overlap reference calls [default %default]", metavar="number"),
		make_option("--passonly",						type="logical", 	default=TRUE,			help="if TRUE, ignore non PASS SVs [default %default]", metavar="logical"),
		make_option("--ignoretype",					type="logical", 	default=TRUE,			help="!TODO! if TRUE, ignore SV types [default %default] [not implemented yet]", metavar="logical")
)
#TODO Add SV type aware code so we can use the ignoretype flag
#-------------------------------------------------------------------------------------------------------------------------#
# FUNCTIONS

# Calculate Manta percentage non refrence
calc_pnr <- function(x) {
	if (lengths(x) <= 1) {	return(0.0)}
	y <- unlist(x)
	if (y[1]==0 && y[2]==0){return(0.0)}
	if (y[1]==0) {					return(1.0)	}
	return (y[2]/sum(y))
}

# Calculate Manta depth
calc_dp <- function(x) {
	if (lengths(x) <= 1) {return(0)}
	return(sum(unlist(x)))
}

# Process Manta file into a more usable format for filtering purposes
process_manta_vcf <- function(vcffile, filter) {
	vcfdf <- data.frame(rowRanges(vcffile))
	if (filter) {
		vcfdf <- subset(vcfdf, Filter=="PASS")
	}
	vcfdf$end <- info(vcffile)$END
	translocations <- which(is.na(info(vcffile)$END))
	vcfdf$end[translocations] <- vcfdf$start[translocations]+1
	vcfdf$width <- vcfdf$end-vcfdf$start
	vcfdf$type <- vcfdf$ALT
	vcfdf$type[translocations] <- "<TRA>"
	indels <- which(! vcfdf$type %in% c("<DUP:TANDEM>","<TRA>","<DEL>","<INS>","<INV>"))
	sizes <- (unlist(lapply(vcfdf$ALT, nchar)) - unlist(lapply(vcfdf$REF, nchar)))
	vcfdf$size <- sizes
	vcfdf$type[indels][sizes[indels]<0] <- "<DEL>"
	vcfdf$type[indels][sizes[indels]>0] <- "<INS>"

	vcfgr <- GRanges(seqnames=vcfdf$seqnames, ranges=IRanges(start=vcfdf$start, end=vcfdf$end), strand="+", type=vcfdf$type)
	rm(vcfdf)
	return(vcfgr)
}

#-------------------------------------------------------------------------------------------------------------------------#
parser <- OptionParser(usage = "%prog [options]", option_list=options)
arguments <- parse_args(parser, args=commandArgs(trailingOnly=TRUE),	positional_arguments=FALSE)

samplevcf <- readVcf(arguments$sample, arguments$reference)
sample <- process_manta_vcf(samplevcf, arguments$passonly)

control <- process_manta_vcf(readVcf(arguments$control, arguments$reference), FALSE)
#interesting_events <- which(countOverlaps(sample, control, minoverlap=1)==0)

#-------------------------------------------------------------------------------------------------------------------------#

# DETERMIN REGIONS UNIQUE TO SAMPLE
sampleunique <- setdiff(sample, control)
rm(control)

# DETERMINE IF OVERLAP OF SVs WITH UNIQUE REGIONS
hits <- findOverlaps(sample, sampleunique, minoverlap=1)
hitsdf <- data.frame(hits)
hitsdf$origwidth <- width(sample[queryHits(hits)])
# DETERMINE WHAT OVERLAPS SVs WITH UNIQUE REGIONS
overlaps <- pintersect(sample[queryHits(hits)], sampleunique[subjectHits(hits)])
hitsdf$hitwidth <- width(overlaps)
# DETERMINE PERCENTAGE OVERLAP OF SVs WITH UNIQUE REGIONS
hitsdf$perc <-round(hitsdf$hitwidth/hitsdf$origwidth,3)
uniqueness <- aggregate(perc ~ queryHits, data=hitsdf, FUN=sum)

# SELECT EVENTS THAT MATCH THE CRITERIA
selected_events <- subset(uniqueness, perc>=arguments$overlap)$queryHits
outvcf <- paste0(gsub("vcf", "", arguments$sample),"_FilteredFor_",arguments$control)
writeVcf(samplevcf[selected_events,], outvcf)

#-------------------------------------------------------------------------------------------------------------------------#
# GATHER RELEVANT DATA FOR PLOTTING
toplot <- data.frame(Type=as.character(sample$type))
toplot$PairPNR <- 	apply(data.frame(geno(samplevcf)$PR), 1, function(x) calc_pnr(x))
toplot$SplitPNR <-	apply(data.frame(geno(samplevcf)$SR), 1, function(x) calc_pnr(x))
toplot$PairDP <- 		apply(data.frame(geno(samplevcf)$PR), 1, function(x) calc_dp(x))
toplot$SplitDP <- 	apply(data.frame(geno(samplevcf)$SR), 1, function(x) calc_dp(x))
toplot$Unique <- FALSE
toplot$Unique[selected_events] <- TRUE
toplot$DP <-	apply(toplot[,c("PairDP","SplitDP")],   1, max)
toplot$PNR <- apply(toplot[,c("PairPNR","SplitPNR")], 1, max)
#-------------------------------------------------------------------------------------------------------------------------#
# PLOT THE FILTERING OVERVIEW
plotfile <- paste0("SVfiltering_",gsub("vcf", "pdf", arguments$sample))
pdf(file=plotfile, width=15, height=15, pointsize=12, bg="white")
	print(ggplot(toplot, aes(DP, PNR)) + geom_jitter(aes(colour=Type, alpha=Unique), width=.01, height=.05))
dev.off()
#+ geom_vline(xintercept=10)

#-------------------------------------------------------------------------------------------------------------------------#
