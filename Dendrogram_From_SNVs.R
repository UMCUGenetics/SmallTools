#-------------------------------------------------------------------------------------------------------------------------#
#
# Rscript Dendrogram_From_SNVs.R [options]
# for help use:
#           Rscript Dendrogram_From_SNVs.R --help
#
#-------------------------------------------------------------------------------------------------------------------------#
#
# Details for parameters are found here:
#           http://cran.r-project.org/web/packages/SNPRelate/SNPRelate.pdf
#
#-------------------------------------------------------------------------------------------------------------------------#
require(optparse)
require(SNPRelate)

require(ape)
require(sparcl)
require(ggdendro)
source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
#-------------------------------------------------------------------------------------------------------------------------#
options <- list(
    make_option(c("-v", "--verbose"), action="store_true",  default=TRUE,     help="Print extra output [default]"),
    make_option(c("-q", "--quietly"), action="store_false", dest="verbose",   help="Print little output"),

    make_option(c("-i", "--input"),   type="character",                       help="Input file with variants", metavar="vcf"),
    make_option(c("-o", "--output"),  type="character",                       help="Output file for plotting", metavar="pdf"),
    make_option(c("-t", "--title"),   type="character", default="Dendrogram", help="Title for the plot [default %default]", metavar="title"),

    # Dissimilarity matrix computation options (snpgdsDiss)
    make_option("--num.thread",       type="integer", default=4,    help="Number of threads to use [default %default]", metavar="number"),
    make_option("--autosome.only",    type="logical", default=TRUE, help="if TRUE, use autosomal SNVs only [default %default]", metavar="logical"),
    make_option("--remove.monosnp",   type="logical", default=TRUE, help="if TRUE, remove monomorphic SNVs [default %default]", metavar="logical"),
    make_option("--maf",              type="double",  default=0.00, help="SNVs with minor allele frequency >= x are used [default %default]", metavar="number"),
    make_option("--missing.rate",     type="double",  default=0.00, help="SNVs with missing rate <= x are used [default %default]", metavar="number"),

    # Clustering tree cutting options (snpgdsCutTree)
    make_option("--permut",           type="integer", default=5000, help="Number of permutations [default \"%default\"]", metavar="number"),
    make_option("--outlier.n",        type="integer", default=1,    help="Clusters with size <= x are outliers [default \"%default\"]", metavar="number"),
    make_option("--z.threshold",      type="integer", default=6,    help="Z-score to determine when to split nodes [default \"%default\"]", metavar="number")

)
#-------------------------------------------------------------------------------------------------------------------------#
cat("RUNNING\n\n")
cat("-------------------------------------\n")

parser <- OptionParser(usage = "%prog [options]", option_list=options)
arguments <- parse_args(parser, args=commandArgs(trailingOnly=TRUE),  positional_arguments=FALSE)


if(is.null(arguments$input) | is.null(arguments$output)) {
    if (is.null(arguments$input)) {
        stop("Input file not specified\nUse --help for options")
    } else {
        stop("Output file not specified\nUse --help for options")
    }
} else {
    if( file.access(arguments$input) == -1) {
        stop(sprintf("Specified file ( %s ) does not exist", arguments$input))
    }
}

#-------------------------------------------------------------------------------------------------------------------------#

vcffile  <- arguments$input
plotfile <- arguments$output
gdsfile  <- sub(".vcf", ".gds", vcffile)

set.seed(100)

if ( ! file.exists(gdsfile)) {
    # Generate GDS file
    if ( arguments$verbose ){ write("Generating GDS file", stderr()) }
    snpgdsVCF2GDS(vcffile, gdsfile)
}


# Load GDS file
genotypes <- openfn.gds(gdsfile)

# Compute dissimilarity matrix
if ( arguments$verbose ){ write("Performing individual dissimilarity analysis", stderr()) }
diss <- snpgdsDiss(genotypes, maf=arguments$maf, verbose=arguments$verbose, num.thread=arguments$num.thread, autosome.only=arguments$autosome.only, missing.rate=arguments$missing.rate, remove.monosnp=arguments$remove.monosnp)

# Cluster samples
if ( arguments$verbose ){ write("Clustering samples", stderr()) }
clustering <- snpgdsHCluster(diss)

# Determine clusters of samples
if ( arguments$verbose ){ write("Determining branches", stderr()) }
grouping <- snpgdsCutTree(clustering, label.H=FALSE, label.Z=FALSE, verbose=arguments$verbose, n.perm=arguments$permut, outlier.n=arguments$outlier.n, z.threshold=arguments$z.threshold, col.list=NULL, col.outlier="red", pch.outlier='X')

mar.default <- c(0,0,0,0) + 0.1
# Make plots
if ( arguments$verbose ){ write("Generating plots", stderr()) }
pdf(file=plotfile, width=15, height=15, pointsize=12, bg="white")
  
    par(mar = mar.default + c(2, 2, 2, 2))
    # vector of colors
    groupcolors = rainbow(length(grouping$clust.count))

    # Plot the dendrogram
    snpgdsDrawTree(grouping, main=arguments$title, shadow.col=NULL, outlier.col=NULL, leaflab="perpendicular", labels=clustering$sample.id, y.label=-0.2)
#edgePar=list(col=rgb(0.5,0.5,0.5, 0.75), t.col="black"),

    # Plot the distribution of Z scores
    snpgdsDrawTree(grouping, type="z-score", main=arguments$title)

    
    ggdendrogram(clustering$hclust)
    ggdendrogram(clustering$hclust, rotate=TRUE, size=4, theme_dendro=TRUE, color="red")
    
    op=par(bg="darkgrey")
    A2Rplot(clustering$hclust, k=length(grouping$clust.count), boxes=FALSE, col.up="black", col.down=groupcolors)

    op=par(bg="white")
    plot(clustering$hclust)
    plot(clustering$hclust, hang=-1)

    plot(clustering$dendrogram)
    plot(clustering$dendrogram, type="triangle")
    
    plot(as.phylo(clustering$hclust))
    plot(as.phylo(clustering$hclust), type="cladogram")
    plot(as.phylo(clustering$hclust), type="unrooted")
    plot(as.phylo(clustering$hclust), type="fan")
    plot(as.phylo(clustering$hclust), type="radial")

    plot(as.phylo(clustering$hclust), tip.color=groupcolors[grouping$samp.group], col="red")
    plot(as.phylo(clustering$hclust), type="cladogram", tip.color=groupcolors[grouping$samp.group], col="red")
    plot(as.phylo(clustering$hclust), type="unrooted", tip.color=groupcolors[grouping$samp.group], col="red")
    plot(as.phylo(clustering$hclust), type="fan", tip.color=groupcolors[grouping$samp.group], col="red")
    plot(as.phylo(clustering$hclust), type="radial", tip.color=groupcolors[grouping$samp.group], col="red")

dev.off()

cat("-------------------------------------\n")

cat("DONE\n\n")
#-------------------------------------------------------------------------------------------------------------------------#

