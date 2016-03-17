library(ggplot2)
library(reshape2)
library(RColorBrewer)

scriptversion <- "Mutation_Profile_Reconstructer-v1.0.0"
file <- "96_count_matrix.txt"

# TODO
# - Add argument parsing
# - Make min precentage and iteration limit variable

# -----------------------------------------------------------------
# COLORS
CtoA <- rgb(0,170,217,maxColorValue=255)
CtoG <- rgb(5,7,8,maxColorValue=255)
CtoT <- rgb(233,36,24,maxColorValue=255)

TtoA <- rgb(203,202,203,maxColorValue=255)
TtoC <- rgb(160,209,93,maxColorValue=255)
TtoG <- rgb(238,200,196,maxColorValue=255)

mutcolors <- c(rep(CtoA,16),rep(CtoG,16),rep(CtoT,16), rep(TtoA,16),rep(TtoC,16),rep(TtoG,16))
names(mutcolors) = c(rep("C>A",16),rep("C>G",16),rep("C>T",16),rep("T>A",16),rep("T>C",16),rep("T>G",16))
# -----------------------------------------------------------------
# FUNCTIONS
# -----------------------------------------------------------------
reconstruct_profile <- function(to_match, reconstruc) {
  
  if (length(reconstruc) == 0) {
    return(rep(0,nrow(to_match)))
  }
  reconstructer <- unlist(reconstruc)
  if (length(reconstructer) == 1) {
    return(to_match[,names(reconstructer)] * reconstructer)
  }
  
  return(rowSums(as.matrix(to_match[,names(reconstructer)]) %*% diag(reconstructer)))
}


fit_model <- function(to_match, test_sample, reconstructer) {
  df <- to_match
  df$TEST <- test_sample-reconstruct_profile(to_match, reconstructer)
  # MODEL data on available signatures
  model <- glm(formula=TEST ~ ., data=df[,!colnames(df)%in%names(reconstructer)])
  # STORE coefficient data
  coefdata <- data.frame(coef(summary(model)))
  # REMOVE intercept
  coefdata <- coefdata[2:nrow(coefdata),]
  colnames(coefdata) <- c("Estimate", "StdError", "tValue", "pValue")
  # Select positive contributions
  potentials <- rownames(subset(coefdata, pValue<=0.05 & Estimate>=0.05))
  if (length(potentials) == 0 ){
    # NO more candiadtes
    return(reconstructer)
  }

  # RE MODEL with positive interactors  
  remodel <- glm(formula=TEST ~ ., data=df[,c("TEST",potentials)])
  coefdata <- data.frame(coef(summary(remodel)))
  coefdata <- coefdata[2:nrow(coefdata),]
  colnames(coefdata) <- c("Estimate", "StdError", "tValue", "pValue")
  
  hit <- coefdata[which.min(coefdata$pValue),]
  hitname <- rownames(hit)
  hitcontrib <- hit$Estimate
  #determine_contrinbution(df$TEST, to_match[,hitname])
  
  if (hitcontrib <= .1) {
    # To small (<10%%) contribution
    return(reconstructer)
  }
  
  if ((hitcontrib + sum(unlist(reconstructer)))>=1.0) {
    # Would reach >100% explained
    hitcontrib <- 1 - sum(unlist(reconstructer))
  }
  
  reconstructer[[hitname]] <- hitcontrib
  return(reconstructer)
}


plot_spectrum <- function(to_match, to_test, reconstruc, sname) {
  reconstructer <- sort(unlist(reconstruc), decreasing=T)
  toplot <- data.frame(cbind(to_match[,c(names(reconstructer))], to_test[,sname]))
  colnames(toplot) <- c(names(reconstructer), sname)
  
  toplot <- data.frame(cbind(toplot, reconstruct_profile(to_match, reconstructer)))
  colnames(toplot)[ncol(toplot)] <- "Reconstructed"
  
  #for (i in c(1:length(reconstructer))) {
  #  toplot <- data.frame(cbind(toplot, reconstruct_profile(to_match, reconstructer[1:i])))
  #  colnames(toplot)[ncol(toplot)] <- paste0("Recon.",i)
  #}
  toplot$Position <- c(1:96)
  toplot$Mutation <- names(mutcolors)
  
  mtoplot <- melt(toplot, id.vars=c("Position","Mutation"))
  colnames(mtoplot) <- c("Position","Mutation","Signature","Proportion")
  
  plottitle <- paste0(sname,"  |  ", paste(paste(names(reconstructer), round(reconstructer,2) ,sep=" : "), collapse=' | ')) 
  
  a <- ggplot(mtoplot, aes(Position, Proportion, group=Signature)) + geom_bar(aes(fill=Mutation),stat="identity", position="dodge") 
  a <- a + ggtitle(plottitle)
  a <- a + scale_fill_manual(values=mutcolors[mtoplot$Mutation]) 
  a <- a + facet_wrap(~Signature, ncol=1)
  a <- a + theme(axis.text.x = element_text(color="black", size=5, angle=90, vjust=.5)) + scale_x_discrete("Mutation context", c(1:96), rownames(toplot), limits=c(1:96))
  return(a)
}

# -----------------------------------------------------------------
# READ input data
to_test <- t(data.frame(read.table(file, sep='\t', row.names=1, header=TRUE)))
# ADJUST to 0-1 range
to_test <- data.frame(apply(to_test, 2, function(x)x/sum(x)))

# -----------------------------------------------------------------
# READ known signature data
# http://cancer.sanger.ac.uk/cosmic/signatures
to_match <- data.frame(read.table("COSMIC_mutational_signatures.tsv", sep='\t', row.names=1, header=TRUE))

# -----------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------

pdf(file=paste0("Classifications_",scriptversion,".pdf"), width=15, height=10, pointsize=13)
overview_df <- data.frame(Sample=character(), Signature=character(), Contribution=numeric())

for (i in colnames(to_test)) {
  print(i)
  
  # RECONSTRUCT profile
  iter <- 0
  reconstructer <- list()
  while (sum(unlist(reconstructer)) <= 0.8 & iter <= 3) {
    reconstructer <- fit_model(to_match, to_test[,i], reconstructer)
    iter <- iter+1
  }

  # STORE contributions
  for (j in names(reconstructer)) {
    overview_df <- rbind(overview_df, data.frame(Sample=i, Signature=j, Contribution=reconstructer[[j]]))
  }

  # PLOT used signatures and Reconstructed profile
  print(plot_spectrum(to_match, to_test, reconstructer, i))
}
dev.off()

# PLOT contribution overview
overview_df$Signature <- factor(overview_df$Signature, levels=paste0("Signature.",c(1:30)))

library(gtools)
overview_df$Sample <- factor(overview_df$Sample, levels=mixedsort((unique(as.character(overview_df$Sample)))))


pdf(file=paste0("Contributions_",scriptversion,".pdf"), width=5, height=10, pointsize=10)
  print( ggplot(overview_df, aes(x=Sample, y=Contribution)) + geom_bar(aes(fill=Signature), colour='black', stat="identity") + theme(axis.text.x=element_text(angle=90, hjust=.1)) + coord_flip() )
dev.off()

# -----------------------------------------------------------------



