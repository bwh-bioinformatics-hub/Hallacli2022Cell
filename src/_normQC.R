###########################################
# Rscript to generate the RLE and NUSe plot to detect outlier
# Usage: Rscript $PATH/_normQC.R genes.fpkm.allSamples.uniq.xls
# Reference: https://www.genevestigator.com/userdocs/manual/qc.html
# Author: Xianjun Dong
# Version: 0.0
# Date: 2014-Aug-8
###########################################

if(!require(ape)) install.packages('ape'); 
library(ape)
if(!require(tidyverse)) install.packages('tidyverse'); 
library(tidyverse)

FPKMfile="genes.fpkm.cufflinks.allSamples.xls";
outputfile=paste0(sub("(.*)\\..*","\\1",FPKMfile),".QC.pdf")

message("loading data...")

fpkm=read.table(file(FPKMfile), header=TRUE, row.names = 1, stringsAsFactors = F, sep = "\t", check.names = F);  # table with header (1st row) and ID (1st column)
fpkm = fpkm[, grep("_[Rr]ep", names(fpkm))]  # only sample columns 
colnames(fpkm)=gsub("FPKM.","",colnames(fpkm))
colnames(fpkm)=gsub("_0$","",colnames(fpkm))
dim(fpkm); head(fpkm)

# normalized to yeast genes
# readcount = read.table("genes.htseqcount.cufflinks.allSamples.xls", header=TRUE, row.names = 1, stringsAsFactors = F, sep = "\t", check.names = F);
# head(readcount)
# spikein_controls_fpkm = colSums(fpkm[grep("^ENSG", rownames(fpkm), invert = T),])
# spikein_controls_count = colSums(readcount[grep("^ENSG", rownames(readcount), invert = T),])
# plot(spikein_controls_fpkm, spikein_controls_count)

# only human genes
fpkm = fpkm[grep("^ENSG", rownames(fpkm)),]; 
dim(fpkm)
# fpkm_norm2spikein = sweep(fpkm,2,spikein_controls_fpkm,"/")
# dim(fpkm_norm2spikein)

pdf(outputfile)

# ========================================
# RLE plot to detect outliers
# ========================================

message("generating RLE plot...")

# RLE: For each gene and each sample, ratios are calculated between the expression of a gene and the median expression of this gene across all samples of the experiment. For each sample, these relative expression values are displayed as a box plot. Since it is assumed that in most experiments only relatively few genes are differentially expressed, the boxes should be similar in range and be centered close to 0.

#Two effects may characterize arrays with lower quality: 1) the spread is greater than that of other arrays from this experiment, and 2) the box is not centered near 0.

# filter genes with 0 in >90% samples
notAllZero <- (rowMeans(fpkm>0)>0.1)
logfpkm=fpkm[notAllZero,]
logfpkm=log10(logfpkm + 1e-4)  # so row value of 0 will be -2 in the transformed value
rle=logfpkm-apply(logfpkm, 1, median) # change "/" to "-" so that we got log(fold-change) which centered on 0 on the RLE plot.
rle=pivot_longer(cbind(ID=rownames(rle), rle), !ID, names_to = "Sample", values_to ="FPKM")
bymedian <- with(rle, reorder(Sample, FPKM, IQR))  # sort by IQR
op=par(mar=c(7,3,3,1))
boxplot(FPKM ~ bymedian, data=rle, outline=F, las=2, boxwex=1, col='gray', cex.axis=0.3, main="Relative Log Expression", xlab="", ylab="RLE", frame=F)

#The other graphical representation (NUSE) represents normalized standard error (SE) estimates from the PLM fit. The SE estimates are normalized such that for each probe set, the median standard error across all arrays is equal to 1. A box plot of NUSE values is drawn for each array. On the NUSE plot, arrays with lower quality will have boxes that are centered higher and/or have a larger spread than the other good quality arrays from the same experiment. Typically, boxes centered above 1.1 represent arrays that have quality problems which are often detected in several of the other QC measures presented in this chapter.

# ========================================
# clustering plot to detect outliers
# ========================================

message("generating clustering plot...")

# clustering on columns
sampleDists = 1 - cor(fpkm, method='spearman')
hc=hclust(as.dist(sampleDists),method = "complete")

genotype=as.numeric(gsub("SNCA(.*)_.*_.*_.*","\\1",hc$labels))-1
timepoint=gsub(".*_.*_h(.*)_.*","\\1",hc$labels)
PFF=gsub(".*_PFF(.*)_.*_.*","\\1",hc$labels)

PFF.colors=topo.colors(length(unique(PFF)))

tree=as.phylo(hc)

myLabels <- c('node', sort(unique(timepoint)))
myColors <- c("black", rainbow(length(unique(timepoint))))
timepoint.colors <- myColors[match(timepoint[tree$edge[,2]], myLabels, nomatch=1)]


par(mar=c(1,1,1,1))
plot(tree, type = "unrooted", 
     cex=.5, lab4ut='axial',underscore = T, 
     tip.color=PFF.colors, 
     edge.lty=genotype,
     edge.color= timepoint.colors, 
     main="Clustering of samples based on Spearman correlation")
legend("bottomleft", 
       c("-- edge line type --",unique(gsub("(SNCA.*)_.*_.*_.*","\\1",hc$labels)), "-- PFF treatment --",unique(PFF), "-- timepoint --",paste("hour",sort(unique(timepoint)))),
       text.col=c('black', 'black', 'black', 'black',PFF.colors, myColors), 
       lty = c(0, 1, 3, rep(0, 9)),
       bty='n', cex=.5)

# ========================================
# PCA
# ========================================
library(tidyr)

message("## PCA based on the most variable genes")

# calculate the variance for each gene
rv <- apply(logfpkm, 1, function(x) sd(x)^2)  # variance is square of sd
cv <- apply(logfpkm, 1, function(x) sd(x) / mean(x))  # CV = (Standard Deviation) / Mean = (Square root of Variance)  / Mean
#plot(rv, cv, xlab="Variance", ylab="Coefficient Variation")
#plot(rv, rowMeans(logfpkm), xlab="Variance", ylab="Mean")
#plot(cv, rowMeans(logfpkm), xlab="Coefficient Variation", ylab="Mean")
# select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(500, nrow(logfpkm)))]
# select the top 25% genes with highest CV.
#select <- rank(cv) / nrow(logfpkm) > 1 - cutoff
# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(logfpkm[select,]))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

d=rownames_to_column(as.data.frame(pca$x)) %>% separate(rowname, c("CONDITION","ID","timepoint"), sep = "_", remove = F)
p=ggplot(d,aes(x=PC1, y=PC2, color=CONDITION, shape=timepoint, label=rowname)) + geom_point(size=3) + 
  geom_text(aes(color=ID), vjust = 0, nudge_y = 0.5, size=3) +
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
  ggtitle("PCA based on logarith fpkm and top 500 most variable genes") + 
  coord_fixed()
print(p)

# ========================================
# D-statistic plot to detect outliers
# ========================================
message("generating D-statistic plot...")

# D-statistic
par(op)
D=apply(1-sampleDists, 1, median)
hist(D, breaks=100, ylab="Number of samples", xlab="D-statistic", main="Histogram of D-statistic")
legend("topleft", paste(names(sort(D[which(D<0.9)])), round(sort(D[which(D<0.9)]),2)),title="Likely outlier (D<0.6)", bty='n')

# ========================================
# correlation between replicates
# ========================================
message("correlation betwee replicates...")

replicated=unique(gsub("(.*)_Rep.*","\\1", colnames(fpkm)))
message(paste(length(replicated),"samples have found with replicates!"))

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- abs(cor(log10(x+1e-3), log10(y+1e-3),method='pearson')) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt) 
  #text(0.5, 0.5, txt, cex = cex.cor * r) 
  legend("center", txt, cex = cex.cor * r, bty='n', xjust=0, adj =c(0.5,0.5)) 
}

for(i in replicated){
  message(paste("plot for sample", i, "..."))
  ii=grep(paste0(i,"_"), colnames(fpkm))
  n=length(ii)
  if(n>2){
    pairs(fpkm[,ii] + 1e-3, log='xy',
          main=paste("Replicates for",i),  pch='.',upper.panel=panel.cor)
  } else {
    par(pty="s")
    plot(fpkm[,ii] + 1e-3, log='xy', pch=20, cex=0.3,
         main=paste("Replicates for",i),
         xlab=colnames(fpkm)[ii][1],
         ylab=colnames(fpkm)[ii][2])
    legend("topleft", paste("r =", round(cor(log10(fpkm[,ii]+ 1e-3),method='pearson')[2], 3)), bty='n', cex=1)
  }
}

message("correlation betwee all vs. all")
pairs(fpkm[,grep("Rep1", names(fpkm))] + 1e-3, log='xy', main="All vs. all correlation",  pch='.',upper.panel=panel.cor)

dev.off()

message("QC done")
