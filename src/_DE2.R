###########################################
# R script for running differntial expression analysis using DESeq2
# author: Xianjun Dong
# email: xdong@rics.bwh.harvard.edu
# date: 2/27/2018
# version: 1.0
# Usage: Rscript /data/rnaseq/src/_DE.R

# Q1: Is there a difference between certain mRNAs (or classes) in steady state between 2copy vs 4 copy synuclein levels?
# Q1: Is there a difference between certain mRNAs (or classes) in steady state  between 2copy vs 4 copy synuclein levels upon PFF treatment? 
# Q3: Is there a global mRNA decay rate difference between 2copy vs 4 copy? If not, are there specific mRNA classes or mRNA motifs in UTRs that show a difference?
# Q4: Is there any effect of PFF on the mRNA decay rate(locally or globally) and if so, is this effect exacerbated by PFF treatment?
## Logic to compute decay rate for each mRNA:  
# 1- Find normalized counts of each mRNA
# 2- Fit a exponential decay curve
# 3- Find Lambda (Nt = N0 e-Lt)
# 4- Find half life 
# 5- Compare halflife of each mRNA from 2 copy synuclein and 4 copy synuclein. 

###########################################
# install packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, googlesheets4, hexbin, pheatmap, RColorBrewer, hwriter, ggforce)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pacman::p_load(vsn, DESeq2, ReportingTools, BiocParallel, limma, EnhancedVolcano)

input_dir="~/projects/vik2021/results/merged"
output_dir="~/projects/vik2021/results/DE"
output_additonal_columns='MI'

# Create folder if the directory doesn't exist
dir.create(file.path(output_dir,'report/figures'), recursive =T, showWarnings = FALSE)

setwd(output_dir)

if(file.exists(file.path(output_dir,"DESeq2.RData"))) load(file.path(output_dir,"DESeq2.RData")) else {
  
  
  ###########################################
  # step1: load data
  ###########################################
  #annotation
  genes_annotation = read.table("~/bioinformatics/referenceGenome/hg38sacCer3_combined/WholeGenomeFasta/combined.gene.bed+2", header = F, stringsAsFactors = F, col.names = c("chr","start","end","geneID","score","strand","geneSymbol","geneType"));
  
  # raw reads count
  cts=read.delim(file.path(input_dir, "genes.htseqcount.cufflinks.allSamples.xls"), row.names = 1)
  head(cts)
  # remove those non-geneID rows, e.g. __no_feature (pre-mRNA reads) and __ambiguous (see http://htseq.readthedocs.io/en/master/count.html)
  dim(cts); cts=cts[grep("_", rownames(cts), invert = T),]; dim(cts);
  
  # yeast genes
  spikein_control=cts[grep("^ENSG", rownames(cts), invert = T),]; dim(spikein_control);
  
  plot(colSums(spikein_control))
  # human genes
  cts=cts[grep("^ENSG", rownames(cts)),]; dim(cts);
  
  ## Normalized to the yeast genes
  # plot(colSums(cts), colSums(spikein_control), xlab="total reads of human gens", ylab="total reads of yeast gens")
  # cts = sweep(cts,2,colSums(spikein_control),"/")
  
  # covariance table
  gsurl="https://docs.google.com/spreadsheets/d/15TN7O0ad9IrdK_PLYNo3EgXbKlqGnsH0igh3I5ThODU/edit#gid=2020811866"
  gs4_auth(email = "xianjun.dong.harvard@gmail.com")
  covarianceTable = read_sheet(gsurl, sheet = 'RNAseq') %>% 
    mutate(Timepoint=factor(Timepoint, levels = paste0("h",c(0,1,3,6,12))),
           CONDITION = factor(Mutant, levels = c("SNCA2","SNCA4")),
           PFF.status = factor(PFF.status, levels = c("PFFno","PFFyes")),
           Replicate = factor(Replicate, levels = paste0("Rep",1:3))) %>% 
    select(SAMPLE_NAME, CONDITION, PFF.status, Replicate, Timepoint) %>% as.data.frame()
  rownames(covarianceTable) = covarianceTable$SAMPLE_NAME;
  str(covarianceTable); dim(covarianceTable); 
  write.table(covarianceTable,file.path(output_dir,"covariance.txt"), sep="\t", col.names = T, row.names = F, quote = F)
  
  # subset and re-order
  all(rownames(covarianceTable) %in% colnames(cts))
  dim(cts); cts = cts[, rownames(covarianceTable)]; dim(cts)
  all(rownames(covarianceTable) == colnames(cts))
  
  ###########################################
  # step2: load data to DEseq
  ###########################################
  
  # Note: With no arguments to results, the results will be for the last variable in the design formula, and if this is a factor, the comparison will be the last level of this variable over the first level.
  # Ref: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
  
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = covarianceTable,
                                design = ~ CONDITION)
  
  #colData(dds)$Diagnosis <- factor(colData(dds)$condition, levels=c("HC","MS"))
  
  ## pre-filtering
  # hist(rowMeans(counts(dds)))
  #keep <- rowSums(counts(dds)) >= 10
  keep <- rowSums(counts(dds) >= 5) >= 3
  dim(dds); dds <- dds[keep,]; dim(dds);
  head(sort(rowMeans(counts(dds)), decreasing = T))
  
  ###########################################
  # step3: QA of the data [optional]
  ###########################################
  # Note: This part is not necessary for DEseq, but important for data QA
  
  ##--------------------------------------
  ## 3.1: compare different variance stabilization methods
  ##--------------------------------------
  
  ntd <- normTransform(dds) # log2(x+1)
  vsd <- vst(dds, blind=T) # Note: blind to the design, equal to design = ~ 1
  
  # using limma to remove covariates, it returns adjusted values in log2 scale
  library(limma)
  vsd_adjusted_log2 <- removeBatchEffect(assay(vsd), batch=vsd$Timepoint, batch2=vsd$PFF.status)
  
  pdf(file.path(output_dir,"diagnosis.pdf"))
  msd <- meanSdPlot(counts(dds), ranks = FALSE); msd$gg + ggtitle("no transformation")
  msd <- meanSdPlot(assay(ntd), ranks = FALSE); msd$gg + ggtitle("log2(x+1) transform")
  msd <- meanSdPlot(assay(vsd), ranks = FALSE); msd$gg + ggtitle("VST")
  msd <- meanSdPlot(vsd_adjusted_log2, ranks = FALSE); msd$gg + ggtitle("vsd_adjusted_log2")
  dev.off()
  
  ##--------------------------------------
  ## 3.2: save normalized reads count
  ##--------------------------------------
  
  # save workspace into ".RData" file
  save.image(file.path(output_dir,"DESeq2.RData"))
  
  ## save the raw reads count 
  write.table(counts(dds), file.path(output_dir,"htseqcount.raw.allsamples.xls"), sep="\t", quote = F, col.names = NA, row.names = TRUE)
  ## save the variance-stabilized data
  write.table(assay(vsd), file.path(output_dir,"htseqcount.vst.allsamples.xls"), sep="\t", quote = F, col.names = NA, row.names = TRUE)
  ## save the variance-stabilized data with covariates adjusted
  write.table(vsd_adjusted_log2, file.path(output_dir,"htseqcount.vsd_adjusted_log2.allsamples.xls"), sep="\t", quote = F, col.names = NA, row.names = TRUE)
  
  ##--------------------------------------
  ## 3.3: clustering of samples
  ##--------------------------------------
  
  #sampleDists <- as.dist((1 - cor(assay(vsd))))
  sampleDists <- dist(t(assay(vsd)))
  
  pdf(file.path(output_dir,"clustering.tree.pdf"))
  
  ## heatmap
  par(cex=0.5, mar=c(5, 8, 4, 1))
  
  sampleDistMatrix <- as.matrix( sampleDists )
  #rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep = " - " )
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           clustering_method = 'ward.D', fontsize = 8,
           col = colors)
  
  ## tree
  
  plot(hclust(sampleDists,method = "ward.D"), xlab='', main="Cluster Dendrogram")
  
  # ## PCA on raw
  # se=SummarizedExperiment(counts(dds), colData=colData(vsd))
  # pcaData <- plotPCA(DESeqTransform(se), intgroup = c("CONDITION"), returnData = TRUE)
  # percentVar <- round(100 * attr(pcaData, "percentVar"))
  # p=ggplot(pcaData, aes(x = PC1, y = PC2, color = CONDITION, label=name)) +
  #   geom_point(size =3) +
  #   geom_text(nudge_x = 0.5, nudge_y=0.5, size=3) +
  #   xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  #   ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  #   ggtitle("PCA based on raw count") 
  # print(p)
  
  ## PCA
  pcaData <- plotPCA(vsd, intgroup = c("CONDITION", "PFF.status", "Timepoint"), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  library(ggrepel)
  p=ggplot(pcaData, aes(x = PC1, y = PC2, color = CONDITION, shape=PFF.status, alpha=Timepoint, label=name)) +
    geom_point(size =3) +
    geom_text_repel(nudge_x = 0.5, nudge_y=0.5, size=2, alpha=1, max.overlaps =100) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    xlim(c(-50,40)) +
    ggtitle("PCA based on vsd") 
  print(p)
  
  dev.off()
  
  save.image(file.path(output_dir,"DESeq2.RData"))
}

# script to generate html report page
makeNewImages <- function(df,...){
  imagename <- c()
  tablename <- c()
  for (i in 1:nrow(df)){
    ensId <- rownames(df)[i]
    symbol <- df$symbol[i]    
    imagename[i] <- paste('plot', ensId, symbol, 'pdf', sep = ".")
    tablename[i] <- paste('plot', ensId, symbol, 'txt', sep = ".")
    
    d <- data.frame(samples=colnames(assay(ntd)), 
                    expression_ntd=assay(ntd)[ensId,], 
                    expression_raw=assay(dds)[ensId,], 
                    condition=colData(dds)$CONDITION)
    
    N=length(levels(colData(dds)$CONDITION))
    p=ggplot(d, aes(x=condition, y=expression_ntd)) + 
      geom_boxplot(position=position_dodge(.8), width=.5, outlier.shape = NA) +
      geom_jitter(size=1.5, position = position_jitter(width=.15)) +
      theme_bw() +
      xlab("CONDITION") + ylab("log2(counts+1)") + ggtitle(symbol,subtitle =ensId)
    
    if(!file.exists(file.path('report/figures',tablename[i]))) {
      write.table(d,file.path('report/figures',tablename[i]),sep="\t", quote =F, row.names=F, col.names = T)
    }
    if(!file.exists(file.path('report/figures',imagename[i]))) {
      #png(file.path('report/figures', imagename[i]), height = 250, width = 600)
      pdf(file.path('report/figures', imagename[i]), height = 4, width = 3*N)
      print(p)
      dev.off()
    }
  }
  ## Using the following code to show thumb figures. It's slow to display if many
  # df$Boxplot <- hwriteImage(paste0('figures/', imagename), 
  #                           link=paste0('figures/', imagename), 
  #                           table=FALSE, width=100)
  df$Boxplot <- hwrite('boxplot', link = paste0('figures/', imagename), table=F)
  df$Rawdata <- hwrite("data", link = paste0('figures/', tablename), table=F)
  df$symbol <- hwrite(as.character(df$symbol), 
                      link = paste0("http://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",as.character(rownames(df))), 
                      table=F)
  return(df)
}

# Computes the half-life of a gene in a transcription inhibition experiment
HalfLivesBlock <- function(y, x=c(0,1,3,6,12), doubling.time = 0) {
  # Computes the half-life of a gene in a transcription inhibition experiment.
  #
  # Args:
  #   x: Vector of experimental time points.
  #   y: Vector of normalized gene abundances.
  #   doubling.time: The cell cycle length of the system under study, used for
  #     correcting determined half-life. No cell cycle correction is performed
  #     if set to 0. Default is 0.
  #
  # Returns:
  #   The determined half-life in the same units used for the time vector x. NA
  #     if the algorithm fails to converge or the determined half-life is
  #     negative.
  suppressWarnings(
    fit <- nls(
      y ~ C * (exp(k * x)),
      start = c(C = max(y), k = -0.25),
      algorithm = "port",
      lower = c(C = 0, k = -Inf), upper = c(C = Inf, k = 0),
      control = list(warnOnly = TRUE)
    )
  )
  
  # check to see if it converged, return if it didn't
  if( !fit$convInfo$isConv ) {
    return( NA )
  }
  
  k <- coef(fit)['k']
  C <- coef(fit)['C']
  hl <- NA
  
  # calculate half-life performing cell cycle correction if required
  if (doubling.time) {
    hl <- log(2) / -(k + log(2) / doubling.time) 
  } else {
    hl <- log(2) / -(k) 
  }
  
  if (hl < 0) {  # we don't want any negative half-lives
    hl <- NA
  }
  
  return(hl)
}

###########################################
# run DE among CONDITION, Timepoint jointly
# to test any genes that react in a condition-specific manner over time
# see http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments
###########################################

for(pff in levels(covarianceTable$PFF.status)){
  message(paste("# Running DEseq on PFF status =", pff, "..."))
  # subsetting
  # debug: pff="PFFno"
  dds_subset=dds[,dds$PFF.status==pff]; 
  dim(dds_subset);
  
  # change Timepoint to numeric variable
  colData(dds_subset)$Timepoint = as.numeric(sub("h","",colData(dds_subset)$Timepoint))
  
  design(dds_subset) <- as.formula(" ~ CONDITION + Timepoint + CONDITION:Timepoint")
  
  # a likelihood ratio test, where we remove the condition-specific differences over time. 
  # Genes with small p values from this test are those which at one or more time points after time 0 showed a condition-specific effect. 
  # Note therefore that this will not give small p values to genes that moved up or down over time in the same way in both strains.
  dds_subset <- DESeq(dds_subset, test="LRT", reduced = ~ CONDITION + Timepoint,
                      parallel=TRUE, BPPARAM=MulticoreParam(4))
  
  resultsNames(dds_subset)
  
  # This will give p-values for LRT test: genes which at one or more time points after time 0 showed a strain-specific effect. 
  # Note that the log2FC is for the interaction term (e.g. "CONDITIONSNCA4.Timepointh12" by default or by the designed name, e.g. name = "CONDITIONSNCA4.Timepointh3"). 
  res <- results(dds_subset, alpha = 0.05, parallel=TRUE, BPPARAM=MulticoreParam(4))  ## Note that the alpha for indepedent filtering has to be same as the number to threshold res$padj (https://support.bioconductor.org/p/71732/#71771)
  
  # extract a matrix of the log2 fold changes for all comparisons using the coef function.
  betas <- coef(dds_subset)
  colnames(betas)  # same as the resultsNames(dds_subset)
  betas = betas[rownames(res), -c(1,2)]
  ## Note that CONDITIONSNCA4.Timepointh1 means difference between SNCA4 vs SNCA2 at h1, controlling for baseline. 
  ## To get the log2 fold change of h1 vs h0 for the SNCA2, Timepoint_h1_vs_h0
  ## To get the log2 fold change of h1 vs h0 for the SNCA4, Timepoint_h1_vs_h0 +  CONDITIONSNCA4.Timepointh1
  betas = as.data.frame(betas) %>% mutate(.keep = "none",
                                          SNCA2  = Timepoint, 
                                          SNCA4  = Timepoint + CONDITIONSNCA4.Timepoint)
  
  colnames(betas) = paste0("log2FC.",colnames(betas)) 
  
  # combine the overall LRT tesst (for any timepoint) and the log2FC for each timepoint together
  res = cbind(as.data.frame(res)[,c('baseMean', 'pvalue','padj','log2FoldChange')], betas)
  # decimal value of Fold-change
  res$log2FC.SNCA2_minus_SNCA4 <- res$log2FoldChange
  
  summary(res)

  # add annotation
  res$symbol <- genes_annotation$geneSymbol[match(row.names(res), genes_annotation$geneID)]
  res$geneType <- genes_annotation$geneType[match(row.names(res), genes_annotation$geneID)]

  ## save the data
  res <- res[order(res$padj),]
  head(res); dim(res)
  
  write.table(as.data.frame(res), 
              file=gzfile(file.path(output_dir, paste0("DEresult.interactionLRT.", pff ,".linear.all.xls.gz"))), 
              sep="\t", quote =F, na="", row.names=T, col.names = NA)
  
  write.table(subset(res, padj<=0.05), 
              file=file.path(output_dir, paste0("DEresult.interactionLRT.", pff ,".linear.padj05.xls")), 
              sep="\t", quote =F, na="", row.names=T, col.names = NA)
  
  ## scatterplot of log2FC.SNCA2 and log2FC.SNCA4
  # res %>% mutate(geneType=ifelse(geneType %in% c("lncRNA","protein_coding"), geneType, "others")) %>% 
  #   mutate(geneType=factor(geneType, levels = c("protein_coding","lncRNA","others"))) %>% 
  #   ggplot(aes(x=log2FC.SNCA2, y=log2FC.SNCA4)) +
  #   geom_point(aes(size=-log10(pvalue), alpha=-log10(pvalue), col=geneType)) +
  #   geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + geom_abline(intercept = 0, slope = 1, linetype=2) + 
  #   coord_fixed(ratio = 1) +
  #   ggtitle(paste("Timepoints regression slopes in SNCA2 vs. SNCA4 for all genes")) 
  # ggsave(file.path(output_dir, paste0("DEresult.interactionLRT.", pff ,".linear.all.scatterplot.pdf")), width = 8, height = 7)
  
  filter(res, pvalue<=.05) %>% mutate(geneType=ifelse(geneType %in% c("lncRNA","protein_coding"), geneType, "others")) %>% 
    mutate(geneType=factor(geneType, levels = c("protein_coding","lncRNA","others"))) %>% 
    ggplot(aes(x=log2FC.SNCA2, y=log2FC.SNCA4)) +
    geom_point(aes(size=-log10(pvalue), alpha=-log10(pvalue), col=geneType)) +
    geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + geom_abline(intercept = 0, slope = 1, linetype=2) + 
    coord_fixed(ratio = 1) +
    ggtitle(paste("Timepoints regression slopes in SNCA2 vs. SNCA4 for genes with LRT p < 0.05")) 
  ggsave(file.path(output_dir, paste0("DEresult.interactionLRT.", pff ,".linear.pvalue05.scatterplot.pdf")), width = 8, height = 7)
  
  ## count plot for individual significant genes ## TOFINISH: LEFT HERE
  if(sum(res$padj<=0.05)>0){
    cnts <- t(counts(dds_subset, normalized = T, replaced = F)[rownames(res)[which(res$padj<=0.05)], ] + 0.5)
    intgroup = c("Timepoint","CONDITION")
    dd <- data.frame(cnts, colData(dds_subset)[rownames(cnts), intgroup])
    dd <- pivot_longer(dd, cols=starts_with("ENSG"), names_to = "geneID", values_to = "count") %>% 
      mutate(geneSymbol=genes_annotation$geneSymbol[match(geneID, genes_annotation$geneID)])
    ggplot(dd, aes(x = Timepoint, y = count, color = CONDITION, group = CONDITION)) + 
      geom_point() + stat_summary(fun.y=mean, geom="line") +
      labs(title=paste0("Normalized counts for genes with significant condition-specific changes over time (LRT padj < 0.05)"), 
           subtitle=file.path(output_dir, paste0("DEresult.interactionLRT.", pff ,".padj05.pdf")), y="Normalized count (log10 scalue)") +
      scale_y_log10() + facet_wrap_paginate(vars(geneSymbol), scales = "free_y", ncol = 3, nrow = 5, page = 1, strip.position = "top") +
      theme(strip.background = element_blank(), strip.placement = "outside") -> p
    required_n_pages <- n_pages(p)
    pdf(file.path(output_dir, paste0("DEresult.interactionLRT.", pff ,".linear.padj05.plotCounts.pdf")), width=8, height = 10, paper = 'letter')
    for(i in 1:required_n_pages){
      ggplot(dd, aes(x = Timepoint, y = count, color = CONDITION, group = CONDITION)) + 
        geom_point() + stat_summary(fun.y=mean, geom="line") +
        labs(title=paste0("Genes with significant condition-specific changes over time (LRT padj < 0.05)"), caption = paste("page",i,"of",required_n_pages),
             subtitle=file.path(output_dir, paste0("DEresult.interactionLRT.", pff ,".linear.padj05.pdf")), y="Normalized count (log10 scalue)") +
        scale_y_log10() + facet_wrap_paginate(vars(geneSymbol), scales = "free_y", ncol = 3, nrow = 5, page =i, strip.position = "top") +
        theme(strip.background = element_blank(), strip.placement = "outside") -> p
      print(p)
    }
    dev.off()
    
    ## now we can cluster the top significant genes (in any timepoint) by their Log2FC profiles.
    topN=100
    topDE = filter(res[1:topN,], padj<=0.05) %>% as_tibble() %>% column_to_rownames('symbol') %>%  select(contains("log2FC"))
    colnames(topDE) = gsub("log2FC.","",colnames(topDE))
    ## trim max and min
    MAX=3; MIN=-3; topDE[topDE>MAX]=MAX;topDE[topDE<MIN]=MIN
    
    ## to make sure the 0-point on the color scale is white (ref: https://stackoverflow.com/a/31707976)
    paletteLength <- 50
    myColor = colorRampPalette(c("blue", "white", "red"))(paletteLength)
    # use floor and ceiling to deal with even/odd length pallettelengths
    myBreaks = c(seq(MIN, 0, length.out=ceiling(paletteLength/2) + 1), 
                 seq(MAX/paletteLength, MAX, length.out=floor(paletteLength/2)))
    
    
    annotation_row = filter(res[1:topN,], padj<=0.05) %>% as_tibble() %>% select(geneType, symbol)  %>% column_to_rownames('symbol')
    annotation_col = data.frame(rowname=names(topDE), CONDITION=sub("\\..*","",names(topDE))) %>% column_to_rownames()
    ann_colors = list(
      #geneType = c(protein_coding = "darkblue", lincRNA = "orange", Mt_rRNA='pink', antisense='yellow', pseudogene='gray', processed_transcript='lightgray',Mt_tRNA='purple',misc_RNA='red',snoRNA='lightblue', polymorphic_pseudogene="gray", processed_pseudogene="gray", TEC="lightred"),
      CONDITION = c(SNCA2 = "#F8766D", SNCA4 = "#00BFC4")
    )
    par(cex=0.5, mar=c(5, 8, 4, 1))
    pheatmap(topDE,
             fontsize = 8,
             main =paste("Heatmap of log2 fold changes for genes with adjusted p < 0.05\n",
                         file.path(output_dir, paste0("DEresult.interactionLRT.", pff ,".linear.padj05.xls"))),
             width=6, height=0.15*nrow(topDE),
             filename=file.path(output_dir, paste0("DEresult.interactionLRT.", pff ,".linear.topP.heatmap.pdf")),
             #border_color = NA,
             color=myColor, breaks=myBreaks,
             annotation_row = annotation_row,
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             drop_levels = TRUE,
             scale = "none", 
             #clustering_method = 'ward.D', 
             cluster_rows = TRUE,
             #clustering_distance_rows = "correlation",
             cluster_cols = F)  
  }
  
  
}

###########################################
# run DE among PFF, Timepoint jointly
# to test any genes that react in a condition-specific manner over time
# see http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments
###########################################

for(condition in levels(covarianceTable$CONDITION)){
  message(paste("# Running DEseq on CONDITION status =", condition, "..."))
  # subsetting
  # debug: condition="SNCA4"
  dds_subset=dds[,dds$CONDITION==condition]; 
  dim(dds_subset);
  
  design(dds_subset) <- as.formula(" ~ PFF.status + Timepoint + PFF.status:Timepoint")
  
  # a likelihood ratio test, where we remove the condition-specific differences over time. 
  # Genes with small p values from this test are those which at one or more time points after time 0 showed a condition-specific effect. 
  # Note therefore that this will not give small p values to genes that moved up or down over time in the same way in both strains.
  dds_subset <- DESeq(dds_subset, test="LRT", reduced = ~ PFF.status + Timepoint,
                      parallel=TRUE, BPPARAM=MulticoreParam(4))
  
  resultsNames(dds_subset)
  
  # This will give p-values for LRT test: genes which at one or more time points after time 0 showed a strain-specific effect. 
  # Note that the log2FC is for the interaction term (e.g. "CONDITIONSNCA4.Timepointh12" by default or by the designed name, e.g. name = "CONDITIONSNCA4.Timepointh3"). 
  res <- results(dds_subset, alpha = 0.05, parallel=TRUE, BPPARAM=MulticoreParam(4))  ## Note that the alpha for indepedent filtering has to be same as the number to threshold res$padj (https://support.bioconductor.org/p/71732/#71771)
  
  # extract a matrix of the log2 fold changes for all comparisons using the coef function.
  betas <- coef(dds_subset)
  colnames(betas)  # same as the resultsNames(dds_subset)
  betas = betas[rownames(res), -c(1,2)]
  ## Note that CONDITIONSNCA4.Timepointh1 means difference between SNCA4 vs SNCA2 at h1, controlling for baseline. 
  ## To get the log2 fold change of h1 vs h0 for the SNCA2, Timepoint_h1_vs_h0
  ## To get the log2 fold change of h1 vs h0 for the SNCA4, Timepoint_h1_vs_h0 +  CONDITIONSNCA4.Timepointh1
  betas = as.data.frame(betas) %>% mutate(.keep = "none",
                                          PFFno.h1_vs_h0  = Timepoint_h1_vs_h0, 
                                          PFFno.h3_vs_h0  = Timepoint_h3_vs_h0, 
                                          PFFno.h6_vs_h0  = Timepoint_h6_vs_h0, 
                                          PFFno.h12_vs_h0 = Timepoint_h12_vs_h0, 
                                          PFFyes.h1_vs_h0  = Timepoint_h1_vs_h0 + PFF.statusPFFyes.Timepointh1,
                                          PFFyes.h3_vs_h0  = Timepoint_h3_vs_h0 + PFF.statusPFFyes.Timepointh3,
                                          PFFyes.h6_vs_h0  = Timepoint_h6_vs_h0 + PFF.statusPFFyes.Timepointh6,
                                          PFFyes.h12_vs_h0 = Timepoint_h12_vs_h0 + PFF.statusPFFyes.Timepointh12)
  
  colnames(betas) = paste0("log2FC.",colnames(betas)) 
  
  # combine the overall LRT tesst (for any timepoint) and the log2FC for each timepoint together
  res = cbind(as.data.frame(res)[,c('baseMean', 'stat','pvalue','padj')], betas)
  
  # add annotation
  res$symbol <- genes_annotation$geneSymbol[match(row.names(res), genes_annotation$geneID)]
  res$geneType <- genes_annotation$geneType[match(row.names(res), genes_annotation$geneID)]
  
  ## save the data
  res <- res[order(res$padj),]
  head(res); dim(res)
  
  write.table(as.data.frame(res), 
              file=gzfile(file.path(output_dir, paste0("DEresult.interactionLRT.", condition ,".all.xls.gz"))), 
              sep="\t", quote =F, na="", row.names=T, col.names = NA)
  
  write.table(subset(res, padj<=0.05), 
              file=file.path(output_dir, paste0("DEresult.interactionLRT.", condition ,".padj05.xls")), 
              sep="\t", quote =F, na="", row.names=T, col.names = NA)
  
  ## any globe difference in slope between PFFno vs. PFFyes
  ttest_pvalue = rownames_to_column(res) %>% 
    pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
    pivot_wider(names_from = PFF.status, values_from = log2FC) %>% 
    filter(pvalue<=1) %>%
    group_by(timepoint_comparison) %>% 
    summarise(PFFno=list(PFFno), PFFyes=list(PFFyes)) %>% rowwise() %>% 
    mutate(pvalue=wilcox.test(unlist(PFFno), unlist(PFFyes), paired = T, alternative = "less")$p.value) %>% # alternative = 'less'?
    mutate(timepoint_comparison=factor(timepoint_comparison, levels = paste0("h",c(1,3,6,12),"_vs_h0"))) 
  
  nDecay = rownames_to_column(res) %>% 
    pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
    pivot_wider(names_from = PFF.status, values_from = log2FC) %>% 
    filter(pvalue<=1) %>%
    group_by(timepoint_comparison) %>% 
    summarise(PFFno.percent=mean(PFFno<0), PFFyes.percent=mean(PFFyes<0), PFFno.n=sum(PFFno<0), PFFyes.n=sum(PFFyes<0)) %>% 
    pivot_longer(cols = contains('PFF'), names_to = c('PFF.status',".value"), names_pattern = "(.*)\\.(.*)") %>% 
    mutate(timepoint_comparison=factor(timepoint_comparison, levels = paste0("h",c(1,3,6,12),"_vs_h0"))) 
  
  pivot_longer(as_tibble(res), starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
    mutate(timepoint_comparison=factor(timepoint_comparison, levels = paste0("h",c(1,3,6,12),"_vs_h0"))) %>%
    filter(pvalue<=1) %>%
    ggplot(aes(x=log2FC, fill=PFF.status, color=PFF.status)) + 
    #geom_histogram(fill="white", alpha=0.5, position="identity", bins = 100) +
    geom_histogram(aes(y=..density..), colour="black", position="identity", alpha=0.3, bins=100)+
    geom_density(alpha=.2) +
    geom_vline(xintercept =0, color='blue', linetype=2) +
    facet_wrap(vars(timepoint_comparison), ncol = 1, nrow = 4, strip.position = "top") + 
    geom_text(data = ttest_pvalue, mapping = aes(x = Inf, y = Inf, label = paste("Paired one-sided Wilcox test p = ",formatC(pvalue, format = "g", digits = 2)), color = NULL,fill= NULL), hjust   = 1.05, vjust   = 1.5) +
    geom_text(data = nDecay, mapping = aes(x = -Inf, y = Inf, label = paste0("Decreasing genes in ",PFF.status, ": n=",n, " (",round(100*percent, 2), "%)"), color = PFF.status), hjust   = -0.05, vjust   = rep(c(1.5, 3),4)) +
    ggtitle(paste0("DEresult.decayDiff.histogram.grouped.", condition ,".all.pdf")) + 
    ggsave(file.path(output_dir, paste0("DEresult.decayDiff.histogram.grouped.", condition ,".all.pdf")), width = 11, height = 10)
  
  # histogram for negative slope genes @ h1
  res_down_down = filter(res, log2FC.PFFno.h1_vs_h0<0, log2FC.PFFyes.h1_vs_h0<0)
  ttest_pvalue = rownames_to_column(res_down_down) %>% 
    pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
    pivot_wider(names_from = PFF.status, values_from = log2FC) %>% 
    filter(pvalue<=1) %>%
    group_by(timepoint_comparison) %>% 
    summarise(PFFno=list(PFFno), PFFyes=list(PFFyes)) %>% rowwise() %>% 
    mutate(pvalue=wilcox.test(unlist(PFFno), unlist(PFFyes), paired = T, alternative = "less")$p.value) %>% # alternative = 'less'?
    mutate(timepoint_comparison=factor(timepoint_comparison, levels = paste0("h",c(1,3,6,12),"_vs_h0"))) 
  
  nDecay = rownames_to_column(res_down_down) %>% 
    pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
    pivot_wider(names_from = PFF.status, values_from = log2FC) %>% 
    filter(pvalue<=1) %>%
    group_by(timepoint_comparison) %>% 
    summarise(PFFno.percent=mean(PFFno<0), PFFyes.percent=mean(PFFyes<0), PFFno.n=sum(PFFno<0), PFFyes.n=sum(PFFyes<0)) %>% 
    pivot_longer(cols = contains('PFF'), names_to = c('PFF.status',".value"), names_pattern = "(.*)\\.(.*)") %>% 
    mutate(timepoint_comparison=factor(timepoint_comparison, levels = paste0("h",c(1,3,6,12),"_vs_h0"))) 
  
  pivot_longer(as_tibble(res_down_down), starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
    mutate(timepoint_comparison=factor(timepoint_comparison, levels = paste0("h",c(1,3,6,12),"_vs_h0"))) %>%
    filter(pvalue<=1) %>%
    ggplot(aes(x=log2FC, fill=PFF.status, color=PFF.status)) + 
    #geom_histogram(fill="white", alpha=0.5, position="identity", bins = 100) +
    geom_histogram(aes(y=..density..), colour="black", position="identity", alpha=0.3, bins=100)+
    geom_density(alpha=.2) +
    geom_vline(xintercept =0, color='blue', linetype=2) +
    facet_wrap(vars(timepoint_comparison), ncol = 1, nrow = 4, strip.position = "top") + 
    geom_text(data = ttest_pvalue, mapping = aes(x = Inf, y = Inf, label = paste("Paired one-sided Wilcox test p = ",formatC(pvalue, format = "g", digits = 2)), color = NULL,fill= NULL), hjust   = 1.05, vjust   = 1.5) +
    geom_text(data = nDecay, mapping = aes(x = -Inf, y = Inf, label = paste0("Decreasing genes in ",PFF.status, ": n=",n, " (",round(100*percent, 2), "%)"), color = PFF.status), hjust   = -0.05, vjust   = rep(c(1.5, 3),4)) +
    ggtitle(paste0("DEresult.decayDiff.histogram.grouped.", condition ,".res_down_down.pdf")) + 
    ggsave(file.path(output_dir, paste0("DEresult.decayDiff.histogram.grouped.", condition ,".subset_negative_h1.pdf")), width = 11, height = 10)
  
  ## monotonic ones (i.e. keep decreasing)
  # all(x == cummin(x)) # see https://stackoverflow.com/a/13094801
  monotonic_decreasing <- function(x, na.rm = T) all(x == cummin(x))
  df=rownames_to_column(res) %>% rowwise() %>% 
    filter(max(log2FC.PFFno.h1_vs_h0, log2FC.PFFno.h3_vs_h0, log2FC.PFFno.h6_vs_h0, log2FC.PFFno.h12_vs_h0)<0,
           max(log2FC.PFFyes.h1_vs_h0, log2FC.PFFyes.h3_vs_h0, log2FC.PFFyes.h6_vs_h0, log2FC.PFFyes.h12_vs_h0)<0) %>% 
    mutate(monotonic_PFFno = monotonic_decreasing(c(log2FC.PFFno.h1_vs_h0, log2FC.PFFno.h3_vs_h0, log2FC.PFFno.h6_vs_h0, log2FC.PFFno.h12_vs_h0)),
           monotonic_PFFyes = monotonic_decreasing(c(log2FC.PFFyes.h1_vs_h0, log2FC.PFFyes.h3_vs_h0, log2FC.PFFyes.h6_vs_h0, log2FC.PFFyes.h12_vs_h0))) %>% 
    filter(monotonic_PFFno | monotonic_PFFyes) %>% print(n_extra=17) 
  with(df, table(monotonic_PFFno, monotonic_PFFyes))
  write.table(df, file=file.path(output_dir, paste0("DEresult.interactionLRT.", condition ,".monotonic.xls")), 
              sep="\t", quote =F, na="", row.names = F)
  
  # if the monotonic ones have bigger decay rate in PFFno vs. PFFyes
  df %>% filter(monotonic_PFFno, monotonic_PFFyes) %>%  #dim()
    #mutate(PFFno_lt_PFFyes = log2FC.PFFno.h1_vs_h0 < log2FC.PFFyes.h1_vs_h0) %>% pull(PFFno_lt_PFFyes) %>% table()
    #mutate(PFFno_minus_PFFyes = log2FC.PFFno.h1_vs_h0-log2FC.PFFyes.h1_vs_h0) %>% pull(PFFno_minus_PFFyes) %>% t.test()
    mutate(PFFno_minus_PFFyes = log2FC.PFFno.h1_vs_h0-log2FC.PFFyes.h1_vs_h0) %>% 
    ggplot(aes(x=PFFno_minus_PFFyes, fill=PFFno_minus_PFFyes<0)) + 
    geom_histogram(colour="black", position="identity", alpha=0.3, bins=100) + xlab("log2FC.PFFno.h1_vs_h0 - log2FC.PFFyes.h1_vs_h0") + 
    ggtitle(paste("Slope difference between PFFno vs. PFFyes at h1 for the monotonic decreasing genes (n =", nrow(df %>% filter(monotonic_PFFno, monotonic_PFFyes)),")")) +
    ggsave(file.path(output_dir, paste0("DEresult.decayDiff.histogram.grouped.", condition ,".monotonic.h1diff.pdf")), width = 10, height = 4)
  
  # # paired difference in slope
  # rownames_to_column(res) %>% 
  #   pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
  #   separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
  #   pivot_wider(names_from = PFF.status, values_from = log2FC) %>%
  #   filter(PFFno<0, PFFyes<0) %>%
  #   mutate(delta = PFFno - PFFyes) %>%
  #   mutate(timepoint_comparison=factor(timepoint_comparison, levels = paste0("h",c(1,3,6,12),"_vs_h0"))) %>%
  #   ggplot(aes(x=delta)) + 
  #   #geom_histogram(fill="white", alpha=0.5, position="identity", bins = 100) +
  #   geom_histogram(aes(y=..density..), colour="black", position="identity", alpha=0.3, bins=200)+
  #   geom_density(alpha=.2) +
  #   facet_wrap(vars(timepoint_comparison), ncol = 1, nrow = 4, strip.position = "top") +
  #   ggtitle(paste0("DEresult.decayDiff.histogram.paired.", condition ,".all.pdf")) + 
  #   ggsave(file.path(output_dir, paste0("DEresult.decayDiff.histogram.paired.", condition ,".all.pdf")))
  
  ## count plot for individual significant genes
  cnts <- t(counts(dds_subset, normalized = T, replaced = F)[rownames(res)[which(res$padj<=0.05)], ] + 0.5)
  intgroup = c("Timepoint","PFF.status")
  dd <- data.frame(cnts, colData(dds_subset)[rownames(cnts), intgroup])
  dd <- pivot_longer(dd, cols=starts_with("ENSG"), names_to = "geneID", values_to = "count") %>% 
    mutate(geneSymbol=genes_annotation$geneSymbol[match(geneID, genes_annotation$geneID)])
  ggplot(dd, aes(x = Timepoint, y = count, color = PFF.status, group = PFF.status)) + 
    geom_point() + stat_summary(fun.y=mean, geom="line") +
    labs(title=paste0("Normalized counts for genes with significant condition-specific changes over time (LRT padj < 0.05)"), 
         subtitle=file.path(output_dir, paste0("DEresult.interactionLRT.", condition ,".padj05.pdf")), y="Normalized count (log10 scalue)") +
    scale_y_log10() + facet_wrap_paginate(vars(geneSymbol), scales = "free_y", ncol = 3, nrow = 5, page = 1, strip.position = "top") +
    theme(strip.background = element_blank(), strip.placement = "outside") -> p
  required_n_pages <- n_pages(p)
  pdf(file.path(output_dir, paste0("DEresult.interactionLRT.", condition ,".padj05.plotCounts.pdf")), width=8, height = 10, paper = 'letter')
  for(i in 1:required_n_pages){
    ggplot(dd, aes(x = Timepoint, y = count, color = PFF.status, group = PFF.status)) + 
      geom_point() + stat_summary(fun.y=mean, geom="line") +
      labs(title=paste0("Genes with significant condition-specific changes over time (LRT padj < 0.05)"), caption = paste("page",i,"of",required_n_pages),
           subtitle=file.path(output_dir, paste0("DEresult.interactionLRT.", condition ,".padj05.pdf")), y="Normalized count (log10 scalue)") +
      scale_y_log10() + facet_wrap_paginate(vars(geneSymbol), scales = "free_y", ncol = 3, nrow = 5, page =i, strip.position = "top") +
      theme(strip.background = element_blank(), strip.placement = "outside") -> p
    print(p)
  }
  dev.off()
  
  ## save the normalized count to a file
  df = dd %>% group_by(Timepoint, PFF.status, geneID, geneSymbol) %>% 
    summarise(mean_normalized_count=mean(count)) %>% ungroup() %>%
    pivot_wider(names_from=Timepoint, values_from=mean_normalized_count) %>%
    mutate(log2FC_h1_vs_h0 = log2(h1/h0), 
           log2FC_h3_vs_h1 = log2(h3/h1),
           log2FC_h6_vs_h3 = log2(h6/h3),
           log2FC_h12_vs_h6= log2(h12/h6)) %>% 
    arrange(geneID, PFF.status)
  write.table(df, 
              file=file.path(output_dir, paste0("DEresult.interactionLRT.", condition ,".padj05.consective_log2FC.xls")), 
              sep="\t", quote =F, na="", row.names = F)
  
  write.table(t(cnts),  
              file=file.path(output_dir, paste0("DEresult.interactionLRT.", condition ,".padj05.normalizedCount.xls")), 
              sep="\t", quote =F, na="", row.names=T, col.names = NA)
  
  ## now we can cluster the top significant genes (in any timepoint) by their Log2FC profiles.
  topN=100
  topDE = filter(res[1:topN,], padj<=0.05) %>% as_tibble() %>% column_to_rownames('symbol') %>%  select(contains("log2FC"))
  colnames(topDE) = gsub("log2FC.","",colnames(topDE))
  ## trim max and min
  MAX=3; MIN=-3; topDE[topDE>MAX]=MAX;topDE[topDE<MIN]=MIN
  
  ## to make sure the 0-point on the color scale is white (ref: https://stackoverflow.com/a/31707976)
  paletteLength <- 50
  myColor = colorRampPalette(c("blue", "white", "red"))(paletteLength)
  # use floor and ceiling to deal with even/odd length pallettelengths
  myBreaks = c(seq(MIN, 0, length.out=ceiling(paletteLength/2) + 1), 
               seq(MAX/paletteLength, MAX, length.out=floor(paletteLength/2)))
  
  
  annotation_row = filter(res[1:topN,], padj<=0.05) %>% as_tibble() %>% select(geneType, symbol)  %>% column_to_rownames('symbol')
  annotation_col = data.frame(rowname=names(topDE), PFF.status=sub("\\..*","",names(topDE))) %>% column_to_rownames()
  ann_colors = list(
    #geneType = c(protein_coding = "darkblue", lincRNA = "orange", Mt_rRNA='pink', antisense='yellow', pseudogene='gray', processed_transcript='lightgray',Mt_tRNA='purple',misc_RNA='red',snoRNA='lightblue', polymorphic_pseudogene="gray", processed_pseudogene="gray", TEC="lightred"),
    PFF.status = c(PFFno = "#F8766D", PFFyes = "#00BFC4")
  )
  par(cex=0.5, mar=c(5, 8, 4, 1))
  pheatmap(topDE,
           fontsize = 8,
           main =paste("Heatmap of log2 fold changes for genes with adjusted p < 0.05\n",
                       file.path(output_dir, paste0("DEresult.interactionLRT.", condition ,".all.xls.gz"))),
           width=6, height=0.3*nrow(topDE),
           filename=file.path(output_dir, paste0("DEresult.interactionLRT.", condition ,".topP.heatmap.pdf")),
           #border_color = NA,
           color=myColor, breaks=myBreaks,
           annotation_row = annotation_row,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           drop_levels = TRUE,
           scale = "none", 
           #clustering_method = 'ward.D', 
           cluster_rows = TRUE,
           #clustering_distance_rows = "correlation",
           cluster_cols = F)
  
}

###########################################
# run DE among CONDITION, Timepoint jointly
# to test any genes that react in a condition-specific manner over time
# for consecutive subset 
###########################################

for(X in levels(covarianceTable$PFF.status)){
  message(paste("# Running DEseq on PFF status =", X, "..."))
  # subsetting
  # debug: X="PFFyes"
  Timepoints_sorted_levels = sort(unique(colData(dds)$Timepoint))
  
  # loop thru all levels and test consecutive comparisons
  for(i in 2:length(Timepoints_sorted_levels)){
    time_i = as.character(Timepoints_sorted_levels[i])
    time_i_m_1 = as.character(Timepoints_sorted_levels[i-1])
    message(paste("# Running LRT between consecutive ",time_i,"vs",time_i_m_1, "..."))
    
    dds_subset=dds[,dds$PFF.status==X & dds$Timepoint %in% c(time_i_m_1,time_i)]; 
    dim(dds_subset);
    
    # re-factorize 
    colData(dds_subset)$Timepoint = factor(colData(dds_subset)$Timepoint, levels = c(time_i_m_1,time_i))
    
    design(dds_subset) <- as.formula(" ~ CONDITION + Timepoint + CONDITION:Timepoint")
    
    # a likelihood ratio test, where we remove the condition-specific differences over time. 
    # Genes with small p values from this test are those which at one or more time points after time 0 showed a condition-specific effect. 
    # Note therefore that this will not give small p values to genes that moved up or down over time in the same way in both strains.
    dds_subset <- DESeq(dds_subset, test="LRT", reduced = ~ CONDITION + Timepoint,
                        parallel=TRUE, BPPARAM=MulticoreParam(4))
    
    resultsNames(dds_subset)
    # This will give p-values for LRT test: genes which at one or more time points after time 0 showed a strain-specific effect. 
    # Note that the log2FC is for the interaction term (e.g. "CONDITIONSNCA4.Timepointh12" by default or by the designed name, e.g. name = "CONDITIONSNCA4.Timepointh3"). 
    res <- results(dds_subset, alpha = 0.05, parallel=TRUE, BPPARAM=MulticoreParam(4))  ## Note that the alpha for indepedent filtering has to be same as the number to threshold res$padj (https://support.bioconductor.org/p/71732/#71771)
    
    # extract a matrix of the log2 fold changes for all comparisons using the coef function.
    betas <- coef(dds_subset)
    colnames(betas)  # same as the resultsNames(dds_subset)
    betas = betas[rownames(res), -c(1,2)]
    ## Note that CONDITIONSNCA4.Timepointh1 means difference between SNCA4 vs SNCA2 at h1, controlling for baseline. 
    ## To get the log2 fold change of h1 vs h0 for the SNCA2, Timepoint_h1_vs_h0
    ## To get the log2 fold change of h1 vs h0 for the SNCA4, Timepoint_h1_vs_h0 +  CONDITIONSNCA4.Timepointh1
    var_snca2 =  paste("Timepoint",time_i, "vs", time_i_m_1, sep="_")
    var_delta =  paste0("CONDITIONSNCA4.Timepoint", time_i)
    betas = as.data.frame(betas) %>% mutate(.keep = "none",
                                            "SNCA2.{time_i}_vs_{time_i_m_1}"  := (!!as.name(var_snca2)), # see https://stackoverflow.com/a/26003971 for dynamic variable names in dplyr
                                            "SNCA4.{time_i}_vs_{time_i_m_1}"  := (!!as.name(var_snca2)) + (!!as.name(var_delta))) # see https://stackoverflow.com/a/29678662
    
    colnames(betas) = paste0("log2FC.",colnames(betas)) 
    
    # combine the overall LRT tesst (for any timepoint) and the log2FC for each timepoint together
    res = cbind(as.data.frame(res)[,c('baseMean', 'stat','pvalue','padj')], betas)
    
    # add annotation
    res$symbol <- genes_annotation$geneSymbol[match(row.names(res), genes_annotation$geneID)]
    res$geneType <- genes_annotation$geneType[match(row.names(res), genes_annotation$geneID)]
    
    ## save the data
    res <- res[order(res$padj),]
    head(res); dim(res)
    
    write.table(as.data.frame(res), 
                file=gzfile(file.path(output_dir, paste0("DEresult.interactionLRT.", X , ".consecutive_",time_i,"_vs_",time_i_m_1,".all.xls.gz"))), 
                sep="\t", quote =F, na="", row.names=T, col.names = NA)
    
    write.table(subset(res, padj<=0.05), 
                file=file.path(output_dir, paste0("DEresult.interactionLRT.", X ,".consecutive_", time_i, "_vs_", time_i_m_1, ".padj05.xls")), 
                sep="\t", quote =F, na="", row.names=T, col.names = NA)
    
    ## any globe difference in slope between SNCA2 vs. SNCA4
    tests_pvalue = rownames_to_column(res) %>% 
      pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
      separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
      pivot_wider(names_from = CONDITION, values_from = log2FC) %>% 
      filter(pvalue<=1) %>%
      group_by(timepoint_comparison) %>% 
      summarise(SNCA2=list(SNCA2), SNCA4=list(SNCA4)) %>% rowwise() %>% 
      mutate(pvalue.paired_ttest_less=t.test(unlist(SNCA2), unlist(SNCA4), paired = T, alternative = "less")$p.value) %>% 
      mutate(pvalue.paired_ttest=t.test(unlist(SNCA2), unlist(SNCA4), paired = T)$p.value) %>% 
      mutate(pvalue.ttest_less=t.test(unlist(SNCA2), unlist(SNCA4), alternative = "less")$p.value) %>% 
      mutate(pvalue.ttest=t.test(unlist(SNCA2), unlist(SNCA4))$p.value) %>% 
      mutate(pvalue.paired_wilcox_less=wilcox.test(unlist(SNCA2), unlist(SNCA4), paired = T, alternative = "less")$p.value) %>% 
      mutate(pvalue.paired_wilcox=wilcox.test(unlist(SNCA2), unlist(SNCA4), paired = T)$p.value) %>% 
      mutate(pvalue.wilcox_less=wilcox.test(unlist(SNCA2), unlist(SNCA4), alternative = "less")$p.value) %>% 
      mutate(pvalue.wilcox=wilcox.test(unlist(SNCA2), unlist(SNCA4))$p.value) 
    ttest_pvalue = pivot_longer(tests_pvalue, cols = starts_with("pvalue"), names_to = "tests", names_prefix="pvalue.", values_to="pvalue") %>% filter(grepl("ttest", tests))
    
    nDecay = rownames_to_column(res) %>% 
      pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
      separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
      pivot_wider(names_from = CONDITION, values_from = log2FC) %>% 
      filter(pvalue<=1) %>%
      group_by(timepoint_comparison) %>% 
      summarise(SNCA2.percent=mean(SNCA2<0), SNCA4.percent=mean(SNCA4<0), SNCA2.n=sum(SNCA2<0), SNCA4.n=sum(SNCA4<0)) %>% 
      pivot_longer(cols = contains('SNCA'), names_to = c('CONDITION',".value"), names_pattern = "(.*)\\.(.*)") 
    
    p=pivot_longer(as_tibble(res), starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
      separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
      filter(pvalue<=1) %>%
      ggplot(aes(x=log2FC, fill=CONDITION, color=CONDITION)) + 
      #geom_histogram(fill="white", alpha=0.5, position="identity", bins = 100) +
      geom_histogram(aes(y=..density..), colour="black", position="identity", alpha=0.3, bins=100)+
      geom_density(alpha=.2) +
      geom_vline(xintercept =0, color='blue', linetype=2) +
      #facet_wrap(vars(timepoint_comparison), ncol = 1, nrow = 4, strip.position = "top") + 
      geom_text(data = nDecay, mapping = aes(x = -Inf, y = Inf, label = paste0("Decreasing genes in ",CONDITION, ": n=",n, " (",round(100*percent, 2), "%)"), color = CONDITION), hjust   = -0.05, vjust   = rep(c(1.5, 3), length(unique(nDecay$timepoint_comparison)))) +
      geom_text(data = ttest_pvalue, mapping = aes(x = -Inf, y = Inf, label = paste(tests, "p = ",formatC(pvalue, format = "g", digits = 2)), color = NULL,fill= NULL), hjust   = -0.05, vjust   = rep(c(4.5, 6, 7.5, 9), length(unique(nDecay$timepoint_comparison)))) +
      ggtitle(paste0("DEresult.interactionLRT.", X ,".consecutive_", time_i, "_vs_", time_i_m_1, ".log2FC.histogram.grouped.all.pdf")) 
    ggsave(file.path(output_dir, paste0("DEresult.interactionLRT.", X ,".consecutive_", time_i, "_vs_", time_i_m_1, ".log2FC.histogram.grouped.all.pdf")), width = 11, height = 3)
    
    # histogram for negative slope genes in both conditions
    var_snca2 =  paste0("log2FC.SNCA2.",time_i, "_vs_", time_i_m_1)
    var_snca4 =  paste0("log2FC.SNCA4.",time_i, "_vs_", time_i_m_1)
    res_down_down = filter(res, (!!as.name(var_snca2)) < 0, (!!as.name(var_snca4)) < 0)
    tests_pvalue = rownames_to_column(res_down_down) %>% 
      pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
      separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
      pivot_wider(names_from = CONDITION, values_from = log2FC) %>% 
      filter(pvalue<=1) %>%
      group_by(timepoint_comparison) %>% 
      summarise(SNCA2=list(SNCA2), SNCA4=list(SNCA4)) %>% rowwise() %>% 
      mutate(pvalue.paired_ttest_less=t.test(unlist(SNCA2), unlist(SNCA4), paired = T, alternative = "less")$p.value) %>% 
      mutate(pvalue.paired_ttest=t.test(unlist(SNCA2), unlist(SNCA4), paired = T)$p.value) %>% 
      mutate(pvalue.ttest_less=t.test(unlist(SNCA2), unlist(SNCA4), alternative = "less")$p.value) %>% 
      mutate(pvalue.ttest=t.test(unlist(SNCA2), unlist(SNCA4))$p.value) %>% 
      mutate(pvalue.paired_wilcox_less=wilcox.test(unlist(SNCA2), unlist(SNCA4), paired = T, alternative = "less")$p.value) %>% 
      mutate(pvalue.paired_wilcox=wilcox.test(unlist(SNCA2), unlist(SNCA4), paired = T)$p.value) %>% 
      mutate(pvalue.wilcox_less=wilcox.test(unlist(SNCA2), unlist(SNCA4), alternative = "less")$p.value) %>% 
      mutate(pvalue.wilcox=wilcox.test(unlist(SNCA2), unlist(SNCA4))$p.value) 
    ttest_pvalue = pivot_longer(tests_pvalue, cols = starts_with("pvalue"), names_to = "tests", names_prefix="pvalue.", values_to="pvalue") %>% filter(grepl("ttest", tests))
    
    nDecay = rownames_to_column(res_down_down) %>% 
      pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
      separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
      pivot_wider(names_from = CONDITION, values_from = log2FC) %>% 
      filter(pvalue<=1) %>%
      group_by(timepoint_comparison) %>% 
      summarise(SNCA2.percent=mean(SNCA2<0), SNCA4.percent=mean(SNCA4<0), SNCA2.n=sum(SNCA2<0), SNCA4.n=sum(SNCA4<0)) %>% 
      pivot_longer(cols = contains('SNCA'), names_to = c('CONDITION',".value"), names_pattern = "(.*)\\.(.*)") 
    
    p=pivot_longer(as_tibble(res_down_down), starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
      separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
      filter(pvalue<=1) %>%
      ggplot(aes(x=log2FC, fill=CONDITION, color=CONDITION)) + 
      #geom_histogram(fill="white", alpha=0.5, position="identity", bins = 100) +
      geom_histogram(aes(y=..density..), colour="black", position="identity", alpha=0.3, bins=100)+
      geom_density(alpha=.2) +
      geom_vline(xintercept =0, color='blue', linetype=2) +
      facet_wrap(vars(timepoint_comparison), ncol = 1, nrow = 4, strip.position = "top") + 
      geom_text(data = nDecay, mapping = aes(x = -Inf, y = Inf, label = paste0("Decreasing genes in ",CONDITION, ": n=",n, " (",round(100*percent, 2), "%)"), color = CONDITION), hjust   = -0.05, vjust   = rep(c(1.5, 3), length(unique(nDecay$timepoint_comparison)))) +
      geom_text(data = ttest_pvalue, mapping = aes(x = -Inf, y = Inf, label = paste(tests, "p = ",formatC(pvalue, format = "g", digits = 2)), color = NULL,fill= NULL), hjust   = -.05, vjust   = rep(c(4.5, 6, 7.5, 9), length(unique(nDecay$timepoint_comparison)))) +
      ggtitle(paste0("DEresult.interactionLRT.", X ,".consecutive_", time_i, "_vs_", time_i_m_1, ".log2FC.histogram.grouped.downdown.pdf")) 
    ggsave(file.path(output_dir, paste0("DEresult.interactionLRT.", X ,".consecutive_", time_i, "_vs_", time_i_m_1, ".log2FC.histogram.grouped.downdown.pdf")), width = 11, height = 3)
    
    # if the decreasing ones have bigger decay rate in SNCA2 vs. SNCA4
    p=rownames_to_column(res_down_down) %>% 
      mutate(SNCA2_minus_SNCA4 = (!!as.name(var_snca2)) - (!!as.name(var_snca4)))  
    p=ggplot(p, aes(x=SNCA2_minus_SNCA4, fill=ifelse(SNCA2_minus_SNCA4<0, paste0("SNCA2<SNCA4: ",round(100*mean(SNCA2_minus_SNCA4<0),2),"%"), paste0("SNCA2>SNCA4: ",round(100*mean(SNCA2_minus_SNCA4>0),2),"%")))) + 
      geom_histogram(colour="black", position="identity", alpha=0.3, breaks=seq(-3,3,0.1)) +
      xlab(paste(var_snca2, "-", var_snca4)) + 
      labs(fill=paste("t.test (less) p = ",formatC(t.test(p$SNCA2_minus_SNCA4, alternative = "less")$p.value, format = "g", digits = 2))) + theme(legend.position="top") + 
      ggtitle(paste(time_i, "vs.", time_i_m_1, "slope difference between SNCA2 and SNCA4 for the decreasing genes (n =", nrow(res_down_down),")"))
    ggsave(file.path(output_dir, paste0("DEresult.interactionLRT.", X ,".consecutive_", time_i, "_vs_", time_i_m_1, ".log2FC.histogram.delta.downdown.pdf")), width = 11, height = 3)
    
    ## up up genes
    res_up_up = filter(res, (!!as.name(var_snca2)) > 0, (!!as.name(var_snca4)) > 0)
    tests_pvalue = rownames_to_column(res_up_up) %>% 
      pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
      separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
      pivot_wider(names_from = CONDITION, values_from = log2FC) %>% 
      filter(pvalue<=1) %>%
      group_by(timepoint_comparison) %>% 
      summarise(SNCA2=list(SNCA2), SNCA4=list(SNCA4)) %>% rowwise() %>% 
      mutate(pvalue.paired_ttest_less=t.test(unlist(SNCA2), unlist(SNCA4), paired = T, alternative = "less")$p.value) %>% 
      mutate(pvalue.paired_ttest=t.test(unlist(SNCA2), unlist(SNCA4), paired = T)$p.value) %>% 
      mutate(pvalue.ttest_less=t.test(unlist(SNCA2), unlist(SNCA4), alternative = "less")$p.value) %>% 
      mutate(pvalue.ttest=t.test(unlist(SNCA2), unlist(SNCA4))$p.value) %>% 
      mutate(pvalue.paired_wilcox_less=wilcox.test(unlist(SNCA2), unlist(SNCA4), paired = T, alternative = "less")$p.value) %>% 
      mutate(pvalue.paired_wilcox=wilcox.test(unlist(SNCA2), unlist(SNCA4), paired = T)$p.value) %>% 
      mutate(pvalue.wilcox_less=wilcox.test(unlist(SNCA2), unlist(SNCA4), alternative = "less")$p.value) %>% 
      mutate(pvalue.wilcox=wilcox.test(unlist(SNCA2), unlist(SNCA4))$p.value) 
    ttest_pvalue = pivot_longer(tests_pvalue, cols = starts_with("pvalue"), names_to = "tests", names_prefix="pvalue.", values_to="pvalue") %>% filter(grepl("ttest", tests))
    
    nDecay = rownames_to_column(res_up_up) %>% 
      pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
      separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
      pivot_wider(names_from = CONDITION, values_from = log2FC) %>% 
      filter(pvalue<=1) %>%
      group_by(timepoint_comparison) %>% 
      summarise(SNCA2.percent=mean(SNCA2>0), SNCA4.percent=mean(SNCA4>0), SNCA2.n=sum(SNCA2>0), SNCA4.n=sum(SNCA4>0)) %>% 
      pivot_longer(cols = contains('SNCA'), names_to = c('CONDITION',".value"), names_pattern = "(.*)\\.(.*)") 
    
    p=pivot_longer(as_tibble(res_up_up), starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
      separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
      filter(pvalue<=1) %>%
      ggplot(aes(x=log2FC, fill=CONDITION, color=CONDITION)) + 
      #geom_histogram(fill="white", alpha=0.5, position="identity", bins = 100) +
      geom_histogram(aes(y=..density..), colour="black", position="identity", alpha=0.3, bins=100)+
      geom_density(alpha=.2) +
      geom_vline(xintercept =0, color='blue', linetype=2) +
      facet_wrap(vars(timepoint_comparison), ncol = 1, nrow = 4, strip.position = "top") + 
      geom_text(data = nDecay, mapping = aes(x = Inf, y = Inf, label = paste0("Increasing genes in ",CONDITION, ": n=",n, " (",round(100*percent, 2), "%)"), color = CONDITION), hjust   = 1.05, vjust   = rep(c(1.5, 3), length(unique(nDecay$timepoint_comparison)))) +
      geom_text(data = ttest_pvalue, mapping = aes(x = Inf, y = Inf, label = paste(tests, "p = ",formatC(pvalue, format = "g", digits = 2)), color = NULL,fill= NULL), hjust   = 1.05, vjust   = rep(c(4.5, 6, 7.5, 9), length(unique(nDecay$timepoint_comparison)))) +
      ggtitle(paste0("DEresult.interactionLRT.", X ,".consecutive_", time_i, "_vs_", time_i_m_1, ".log2FC.histogram.grouped.upup.pdf")) 
    ggsave(file.path(output_dir, paste0("DEresult.interactionLRT.", X ,".consecutive_", time_i, "_vs_", time_i_m_1, ".log2FC.histogram.grouped.upup.pdf")), width = 11, height = 3)
    
    # if the decreasing ones have bigger decay rate in SNCA2 vs. SNCA4
    p=rownames_to_column(res_up_up) %>% 
      mutate(SNCA2_minus_SNCA4 = (!!as.name(var_snca2)) - (!!as.name(var_snca4)))  
    p=ggplot(p, aes(x=SNCA2_minus_SNCA4, fill=ifelse(SNCA2_minus_SNCA4<0, paste0("SNCA2<SNCA4: ",round(100*mean(SNCA2_minus_SNCA4<0),2),"%"), paste0("SNCA2>SNCA4: ",round(100*mean(SNCA2_minus_SNCA4>0),2),"%")))) + 
      geom_histogram(colour="black", position="identity", alpha=0.3, breaks=seq(-3,3,0.1)) +
      xlab(paste(var_snca2, "-", var_snca4)) + 
      labs(fill=paste("t.test (less) p = ",formatC(t.test(p$SNCA2_minus_SNCA4, alternative = "less")$p.value, format = "g", digits = 2))) + theme(legend.position="top") + 
      ggtitle(paste(time_i, "vs.", time_i_m_1, "slope difference between SNCA2 and SNCA4 for the increasing genes (n =", nrow(res_up_up),")"))
    ggsave(file.path(output_dir, paste0("DEresult.interactionLRT.", X ,".consecutive_", time_i, "_vs_", time_i_m_1, ".log2FC.histogram.delta.upup.pdf")), width = 11, height = 3)
    
    ## now we can cluster the top significant genes by their Log2FC profiles.
    topN=1000
    topDE = filter(res, padj<=0.05) %>% head(topN) %>% as_tibble() %>% column_to_rownames('symbol') %>%  select(contains("log2FC")) %>% arrange((!!as.name(var_snca2)))
    colnames(topDE) = gsub("log2FC.","",colnames(topDE))
    dim(topDE)
    ## trim max and min
    MAX=2; MIN=-2; topDE[topDE>MAX]=MAX;topDE[topDE<MIN]=MIN
    
    ## to make sure the 0-point on the color scale is white (ref: https://stackoverflow.com/a/31707976)
    paletteLength <- 50
    myColor = colorRampPalette(c("blue", "white", "red"))(paletteLength)
    # use floor and ceiling to deal with even/odd length pallettelengths
    myBreaks = c(seq(MIN, 0, length.out=ceiling(paletteLength/2) + 1), 
                 seq(MAX/paletteLength, MAX, length.out=floor(paletteLength/2)))
    
    
    annotation_row = filter(res[1:topN,], padj<=0.05) %>% as_tibble() %>% select(geneType, symbol)  %>% column_to_rownames('symbol')
    annotation_col = data.frame(rowname=names(topDE), CONDITION=sub("\\..*","",names(topDE))) %>% column_to_rownames()
    ann_colors = list(
      #geneType = c(protein_coding = "darkblue", lincRNA = "orange", Mt_rRNA='pink', antisense='yellow', pseudogene='gray', processed_transcript='lightgray',Mt_tRNA='purple',misc_RNA='red',snoRNA='lightblue', polymorphic_pseudogene="gray", processed_pseudogene="gray", TEC="lightred"),
      CONDITION = c(SNCA2 = "#F8766D", SNCA4 = "#00BFC4")
    )
    if(nrow(topDE)>0){
      par(cex=0.5, mar=c(5, 8, 4, 1))
      pheatmap(t(topDE),
               fontsize = 8,
               main =paste("Heatmap of log2FC for genes with LRT FDR<0.05\n",
                           file.path(output_dir, paste0("DEresult.interactionLRT.", X ,".consecutive_", time_i, "_vs_", time_i_m_1, ".padj05.xls"))),
               width=3+0.15*nrow(topDE), height=2,
               filename=file.path(output_dir, paste0("DEresult.interactionLRT.", X ,".consecutive_", time_i, "_vs_", time_i_m_1, ".topP.heatmap.pdf")),
               #border_color = NA,
               color=myColor, breaks=myBreaks,
               annotation_row = annotation_col,
               annotation_col = annotation_row,
               annotation_colors = ann_colors,
               drop_levels = TRUE,
               scale = "none", border_color = F,
               #clustering_method = 'ward.D', 
               cluster_rows = F,
               #clustering_distance_rows = "correlation",
               cluster_cols = F)
    }
    
  }
}

###########################################
# run DE among PFF, Timepoint jointly
# to test any genes that react in a condition-specific manner over time
# see http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments
# for consecutive subset only
###########################################

for(X in levels(covarianceTable$CONDITION)){
  message(paste("# Running DEseq on CONDITION status =", X, "..."))
  # subsetting
  # debug: X="SNCA2"
  
  Timepoints_sorted_levels = sort(unique(colData(dds)$Timepoint))
  # loop thru all levels and test consecutive comparisons
  for(i in 2:length(Timepoints_sorted_levels)){
    time_i = as.character(Timepoints_sorted_levels[i])
    time_i_m_1 = as.character(Timepoints_sorted_levels[i-1])
    message(paste("# Running LRT between consecutive ",time_i,"vs",time_i_m_1, "..."))
    
    dds_subset=dds[,dds$CONDITION==X & dds$Timepoint %in% c(time_i_m_1,time_i)]; 
    dim(dds_subset);
    
    # re-factorize 
    colData(dds_subset)$Timepoint = factor(colData(dds_subset)$Timepoint, levels = c(time_i_m_1,time_i))
    
    design(dds_subset) <- as.formula(" ~ PFF.status + Timepoint + PFF.status:Timepoint")
    
    # a likelihood ratio test, where we remove the PFF.status-specific differences over time. 
    # Genes with small p values from this test are those which at one or more time points after time 0 showed a PFF.status-specific effect. 
    # Note therefore that this will not give small p values to genes that moved up or down over time in the same way in both strains.
    dds_subset <- DESeq(dds_subset, test="LRT", reduced = ~ PFF.status + Timepoint,
                        parallel=TRUE, BPPARAM=MulticoreParam(4))
    
    resultsNames(dds_subset)
    # This will give p-values for LRT test: genes which at one or more time points after time 0 showed a strain-specific effect. 
    # Note that the log2FC is for the interaction term (e.g. "PFF.statusPFFyes.Timepointh12" by default or by the designed name, e.g. name = "PFF.statusPFFyes.Timepointh3"). 
    res <- results(dds_subset, alpha = 0.05, parallel=TRUE, BPPARAM=MulticoreParam(4))  ## Note that the alpha for indepedent filtering has to be same as the number to threshold res$padj (https://support.bioconductor.org/p/71732/#71771)
    
    # extract a matrix of the log2 fold changes for all comparisons using the coef function.
    betas <- coef(dds_subset)
    colnames(betas)  # same as the resultsNames(dds_subset)
    betas = betas[rownames(res), -c(1,2)]
    ## Note that PFF.statusPFFyes.Timepointh1 means difference between PFFyes vs PFFno at h1, controlling for baseline. 
    ## To get the log2 fold change of h1 vs h0 for the PFFno, Timepoint_h1_vs_h0
    ## To get the log2 fold change of h1 vs h0 for the PFFyes, Timepoint_h1_vs_h0 +  PFF.statusPFFyes.Timepointh1
    var_snca2 =  paste("Timepoint",time_i, "vs", time_i_m_1, sep="_")
    var_delta =  paste0("PFF.statusPFFyes.Timepoint", time_i)
    betas = as.data.frame(betas) %>% mutate(.keep = "none",
                                            "PFFno.{time_i}_vs_{time_i_m_1}"  := (!!as.name(var_snca2)), # see https://stackoverflow.com/a/26003971 for dynamic variable names in dplyr
                                            "PFFyes.{time_i}_vs_{time_i_m_1}"  := (!!as.name(var_snca2)) + (!!as.name(var_delta))) # see https://stackoverflow.com/a/29678662
    
    colnames(betas) = paste0("log2FC.",colnames(betas)) 
    
    # combine the overall LRT tesst (for any timepoint) and the log2FC for each timepoint together
    res = cbind(as.data.frame(res)[,c('baseMean', 'stat','pvalue','padj')], betas)
    
    # add annotation
    res$symbol <- genes_annotation$geneSymbol[match(row.names(res), genes_annotation$geneID)]
    res$geneType <- genes_annotation$geneType[match(row.names(res), genes_annotation$geneID)]
    
    ## save the data
    res <- res[order(res$padj),]
    head(res); dim(res)
    
    write.table(as.data.frame(res), 
                file=gzfile(file.path(output_dir, paste0("DEresult.interactionLRT.", X , ".consecutive_",time_i,"_vs_",time_i_m_1,".all.xls.gz"))), 
                sep="\t", quote =F, na="", row.names=T, col.names = NA)
    
    write.table(subset(res, padj<=0.05), 
                file=file.path(output_dir, paste0("DEresult.interactionLRT.", X ,".consecutive_", time_i, "_vs_", time_i_m_1, ".padj05.xls")), 
                sep="\t", quote =F, na="", row.names=T, col.names = NA)
    
    ## any globe difference in slope between PFFno vs. PFFyes
    tests_pvalue = rownames_to_column(res) %>% 
      pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
      separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
      pivot_wider(names_from = PFF.status, values_from = log2FC) %>% 
      filter(pvalue<=1) %>%
      group_by(timepoint_comparison) %>% 
      summarise(PFFno=list(PFFno), PFFyes=list(PFFyes)) %>% rowwise() %>% 
      mutate(pvalue.paired_ttest_less=t.test(unlist(PFFno), unlist(PFFyes), paired = T, alternative = "less")$p.value) %>% 
      mutate(pvalue.paired_ttest=t.test(unlist(PFFno), unlist(PFFyes), paired = T)$p.value) %>% 
      mutate(pvalue.ttest_less=t.test(unlist(PFFno), unlist(PFFyes), alternative = "less")$p.value) %>% 
      mutate(pvalue.ttest=t.test(unlist(PFFno), unlist(PFFyes))$p.value) %>% 
      mutate(pvalue.paired_wilcox_less=wilcox.test(unlist(PFFno), unlist(PFFyes), paired = T, alternative = "less")$p.value) %>% 
      mutate(pvalue.paired_wilcox=wilcox.test(unlist(PFFno), unlist(PFFyes), paired = T)$p.value) %>% 
      mutate(pvalue.wilcox_less=wilcox.test(unlist(PFFno), unlist(PFFyes), alternative = "less")$p.value) %>% 
      mutate(pvalue.wilcox=wilcox.test(unlist(PFFno), unlist(PFFyes))$p.value) 
    ttest_pvalue = pivot_longer(tests_pvalue, cols = starts_with("pvalue"), names_to = "tests", names_prefix="pvalue.", values_to="pvalue") %>% filter(grepl("ttest", tests))
    
    nDecay = rownames_to_column(res) %>% 
      pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
      separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
      pivot_wider(names_from = PFF.status, values_from = log2FC) %>% 
      filter(pvalue<=1) %>%
      group_by(timepoint_comparison) %>% 
      summarise(PFFno.percent=mean(PFFno<0), PFFyes.percent=mean(PFFyes<0), PFFno.n=sum(PFFno<0), PFFyes.n=sum(PFFyes<0)) %>% 
      pivot_longer(cols = contains('PFF'), names_to = c('PFF.status',".value"), names_pattern = "(.*)\\.(.*)") 
    
    p=pivot_longer(as_tibble(res), starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
      separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
      filter(pvalue<=1) %>%
      ggplot(aes(x=log2FC, fill=PFF.status, color=PFF.status)) + 
      #geom_histogram(fill="white", alpha=0.5, position="identity", bins = 100) +
      geom_histogram(aes(y=..density..), colour="black", position="identity", alpha=0.3, bins=100)+
      geom_density(alpha=.2) +
      geom_vline(xintercept =0, color='blue', linetype=2) +
      #facet_wrap(vars(timepoint_comparison), ncol = 1, nrow = 4, strip.position = "top") + 
      geom_text(data = nDecay, mapping = aes(x = -Inf, y = Inf, label = paste0("Decreasing genes in ",PFF.status, ": n=",n, " (",round(100*percent, 2), "%)"), color = PFF.status), hjust   = -0.05, vjust   = rep(c(1.5, 3), length(unique(nDecay$timepoint_comparison)))) +
      geom_text(data = ttest_pvalue, mapping = aes(x = -Inf, y = Inf, label = paste(tests, "p = ",formatC(pvalue, format = "g", digits = 2)), color = NULL,fill= NULL), hjust   = -0.05, vjust   = rep(c(4.5, 6, 7.5, 9), length(unique(nDecay$timepoint_comparison)))) +
      ggtitle(paste0("DEresult.interactionLRT.", X ,".consecutive_", time_i, "_vs_", time_i_m_1, ".log2FC.histogram.grouped.all.pdf")) 
    ggsave(file.path(output_dir, paste0("DEresult.interactionLRT.", X ,".consecutive_", time_i, "_vs_", time_i_m_1, ".log2FC.histogram.grouped.all.pdf")), width = 11, height = 3)
    
    # histogram for negative slope genes in both conditions
    var_snca2 =  paste0("log2FC.PFFno.",time_i, "_vs_", time_i_m_1)
    var_snca4 =  paste0("log2FC.PFFyes.",time_i, "_vs_", time_i_m_1)
    res_down_down = filter(res, (!!as.name(var_snca2)) < 0, (!!as.name(var_snca4)) < 0)
    tests_pvalue = rownames_to_column(res_down_down) %>% 
      pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
      separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
      pivot_wider(names_from = PFF.status, values_from = log2FC) %>% 
      filter(pvalue<=1) %>%
      group_by(timepoint_comparison) %>% 
      summarise(PFFno=list(PFFno), PFFyes=list(PFFyes)) %>% rowwise() %>% 
      mutate(pvalue.paired_ttest_less=t.test(unlist(PFFno), unlist(PFFyes), paired = T, alternative = "less")$p.value) %>% 
      mutate(pvalue.paired_ttest=t.test(unlist(PFFno), unlist(PFFyes), paired = T)$p.value) %>% 
      mutate(pvalue.ttest_less=t.test(unlist(PFFno), unlist(PFFyes), alternative = "less")$p.value) %>% 
      mutate(pvalue.ttest=t.test(unlist(PFFno), unlist(PFFyes))$p.value) %>% 
      mutate(pvalue.paired_wilcox_less=wilcox.test(unlist(PFFno), unlist(PFFyes), paired = T, alternative = "less")$p.value) %>% 
      mutate(pvalue.paired_wilcox=wilcox.test(unlist(PFFno), unlist(PFFyes), paired = T)$p.value) %>% 
      mutate(pvalue.wilcox_less=wilcox.test(unlist(PFFno), unlist(PFFyes), alternative = "less")$p.value) %>% 
      mutate(pvalue.wilcox=wilcox.test(unlist(PFFno), unlist(PFFyes))$p.value) 
    ttest_pvalue = pivot_longer(tests_pvalue, cols = starts_with("pvalue"), names_to = "tests", names_prefix="pvalue.", values_to="pvalue") %>% filter(grepl("ttest", tests))
    
    nDecay = rownames_to_column(res_down_down) %>% 
      pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
      separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
      pivot_wider(names_from = PFF.status, values_from = log2FC) %>% 
      filter(pvalue<=1) %>%
      group_by(timepoint_comparison) %>% 
      summarise(PFFno.percent=mean(PFFno<0), PFFyes.percent=mean(PFFyes<0), PFFno.n=sum(PFFno<0), PFFyes.n=sum(PFFyes<0)) %>% 
      pivot_longer(cols = contains('PFF'), names_to = c('PFF.status',".value"), names_pattern = "(.*)\\.(.*)") 
    
    p=pivot_longer(as_tibble(res_down_down), starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
      separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
      filter(pvalue<=1) %>%
      ggplot(aes(x=log2FC, fill=PFF.status, color=PFF.status)) + 
      #geom_histogram(fill="white", alpha=0.5, position="identity", bins = 100) +
      geom_histogram(aes(y=..density..), colour="black", position="identity", alpha=0.3, bins=100)+
      geom_density(alpha=.2) +
      geom_vline(xintercept =0, color='blue', linetype=2) +
      facet_wrap(vars(timepoint_comparison), ncol = 1, nrow = 4, strip.position = "top") + 
      geom_text(data = nDecay, mapping = aes(x = -Inf, y = Inf, label = paste0("Decreasing genes in ",PFF.status, ": n=",n, " (",round(100*percent, 2), "%)"), color = PFF.status), hjust   = -0.05, vjust   = rep(c(1.5, 3), length(unique(nDecay$timepoint_comparison)))) +
      geom_text(data = ttest_pvalue, mapping = aes(x = -Inf, y = Inf, label = paste(tests, "p = ",formatC(pvalue, format = "g", digits = 2)), color = NULL,fill= NULL), hjust   = -.05, vjust   = rep(c(4.5, 6, 7.5, 9), length(unique(nDecay$timepoint_comparison)))) +
      ggtitle(paste0("DEresult.interactionLRT.", X ,".consecutive_", time_i, "_vs_", time_i_m_1, ".log2FC.histogram.grouped.downdown.pdf")) 
    ggsave(file.path(output_dir, paste0("DEresult.interactionLRT.", X ,".consecutive_", time_i, "_vs_", time_i_m_1, ".log2FC.histogram.grouped.downdown.pdf")), width = 11, height = 3)
    
    # if the decreasing ones have bigger decay rate in PFFno vs. PFFyes
    p=rownames_to_column(res_down_down) %>% 
      mutate(PFFno_minus_PFFyes = (!!as.name(var_snca2)) - (!!as.name(var_snca4)))  
    p=ggplot(p, aes(x=PFFno_minus_PFFyes, fill=ifelse(PFFno_minus_PFFyes<0, paste0("PFFno<PFFyes: ",round(100*mean(PFFno_minus_PFFyes<0),2),"%"), paste0("PFFno>PFFyes: ",round(100*mean(PFFno_minus_PFFyes>0),2),"%")))) + 
      geom_histogram(colour="black", position="identity", alpha=0.3, breaks=seq(-3,3,0.1)) +
      xlab(paste(var_snca2, "-", var_snca4)) + 
      labs(fill=paste("t.test (less) p = ",formatC(t.test(p$PFFno_minus_PFFyes, alternative = "less")$p.value, format = "g", digits = 2))) + theme(legend.position="top") + 
      ggtitle(paste(time_i, "vs.", time_i_m_1, "slope difference between PFFno and PFFyes for the decreasing genes (n =", nrow(res_down_down),")"))
    ggsave(file.path(output_dir, paste0("DEresult.interactionLRT.", X ,".consecutive_", time_i, "_vs_", time_i_m_1, ".log2FC.histogram.delta.downdown.pdf")), width = 11, height = 3)
    
    ## up up genes
    res_up_up = filter(res, (!!as.name(var_snca2)) > 0, (!!as.name(var_snca4)) > 0)
    tests_pvalue = rownames_to_column(res_up_up) %>% 
      pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
      separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
      pivot_wider(names_from = PFF.status, values_from = log2FC) %>% 
      filter(pvalue<=1) %>%
      group_by(timepoint_comparison) %>% 
      summarise(PFFno=list(PFFno), PFFyes=list(PFFyes)) %>% rowwise() %>% 
      mutate(pvalue.paired_ttest_less=t.test(unlist(PFFno), unlist(PFFyes), paired = T, alternative = "less")$p.value) %>% 
      mutate(pvalue.paired_ttest=t.test(unlist(PFFno), unlist(PFFyes), paired = T)$p.value) %>% 
      mutate(pvalue.ttest_less=t.test(unlist(PFFno), unlist(PFFyes), alternative = "less")$p.value) %>% 
      mutate(pvalue.ttest=t.test(unlist(PFFno), unlist(PFFyes))$p.value) %>% 
      mutate(pvalue.paired_wilcox_less=wilcox.test(unlist(PFFno), unlist(PFFyes), paired = T, alternative = "less")$p.value) %>% 
      mutate(pvalue.paired_wilcox=wilcox.test(unlist(PFFno), unlist(PFFyes), paired = T)$p.value) %>% 
      mutate(pvalue.wilcox_less=wilcox.test(unlist(PFFno), unlist(PFFyes), alternative = "less")$p.value) %>% 
      mutate(pvalue.wilcox=wilcox.test(unlist(PFFno), unlist(PFFyes))$p.value) 
    ttest_pvalue = pivot_longer(tests_pvalue, cols = starts_with("pvalue"), names_to = "tests", names_prefix="pvalue.", values_to="pvalue") %>% filter(grepl("ttest", tests))
    
    nDecay = rownames_to_column(res_up_up) %>% 
      pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
      separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
      pivot_wider(names_from = PFF.status, values_from = log2FC) %>% 
      filter(pvalue<=1) %>%
      group_by(timepoint_comparison) %>% 
      summarise(PFFno.percent=mean(PFFno>0), PFFyes.percent=mean(PFFyes>0), PFFno.n=sum(PFFno>0), PFFyes.n=sum(PFFyes>0)) %>% 
      pivot_longer(cols = contains('PFF'), names_to = c('PFF.status',".value"), names_pattern = "(.*)\\.(.*)") 
    
    p=pivot_longer(as_tibble(res_up_up), starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
      separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
      filter(pvalue<=1) %>%
      ggplot(aes(x=log2FC, fill=PFF.status, color=PFF.status)) + 
      #geom_histogram(fill="white", alpha=0.5, position="identity", bins = 100) +
      geom_histogram(aes(y=..density..), colour="black", position="identity", alpha=0.3, bins=100)+
      geom_density(alpha=.2) +
      geom_vline(xintercept =0, color='blue', linetype=2) +
      facet_wrap(vars(timepoint_comparison), ncol = 1, nrow = 4, strip.position = "top") + 
      geom_text(data = nDecay, mapping = aes(x = Inf, y = Inf, label = paste0("Increasing genes in ",PFF.status, ": n=",n, " (",round(100*percent, 2), "%)"), color = PFF.status), hjust   = 1.05, vjust   = rep(c(1.5, 3), length(unique(nDecay$timepoint_comparison)))) +
      geom_text(data = ttest_pvalue, mapping = aes(x = Inf, y = Inf, label = paste(tests, "p = ",formatC(pvalue, format = "g", digits = 2)), color = NULL,fill= NULL), hjust   = 1.05, vjust   = rep(c(4.5, 6, 7.5, 9), length(unique(nDecay$timepoint_comparison)))) +
      ggtitle(paste0("DEresult.interactionLRT.", X ,".consecutive_", time_i, "_vs_", time_i_m_1, ".log2FC.histogram.grouped.upup.pdf")) 
    ggsave(file.path(output_dir, paste0("DEresult.interactionLRT.", X ,".consecutive_", time_i, "_vs_", time_i_m_1, ".log2FC.histogram.grouped.upup.pdf")), width = 11, height = 3)
    
    # if the increasing ones have bigger decay rate in PFFno vs. PFFyes
    p=rownames_to_column(res_up_up) %>% 
      mutate(PFFno_minus_PFFyes = (!!as.name(var_snca2)) - (!!as.name(var_snca4)))  
    p=ggplot(p, aes(x=PFFno_minus_PFFyes, fill=ifelse(PFFno_minus_PFFyes<0, paste0("PFFno<PFFyes: ",round(100*mean(PFFno_minus_PFFyes<0),2),"%"), paste0("PFFno>PFFyes: ",round(100*mean(PFFno_minus_PFFyes>0),2),"%")))) + 
      geom_histogram(colour="black", position="identity", alpha=0.3, breaks=seq(-3,3,0.1)) +
      xlab(paste(var_snca2, "-", var_snca4)) + 
      labs(fill=paste("t.test (less) p = ",formatC(t.test(p$PFFno_minus_PFFyes, alternative = "less")$p.value, format = "g", digits = 2))) + theme(legend.position="top") + 
      ggtitle(paste(time_i, "vs.", time_i_m_1, "slope difference between PFFno and PFFyes for the increasing genes (n =", nrow(res_down_down),")"))
    ggsave(file.path(output_dir, paste0("DEresult.interactionLRT.", X ,".consecutive_", time_i, "_vs_", time_i_m_1, ".log2FC.histogram.delta.upup.pdf")), width = 11, height = 3)
    
    ## now we can cluster the top significant genes by their Log2FC profiles.
    topN=1000
    topDE = filter(res, padj<=0.05) %>% head(topN) %>% as_tibble() %>% column_to_rownames('symbol') %>%  select(contains("log2FC")) %>% arrange((!!as.name(var_snca2)))
    colnames(topDE) = gsub("log2FC.","",colnames(topDE))
    dim(topDE)
    ## trim max and min
    MAX=2; MIN=-2; topDE[topDE>MAX]=MAX;topDE[topDE<MIN]=MIN
    
    ## to make sure the 0-point on the color scale is white (ref: https://stackoverflow.com/a/31707976)
    paletteLength <- 50
    myColor = colorRampPalette(c("blue", "white", "red"))(paletteLength)
    # use floor and ceiling to deal with even/odd length pallettelengths
    myBreaks = c(seq(MIN, 0, length.out=ceiling(paletteLength/2) + 1), 
                 seq(MAX/paletteLength, MAX, length.out=floor(paletteLength/2)))
    
    
    annotation_row = filter(res[1:topN,], padj<=0.05) %>% as_tibble() %>% select(geneType, symbol)  %>% column_to_rownames('symbol')
    annotation_col = data.frame(rowname=names(topDE), PFF.status=sub("\\..*","",names(topDE))) %>% column_to_rownames()
    ann_colors = list(
      #geneType = c(protein_coding = "darkblue", lincRNA = "orange", Mt_rRNA='pink', antisense='yellow', pseudogene='gray', processed_transcript='lightgray',Mt_tRNA='purple',misc_RNA='red',snoRNA='lightblue', polymorphic_pseudogene="gray", processed_pseudogene="gray", TEC="lightred"),
      PFF.status = c(PFFno = "#F8766D", PFFyes = "#00BFC4")
    )
    if(nrow(topDE)>0){
      par(cex=0.5, mar=c(5, 8, 4, 1))
      pheatmap(t(topDE),
               fontsize = 8,
               main =paste("Heatmap of log2FC for genes with LRT FDR<0.05\n",
                           file.path(output_dir, paste0("DEresult.interactionLRT.", X ,".consecutive_", time_i, "_vs_", time_i_m_1, ".padj05.xls"))),
               width=3+0.15*nrow(topDE), height=2,
               filename=file.path(output_dir, paste0("DEresult.interactionLRT.", X ,".consecutive_", time_i, "_vs_", time_i_m_1, ".topP.heatmap.pdf")),
               #border_color = NA,
               color=myColor, breaks=myBreaks,
               annotation_row = annotation_col,
               annotation_col = annotation_row,
               annotation_colors = ann_colors,
               drop_levels = TRUE,
               scale = "none", border_color = F,
               #clustering_method = 'ward.D', 
               cluster_rows = F,
               #clustering_distance_rows = "correlation",
               cluster_cols = F)
    }
    
  }
  
}
