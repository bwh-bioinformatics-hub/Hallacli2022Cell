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
# run DE among CONDITION per time point per PFF
###########################################

for(timepoint in levels(covarianceTable$Timepoint)){
  for(pff in levels(covarianceTable$PFF.status)){
    # debug: timepoint="h0"; pff="PFFno"
    message(paste("# Running DEseq on timepoint =",timepoint,"and PFF status =", pff, "..."))
    # subsetting
    dds_subset=dds[,dds$Timepoint==timepoint & dds$PFF.status==pff]; 
    dim(dds_subset);
    
    design(dds_subset) <- as.formula(" ~ CONDITION")
    
    # In versions >=1.16, the default is set to FALSE, and shrunken LFCs are obtained afterwards using lfcShrink.
    dds_subset <- DESeq(dds_subset, betaPrior=T, parallel=TRUE, BPPARAM=MulticoreParam(4))  
    resultsNames(dds_subset)
    
    comparisons=data.frame(stringsAsFactors=F);
    comparisons=rbind(comparisons, cbind("CONDITION", t(combn(levels(covarianceTable$CONDITION),2))))
    
    apply(comparisons, 1, function(x){
      # x=c("CONDITION", "SNCA2", "SNCA4")
      com_name= paste(timepoint, pff, x[3], "vs", x[2], sep="_")
      message(paste("processing comparison:",com_name))
      res <- results(dds_subset, contrast = c(x[1], x[3], x[2]),   # contract format: factor name, numerator in the fold change, denominator in the fold change
                     alpha = 0.05, 
                     parallel=TRUE, BPPARAM=MulticoreParam(4))
      ## Shrink log2 fold changes [required for version >=1.16]
      #res2 <- lfcShrink(dds, contrast = c(x[1], x[2], x[3]), res=res, type='normal', parallel=TRUE, BPPARAM=MulticoreParam(4))
      ## You can get the shrunken LFC either with lfcShrink like above or with betaPrior=TRUE. It will be the same shrunken LFC and the same as previously calculated in DESeq2. The difference is that betaPrior=TRUE will give you a p-value for the shrunken LFC, while lfcShrink (at the moment) is only giving you the LFC, and is keeping the p-value for the test of the MLE LFC. 
      ## see https://support.bioconductor.org/p/95695/ and https://support.bioconductor.org/p/98833/#98843
      
      summary(res)
      head(res); dim(res)
      # decimal value of Fold-change
      res$FoldChange <- 2**res$log2FoldChange
      # add annotation
      res$symbol <- genes_annotation$geneSymbol[match(row.names(res), genes_annotation$geneID)]
      res$geneType <- genes_annotation$geneType[match(row.names(res), genes_annotation$geneID)]
      res$geneDescription <- genes_annotation$geneDescription[match(row.names(res), genes_annotation$geneID)]
      
      # add additional columns in the output
      if(!is.null(output_additonal_columns) && grepl("m", output_additonal_columns)){ # mean of raw expression values for each group, 
        baseMeanPerLvl <- sapply( levels(dds_subset$CONDITION), function(lvl) rowMeans( counts(dds_subset,normalized=FALSE)[,dds_subset$CONDITION == lvl] ) )
        colnames(baseMeanPerLvl) = paste0("baseMean_raw.", colnames(baseMeanPerLvl))
        res = cbind(res, baseMeanPerLvl)
      }
      if(!is.null(output_additonal_columns) && grepl("M", output_additonal_columns)){ # mean of normalized expression values for each group, 
        baseMeanPerLvl <- sapply( levels(dds_subset$CONDITION), function(lvl) rowMeans( counts(dds_subset,normalized=TRUE)[,dds_subset$CONDITION == lvl] ) )
        colnames(baseMeanPerLvl) = paste0("baseMean_norm.", colnames(baseMeanPerLvl))
        res = cbind(res, baseMeanPerLvl)
      }
      if(!is.null(output_additonal_columns) && grepl("i", output_additonal_columns)){ # individual raw expression values of each sample 
        individual <- counts(dds_subset,normalized=FALSE)
        colnames(individual) = paste0("ind_raw.", colnames(individual))
        res = cbind(res, individual)
      }
      if(!is.null(output_additonal_columns) && grepl("I", output_additonal_columns)){ # individual normalized expression values of each sample
        individual <- counts(dds_subset,normalized=TRUE)
        colnames(individual) = paste0("ind_norm.", colnames(individual))
        res = cbind(res, individual)
      }
      res <- res[order(res$padj),]
      head(res); dim(res)
      
      write.table(as.data.frame(res), 
                  file=file.path(output_dir, paste0("DEresult.all.", com_name ,".xls.gz")), 
                  sep="\t", quote =F, na="", row.names=T, col.names = NA)
      
      write.table(as.data.frame(subset(res, padj<=0.05 & abs(log2FoldChange)>=1)), 
                  file=file.path(output_dir, paste0("DEresult.padj05_log2FCgt1.", com_name ,".xls")), 
                  sep="\t", quote =F, na="", row.names=T, col.names = NA)
      
      # ## Note: 10% of non-significant(NS) genes are randomly selected in order to increase the performance of the generating graph
      NS=subset(res, padj>0.05 | abs(log2FoldChange)<1)
      n_NS=nrow(NS)
      NS=NS[sample(n_NS,round(n_NS*.20)),]
      dim(res); res=DESeqResults(rbind(NS, subset(res, padj<=0.05 & abs(log2FoldChange)>=1))); dim(res);
      
      ## MAKING PLOTS
      pdf(file.path(output_dir, paste0("DEresult.padj_05.", com_name ,".pdf")), paper = 'USr')
      # scatter plot
      # TOADD: 
      
      ## ==============================
      # MA plot
      ## ==============================
      DESeq2::plotMA(res, alpha = 0.05, colNonSig = "gray", main=paste0(output_dir,": ",com_name))
      
      ## ==============================
      # vocano plot
      ## ==============================
      topT <- as.data.frame(res)
      with(topT, plot(log2FoldChange, -log10(padj), 
                      pch=20, cex=0.5, main=paste0(output_dir,": ",com_name), col='gray',
                      xlab=bquote(~Log[2]~fold~change), 
                      ylab=bquote(~-log[10]~FDR)))
      if(nrow(subset(topT, padj<=0.05 & abs(log2FoldChange)>=1))>0){
        with(subset(topT, padj<=0.05 & log2FoldChange>=1), 
             points(log2FoldChange, -log10(padj), pch=20, col="red", cex=1))
        with(subset(topT, padj<=0.05 & log2FoldChange<=-1), 
             points(log2FoldChange, -log10(padj), pch=20, col="blue", cex=1))
        with(subset(topT, padj<=0.05 & abs(log2FoldChange)>=1), 
             text(log2FoldChange, -log10(padj), labels=symbol, col="gray", cex=0.5, pos=1, offset=0.2))
      }
      abline(v=c(0,-1, 1), lty=c(3,4,4), lwd=c(1,2,2))
      abline(h=-log10(0.05), col="black", lty=4, lwd=2.0)
      
      ## enhanced Vocano plot
      library(EnhancedVolcano)
      EnhancedVolcano(res,
                      lab = res$symbol,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      drawConnectors = TRUE,
                      arrowheads =F)
      
      if(nrow(subset(res, abs(log2FoldChange)>=1 & padj<0.05))>0){    
        ## ==============================
        message("# heatmap for top 10 DE genes")
        ## ==============================
        topT=subset(res, abs(log2FoldChange)>=1 & padj<0.05)
        
        #topT=topT[order(-topT$baseMean),] # sort by baseMean in decreasing order
        #topT=topT[order(-abs(topT$log2FoldChange)),] # sort by abs(log2FoldChange) in decreasing order
        topT=topT[order(topT$padj),] # sort by padj in increasing order
        topT=rbind(head(subset(topT, log2FoldChange<0),10),head(subset(topT, log2FoldChange>0),10)) # top 10 down-regulated and top 10 up-regulated
        topDE=assay(vsd)[rownames(topT),]
        rownames(topDE) = topT$symbol
        annotation_row = dplyr::select(as.data.frame(topT), geneType, log2FoldChange, symbol) %>% mutate(updown=ifelse(log2FoldChange>0,"up","down")) %>% select(updown, geneType)
        rownames(annotation_row) = topT$symbol
        annotation_col = dplyr::select(as.data.frame(colData(dds)), CONDITION)
        ann_colors = list(
          updown = c(up = "red", down = "blue"),
          #geneType = c(protein_coding = "darkblue", lincRNA = "orange", Mt_rRNA='pink', antisense='yellow', pseudogene='gray', processed_transcript='lightgray',Mt_tRNA='purple',misc_RNA='red',snoRNA='lightblue', polymorphic_pseudogene="gray", processed_pseudogene="gray", TEC="lightred"),
          CONDITION = c(SNCA2 = "green", SNCA4 = "red")
        )
        
        ## Scale/center each genes (by rows)
        topDE=t(scale(t(as.matrix(topDE))))
        ## trim max and min
        MAX=2; MIN=-2; topDE[topDE>MAX]=MAX;topDE[topDE<MIN]=MIN
        
        par(cex=0.5, mar=c(5, 8, 4, 1))
        pheatmap(topDE,
                 fontsize = 8,
                 main =paste0(output_dir,": heatmap for top 10 DE genes in ", com_name),
                 #fontsize_col = 5,
                 #width=8, height=5,filename=paste0("DEresult.padj_05.", com_name ,".heatmap.pdf"),
                 border_color = NA,
                 color = colorRampPalette(c("blue", "white", "red"))(50),
                 annotation_row = annotation_row,
                 annotation_col = annotation_col,
                 annotation_colors = ann_colors,
                 drop_levels = TRUE,
                 scale = "none", 
                 clustering_method = 'ward.D', 
                 cluster_rows = TRUE,
                 clustering_distance_rows = "correlation",
                 cutree_rows = 2,cutree_cols=2,
                 cluster_cols = TRUE,
                 clustering_distance_cols = "correlation")
        
        ## ==============================
        message("# heatmap for top 20 DE genes")
        ## ==============================
        topT=subset(res, abs(log2FoldChange)>=1 & padj<0.05)
        topT=topT[order(topT$padj),] # sort by padj in increasing order
        topT=rbind(head(subset(topT, log2FoldChange<0),20),head(subset(topT, log2FoldChange>0),20)) # top 10 down-regulated and top 10 up-regulated
        topDE=assay(vsd)[rownames(topT),]; rownames(topDE) = topT$symbol
        annotation_row = dplyr::select(as.data.frame(topT), geneType, log2FoldChange, symbol) %>% mutate(updown=ifelse(log2FoldChange>0,"up","down")) %>% select(updown, geneType)
        rownames(annotation_row) = topT$symbol
        annotation_col = dplyr::select(as.data.frame(colData(dds)), CONDITION)
        ann_colors = list(
          updown = c(up = "red", down = "blue"),
          #geneType = c(protein_coding = "darkblue", lincRNA = "orange", Mt_rRNA='pink', antisense='yellow', pseudogene='gray', processed_transcript='lightgray',Mt_tRNA='purple',misc_RNA='red',snoRNA='lightblue'),
          CONDITION = c(SNCA2 = "green", SNCA4 = "red")
        )
        
        ## Scale/center each genes (by rows)
        topDE=t(scale(t(as.matrix(topDE))))
        ## trim max and min
        MAX=2; MIN=-2; topDE[topDE>MAX]=MAX;topDE[topDE<MIN]=MIN
        
        par(cex=0.5, mar=c(5, 8, 4, 1))
        pheatmap(topDE,
                 main =paste0(output_dir,": heatmap for top 20 DE genes in ", com_name),
                 fontsize = 8,
                 #fontsize_col = 5,
                 border_color = NA,
                 color = colorRampPalette(c("blue", "white", "red"))(50),
                 annotation_row = annotation_row,
                 annotation_col = annotation_col,
                 annotation_colors = ann_colors,
                 drop_levels = TRUE,
                 scale = "none", 
                 clustering_method = 'ward.D', 
                 cluster_rows = TRUE,
                 clustering_distance_rows = "correlation",
                 cutree_rows = 2,cutree_cols=2,
                 cluster_cols = TRUE,
                 clustering_distance_cols = "correlation")
        
      }
      
      dev.off()     
      
      # # Exporting results to HTML
      # htmlRep <- HTMLReport(shortName=com_name, title=com_name, reportDirectory="./report")
      # publish(makeNewImages(head(as.data.frame(subset(res, !is.na(padj) & pvalue<0.05)),1000)), htmlRep)
      # finish(htmlRep)
    })
  }
}

###########################################
# run DE among PFF per time point per CONDITION
###########################################

for(timepoint in levels(covarianceTable$Timepoint)){
  for(condition in levels(covarianceTable$CONDITION)){
    # debug: timepoint="h0"; condition="SNCA2"
    message(paste("# Running DEseq on timepoint =",timepoint,"and CONDITION status =", condition, "..."))
    # subsetting
    dds_subset=dds[,dds$Timepoint==timepoint & dds$CONDITION==condition]; 
    dim(dds_subset);
    
    design(dds_subset) <- as.formula(" ~ PFF.status")
    
    # In versions >=1.16, the default is set to FALSE, and shrunken LFCs are obtained afterwards using lfcShrink.
    dds_subset <- DESeq(dds_subset, betaPrior=T, parallel=TRUE, BPPARAM=MulticoreParam(4))  
    resultsNames(dds_subset)
    
    comparisons=data.frame(stringsAsFactors=F);
    comparisons=rbind(comparisons, cbind("PFF.status", t(combn(levels(covarianceTable$PFF.status),2))))
    
    apply(comparisons, 1, function(x){
      # x=c("PFF.status", "PFFyes", "PFFno")
      com_name= paste(timepoint, condition, x[3], "vs", x[2], sep="_")
      message(paste("processing comparison:",com_name))
      res <- results(dds_subset, contrast = c(x[1], x[3], x[2]),   # contract format: factor name, numerator in the fold change, denominator in the fold change
                     alpha = 0.05, 
                     parallel=TRUE, BPPARAM=MulticoreParam(4))
      ## Shrink log2 fold changes [required for version >=1.16]
      #res2 <- lfcShrink(dds, contrast = c(x[1], x[2], x[3]), res=res, type='normal', parallel=TRUE, BPPARAM=MulticoreParam(4))
      ## You can get the shrunken LFC either with lfcShrink like above or with betaPrior=TRUE. It will be the same shrunken LFC and the same as previously calculated in DESeq2. The difference is that betaPrior=TRUE will give you a p-value for the shrunken LFC, while lfcShrink (at the moment) is only giving you the LFC, and is keeping the p-value for the test of the MLE LFC. 
      ## see https://support.bioconductor.org/p/95695/ and https://support.bioconductor.org/p/98833/#98843
      
      summary(res)
      head(res); dim(res)
      # decimal value of Fold-change
      res$FoldChange <- 2**res$log2FoldChange
      # add annotation
      res$symbol <- genes_annotation$geneSymbol[match(row.names(res), genes_annotation$geneID)]
      res$geneType <- genes_annotation$geneType[match(row.names(res), genes_annotation$geneID)]
      res$geneDescription <- genes_annotation$geneDescription[match(row.names(res), genes_annotation$geneID)]
      
      # add additional columns in the output
      if(!is.null(output_additonal_columns) && grepl("m", output_additonal_columns)){ # mean of raw expression values for each group, 
        baseMeanPerLvl <- sapply( levels(dds_subset$PFF.status), function(lvl) rowMeans( counts(dds_subset,normalized=FALSE)[,dds_subset$PFF.status == lvl] ) )
        colnames(baseMeanPerLvl) = paste0("baseMean_raw.", colnames(baseMeanPerLvl))
        res = cbind(res, baseMeanPerLvl)
      }
      if(!is.null(output_additonal_columns) && grepl("M", output_additonal_columns)){ # mean of normalized expression values for each group, 
        baseMeanPerLvl <- sapply( levels(dds_subset$PFF.status), function(lvl) rowMeans( counts(dds_subset,normalized=TRUE)[,dds_subset$PFF.status == lvl] ) )
        colnames(baseMeanPerLvl) = paste0("baseMean_norm.", colnames(baseMeanPerLvl))
        res = cbind(res, baseMeanPerLvl)
      }
      if(!is.null(output_additonal_columns) && grepl("i", output_additonal_columns)){ # individual raw expression values of each sample 
        individual <- counts(dds_subset,normalized=FALSE)
        colnames(individual) = paste0("ind_raw.", colnames(individual))
        res = cbind(res, individual)
      }
      if(!is.null(output_additonal_columns) && grepl("I", output_additonal_columns)){ # individual normalized expression values of each sample
        individual <- counts(dds_subset,normalized=TRUE)
        colnames(individual) = paste0("ind_norm.", colnames(individual))
        res = cbind(res, individual)
      }
      res <- res[order(res$padj),]
      head(res); dim(res)
      
      write.table(as.data.frame(res), 
                  file=file.path(output_dir, paste0("DEresult.all.", com_name ,".xls")), 
                  sep="\t", quote =F, na="", row.names=T, col.names = NA)
      
      write.table(as.data.frame(subset(res, padj<=0.05 & abs(log2FoldChange)>=1)), 
                  file=file.path(output_dir, paste0("DEresult.padj05_log2FCgt1.", com_name ,".xls")), 
                  sep="\t", quote =F, na="", row.names=T, col.names = NA)
      
      # ## Note: 10% of non-significant(NS) genes are randomly selected in order to increase the performance of the generating graph
      NS=subset(res, padj>0.05 | abs(log2FoldChange)<1)
      n_NS=nrow(NS)
      NS=NS[sample(n_NS,round(n_NS*.20)),]
      dim(res); res=DESeqResults(rbind(NS, subset(res, padj<=0.05 & abs(log2FoldChange)>=1))); dim(res);
      
      ## MAKING PLOTS
      pdf(file.path(output_dir, paste0("DEresult.padj_05.", com_name ,".pdf")), paper = 'USr')
      # scatter plot
      # TOADD: 
      
      ## ==============================
      # MA plot
      ## ==============================
      DESeq2::plotMA(res, alpha = 0.05, colNonSig = "gray", main=paste0(output_dir,": ",com_name))
      
      ## ==============================
      # vocano plot
      ## ==============================
      topT <- as.data.frame(res)
      with(topT, plot(log2FoldChange, -log10(padj), 
                      pch=20, cex=0.5, main=paste0(output_dir,": ",com_name), col='gray',
                      xlab=bquote(~Log[2]~fold~change), 
                      ylab=bquote(~-log[10]~FDR)))
      if(nrow(subset(topT, padj<=0.05 & abs(log2FoldChange)>=1))>0){
        with(subset(topT, padj<=0.05 & log2FoldChange>=1), 
             points(log2FoldChange, -log10(padj), pch=20, col="red", cex=1))
        with(subset(topT, padj<=0.05 & log2FoldChange<=-1), 
             points(log2FoldChange, -log10(padj), pch=20, col="blue", cex=1))
        with(subset(topT, padj<=0.05 & abs(log2FoldChange)>=1), 
             text(log2FoldChange, -log10(padj), labels=symbol, col="gray", cex=0.5, pos=1, offset=0.2))
      }
      abline(v=c(0,-1, 1), lty=c(3,4,4), lwd=c(1,2,2))
      abline(h=-log10(0.05), col="black", lty=4, lwd=2.0)
      
      ## enhanced Vocano plot
      library(EnhancedVolcano)
      p=EnhancedVolcano(res,
                        lab = res$symbol,
                        x = 'log2FoldChange',
                        y = 'pvalue',
                        drawConnectors = TRUE,
                        arrowheads =F)
      print(p)
      
      if(nrow(subset(res, abs(log2FoldChange)>=1 & padj<0.05))>0){    
        ## ==============================
        message("# heatmap for top 10 DE genes")
        ## ==============================
        topT=subset(res, abs(log2FoldChange)>=1 & padj<0.05)
        
        #topT=topT[order(-topT$baseMean),] # sort by baseMean in decreasing order
        #topT=topT[order(-abs(topT$log2FoldChange)),] # sort by abs(log2FoldChange) in decreasing order
        topT=topT[order(topT$padj),] # sort by padj in increasing order
        topT=rbind(head(subset(topT, log2FoldChange<0),10),head(subset(topT, log2FoldChange>0),10)) # top 10 down-regulated and top 10 up-regulated
        topDE=assay(vsd)[rownames(topT),]
        rownames(topDE) = topT$symbol
        annotation_row = dplyr::select(as.data.frame(topT), geneType, log2FoldChange, symbol) %>% mutate(updown=ifelse(log2FoldChange>0,"up","down")) %>% select(updown, geneType)
        rownames(annotation_row) = topT$symbol
        annotation_col = dplyr::select(as.data.frame(colData(dds)), PFF.status)
        ann_colors = list(
          updown = c(up = "red", down = "blue"),
          #geneType = c(protein_coding = "darkblue", lincRNA = "orange", Mt_rRNA='pink', antisense='yellow', pseudogene='gray', processed_transcript='lightgray',Mt_tRNA='purple',misc_RNA='red',snoRNA='lightblue', polymorphic_pseudogene="gray", processed_pseudogene="gray", TEC="lightred"),
          PFF.status = c(PFFno = "green", PFFyes = "red")
        )
        
        ## Scale/center each genes (by rows)
        topDE=t(scale(t(as.matrix(topDE))))
        ## trim max and min
        MAX=2; MIN=-2; topDE[topDE>MAX]=MAX;topDE[topDE<MIN]=MIN
        
        par(cex=0.5, mar=c(5, 8, 4, 1))
        pheatmap(topDE,
                 fontsize = 8,
                 main =paste0(output_dir,": heatmap for top 10 DE genes in ", com_name),
                 #fontsize_col = 5,
                 #width=8, height=5,filename=paste0("DEresult.padj_05.", com_name ,".heatmap.pdf"),
                 border_color = NA,
                 color = colorRampPalette(c("blue", "white", "red"))(50),
                 annotation_row = annotation_row,
                 annotation_col = annotation_col,
                 annotation_colors = ann_colors,
                 drop_levels = TRUE,
                 scale = "none", 
                 clustering_method = 'ward.D', 
                 cluster_rows = TRUE,
                 clustering_distance_rows = "correlation",
                 cutree_rows = 2,cutree_cols=2,
                 cluster_cols = TRUE,
                 clustering_distance_cols = "correlation")
        
        ## ==============================
        message("# heatmap for top 20 DE genes")
        ## ==============================
        topT=subset(res, abs(log2FoldChange)>=1 & padj<0.05)
        topT=topT[order(topT$padj),] # sort by padj in increasing order
        topT=rbind(head(subset(topT, log2FoldChange<0),20),head(subset(topT, log2FoldChange>0),20)) # top 10 down-regulated and top 10 up-regulated
        topDE=assay(vsd)[rownames(topT),]; rownames(topDE) = topT$symbol
        annotation_row = dplyr::select(as.data.frame(topT), geneType, log2FoldChange, symbol) %>% mutate(updown=ifelse(log2FoldChange>0,"up","down")) %>% select(updown, geneType)
        rownames(annotation_row) = topT$symbol
        annotation_col = dplyr::select(as.data.frame(colData(dds)), PFF.status)
        ann_colors = list(
          updown = c(up = "red", down = "blue"),
          #geneType = c(protein_coding = "darkblue", lincRNA = "orange", Mt_rRNA='pink', antisense='yellow', pseudogene='gray', processed_transcript='lightgray',Mt_tRNA='purple',misc_RNA='red',snoRNA='lightblue'),
          PFF.status = c(PFFno = "green", PFFyes = "red")
        )
        
        ## Scale/center each genes (by rows)
        topDE=t(scale(t(as.matrix(topDE))))
        ## trim max and min
        MAX=2; MIN=-2; topDE[topDE>MAX]=MAX;topDE[topDE<MIN]=MIN
        
        par(cex=0.5, mar=c(5, 8, 4, 1))
        pheatmap(topDE,
                 main =paste0(output_dir,": heatmap for top 20 DE genes in ", com_name),
                 fontsize = 8,
                 #fontsize_col = 5,
                 border_color = NA,
                 color = colorRampPalette(c("blue", "white", "red"))(50),
                 annotation_row = annotation_row,
                 annotation_col = annotation_col,
                 annotation_colors = ann_colors,
                 drop_levels = TRUE,
                 scale = "none", 
                 clustering_method = 'ward.D', 
                 cluster_rows = TRUE,
                 clustering_distance_rows = "correlation",
                 cutree_rows = 2,cutree_cols=2,
                 cluster_cols = TRUE,
                 clustering_distance_cols = "correlation")
        
      }
      
      dev.off()     
      
      # # Exporting results to HTML
      # htmlRep <- HTMLReport(shortName=com_name, title=com_name, reportDirectory="./report")
      # publish(makeNewImages(head(as.data.frame(subset(res, !is.na(padj) & pvalue<0.05)),1000)), htmlRep)
      # finish(htmlRep)
    })
  }
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
  
  # ## Wald tests for the log2 fold changes at individual time points
  # res_h1 = results(dds_subset, name="CONDITIONSNCA4.Timepointh1", test="Wald")
  # res_h1$symbol <- genes_annotation$geneSymbol[match(row.names(res_h1), genes_annotation$geneID)]
  # res_h1$geneType <- genes_annotation$geneType[match(row.names(res_h1), genes_annotation$geneID)]
  # res_h1 <- res_h1[order(res_h1$padj),]
  # head(res_h1)
  
  # extract a matrix of the log2 fold changes for all comparisons using the coef function.
  betas <- coef(dds_subset)
  colnames(betas)  # same as the resultsNames(dds_subset)
  betas = betas[rownames(res), -c(1,2)]
  ## Note that CONDITIONSNCA4.Timepointh1 means difference between SNCA4 vs SNCA2 at h1, controlling for baseline. 
  ## To get the log2 fold change of h1 vs h0 for the SNCA2, Timepoint_h1_vs_h0
  ## To get the log2 fold change of h1 vs h0 for the SNCA4, Timepoint_h1_vs_h0 +  CONDITIONSNCA4.Timepointh1
  betas = as.data.frame(betas) %>% mutate(.keep = "none",
                                          SNCA2.h1_vs_h0  = Timepoint_h1_vs_h0, 
                                          SNCA2.h3_vs_h0  = Timepoint_h3_vs_h0, 
                                          SNCA2.h6_vs_h0  = Timepoint_h6_vs_h0, 
                                          SNCA2.h12_vs_h0 = Timepoint_h12_vs_h0, 
                                          SNCA4.h1_vs_h0  = Timepoint_h1_vs_h0 + CONDITIONSNCA4.Timepointh1,
                                          SNCA4.h3_vs_h0  = Timepoint_h3_vs_h0 + CONDITIONSNCA4.Timepointh3,
                                          SNCA4.h6_vs_h0  = Timepoint_h6_vs_h0 + CONDITIONSNCA4.Timepointh6,
                                          SNCA4.h12_vs_h0 = Timepoint_h12_vs_h0 + CONDITIONSNCA4.Timepointh12)
  
  colnames(betas) = paste0("log2FC.",colnames(betas)) 
  
  # combine the overall LRT tesst (for any timepoint) and the log2FC for each timepoint together
  res = cbind(as.data.frame(res)[,c('baseMean', 'stat','pvalue','padj')], betas)
  
  # add annotation
  res$symbol <- genes_annotation$geneSymbol[match(row.names(res), genes_annotation$geneID)]
  res$geneType <- genes_annotation$geneType[match(row.names(res), genes_annotation$geneID)]
  res$geneDescription <- genes_annotation$geneDescription[match(row.names(res), genes_annotation$geneID)]
  
  ## save the data
  res <- res[order(res$padj),]
  head(res); dim(res)
  
  write.table(as.data.frame(res), 
              file=gzfile(file.path(output_dir, paste0("DEresult.interactionLRT.", pff ,".all.xls.gz"))), 
              sep="\t", quote =F, na="", row.names=T, col.names = NA)
  
  write.table(subset(res, padj<=0.05), 
              file=file.path(output_dir, paste0("DEresult.interactionLRT.", pff ,".padj05.xls")), 
              sep="\t", quote =F, na="", row.names=T, col.names = NA)
  
  ## any globe difference in slope between SNCA2 vs. SNCA4
  ttest_pvalue = rownames_to_column(res) %>% 
    pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
    pivot_wider(names_from = CONDITION, values_from = log2FC) %>% 
    filter(pvalue<=1) %>%
    mutate(timepoint_comparison=factor(timepoint_comparison, levels = paste0("h",c(1,3,6,12),"_vs_h0"))) %>%
    group_by(timepoint_comparison) %>% 
    summarise(SNCA2=list(SNCA2), SNCA4=list(SNCA4)) %>% rowwise() %>% 
    mutate(pvalue.paired_ttest_less=t.test(unlist(SNCA2), unlist(SNCA4), paired = T, alternative = "less")$p.value) %>% 
    mutate(pvalue.paired_ttest=t.test(unlist(SNCA2), unlist(SNCA4), paired = T)$p.value) %>% 
    mutate(pvalue.ttest_less=t.test(unlist(SNCA2), unlist(SNCA4), alternative = "less")$p.value) %>% 
    mutate(pvalue.ttest=t.test(unlist(SNCA2), unlist(SNCA4))$p.value) %>% 
    mutate(pvalue.paired_wilcox_less=wilcox.test(unlist(SNCA2), unlist(SNCA4), paired = T, alternative = "less")$p.value) %>% 
    mutate(pvalue.paired_wilcox=wilcox.test(unlist(SNCA2), unlist(SNCA4), paired = T)$p.value) %>% 
    mutate(pvalue.wilcox_less=wilcox.test(unlist(SNCA2), unlist(SNCA4), alternative = "less")$p.value) %>% 
    mutate(pvalue.wilcox=wilcox.test(unlist(SNCA2), unlist(SNCA4))$p.value) %>% 
    pivot_longer(cols = starts_with("pvalue"), names_to = "tests", names_prefix="pvalue.", values_to="pvalue")
  
  # # pff=yes
  # # Rowwise: 
  # timepoint_comparison SNCA2          SNCA4          pvalue.paired_ttest_less pvalue.paired_ttest pvalue.ttest_less pvalue.ttest pvalue.paired_wilcox_less pvalue.paired_wilcox pvalue.wilcox_less pvalue.wilcox
  # <fct>                <list>         <list>                            <dbl>               <dbl>             <dbl>        <dbl>                     <dbl>                <dbl>              <dbl>         <dbl>
  # 1 h1_vs_h0             <dbl [26,219]> <dbl [26,219]>                 0.0200              0.0400              0.0402       0.0804                     0.202        0.404                     0.190         0.379 
  # 2 h3_vs_h0             <dbl [26,219]> <dbl [26,219]>                 0.999               0.00129             0.991        0.0181                     1.00         0.00000000216             0.825         0.350 
  # 3 h6_vs_h0             <dbl [26,219]> <dbl [26,219]>                 0.242               0.484               0.320        0.640                      0.650        0.700                     0.156         0.312 
  # 4 h12_vs_h0            <dbl [26,219]> <dbl [26,219]>                 0.000102            0.000205            0.0105       0.0209                     0.180        0.360                     0.0220        0.0440
  # 
  
  # # pff=no
  # # Rowwise: 
  # timepoint_comparison SNCA2          SNCA4          pvalue.paired_ttest_less pvalue.paired_ttest pvalue.ttest_less pvalue.ttest pvalue.paired_wilcox_less pvalue.paired_wilcox pvalue.wilcox_less pvalue.wilcox
  # <fct>                <list>         <list>                            <dbl>               <dbl>             <dbl>        <dbl>                     <dbl>                <dbl>              <dbl>         <dbl>
  # 1 h1_vs_h0             <dbl [26,218]> <dbl [26,218]>       0.966               0.0677                  0.944        0.113               0.999                0.00114                   0.987         0.0257     
  # 2 h3_vs_h0             <dbl [26,218]> <dbl [26,218]>       0.0279              0.0558                  0.0754       0.151               1.00                 0.000000621               0.227         0.454      
  # 3 h6_vs_h0             <dbl [26,218]> <dbl [26,218]>       0.0000000000000117  0.0000000000000233      0.0000000350 0.0000000699        0.0000000000000111   0.0000000000000222        0.000000184   0.000000369
  # 4 h12_vs_h0            <dbl [26,218]> <dbl [26,218]>       1.00                0.000581                0.985        0.0307              1.00                 0.000000193               0.806         0.388       
  # 
  nDecay = rownames_to_column(res) %>% 
    pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
    pivot_wider(names_from = CONDITION, values_from = log2FC) %>% 
    filter(pvalue<=1) %>%
    group_by(timepoint_comparison) %>% 
    summarise(SNCA2.percent=mean(SNCA2<0), SNCA4.percent=mean(SNCA4<0), SNCA2.n=sum(SNCA2<0), SNCA4.n=sum(SNCA4<0)) %>% 
    pivot_longer(cols = contains('SNCA'), names_to = c('CONDITION',".value"), names_pattern = "(.*)\\.(.*)") %>% 
    mutate(timepoint_comparison=factor(timepoint_comparison, levels = paste0("h",c(1,3,6,12),"_vs_h0"))) 
  
  pivot_longer(as_tibble(res), starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
    mutate(timepoint_comparison=factor(timepoint_comparison, levels = paste0("h",c(1,3,6,12),"_vs_h0"))) %>%
    filter(pvalue<=1) %>%
    ggplot(aes(x=log2FC, fill=CONDITION, color=CONDITION)) + 
    #geom_histogram(fill="white", alpha=0.5, position="identity", bins = 100) +
    geom_histogram(aes(y=..density..), colour="black", position="identity", alpha=0.3, bins=100)+
    geom_density(alpha=.2) +
    geom_vline(xintercept =0, color='blue', linetype=2) +
    facet_wrap(vars(timepoint_comparison), ncol = 1, nrow = 4, strip.position = "top") + 
    geom_text(data = filter(ttest_pvalue,grepl("ttest", tests)), mapping = aes(x = Inf, y = Inf, label = paste(tests, "p = ",formatC(pvalue, format = "g", digits = 2)), color = NULL,fill= NULL), hjust   = 1.05, vjust   = rep(c(1.5, 3, 4.5, 6),4)) +
    geom_text(data = nDecay, mapping = aes(x = -Inf, y = Inf, label = paste0("Decreasing genes in ",CONDITION, ": n=",n, " (",round(100*percent, 2), "%)"), color = CONDITION), hjust   = -0.05, vjust   = rep(c(1.5, 3),4)) +
    ggtitle(paste0("DEresult.decayDiff.histogram.grouped.", pff ,".all.pdf")) 
  ggsave(file.path(output_dir, paste0("DEresult.decayDiff.histogram.grouped.", pff ,".all.pdf")), width = 11, height = 10)
  
  # histogram for negative slope genes @ h1 (same as the Model1 by Erinc's reference)
  res_negative_h1 = filter(res, log2FC.SNCA2.h1_vs_h0<0, log2FC.SNCA4.h1_vs_h0<0)
  ttest_pvalue = rownames_to_column(res_negative_h1) %>% 
    pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
    pivot_wider(names_from = CONDITION, values_from = log2FC) %>% 
    filter(pvalue<=1) %>%
    mutate(timepoint_comparison=factor(timepoint_comparison, levels = paste0("h",c(1,3,6,12),"_vs_h0"))) %>%
    group_by(timepoint_comparison) %>% 
    summarise(SNCA2=list(SNCA2), SNCA4=list(SNCA4)) %>% rowwise() %>% 
    mutate(pvalue.paired_ttest_less=t.test(unlist(SNCA2), unlist(SNCA4), paired = T, alternative = "less")$p.value) %>% 
    mutate(pvalue.paired_ttest=t.test(unlist(SNCA2), unlist(SNCA4), paired = T)$p.value) %>% 
    mutate(pvalue.ttest_less=t.test(unlist(SNCA2), unlist(SNCA4), alternative = "less")$p.value) %>% 
    mutate(pvalue.ttest=t.test(unlist(SNCA2), unlist(SNCA4))$p.value) %>% 
    mutate(pvalue.paired_wilcox_less=wilcox.test(unlist(SNCA2), unlist(SNCA4), paired = T, alternative = "less")$p.value) %>% 
    mutate(pvalue.paired_wilcox=wilcox.test(unlist(SNCA2), unlist(SNCA4), paired = T)$p.value) %>% 
    mutate(pvalue.wilcox_less=wilcox.test(unlist(SNCA2), unlist(SNCA4), alternative = "less")$p.value) %>% 
    mutate(pvalue.wilcox=wilcox.test(unlist(SNCA2), unlist(SNCA4))$p.value) %>% 
    pivot_longer(cols = starts_with("pvalue"), names_to = "tests", names_prefix="pvalue.", values_to="pvalue")
  
  nDecay = rownames_to_column(res_negative_h1) %>% 
    pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
    pivot_wider(names_from = CONDITION, values_from = log2FC) %>% 
    filter(pvalue<=1) %>%
    group_by(timepoint_comparison) %>% 
    summarise(SNCA2.percent=mean(SNCA2<0), SNCA4.percent=mean(SNCA4<0), SNCA2.n=sum(SNCA2<0), SNCA4.n=sum(SNCA4<0)) %>% 
    pivot_longer(cols = contains('SNCA'), names_to = c('CONDITION',".value"), names_pattern = "(.*)\\.(.*)") %>% 
    mutate(timepoint_comparison=factor(timepoint_comparison, levels = paste0("h",c(1,3,6,12),"_vs_h0"))) 
  
  pivot_longer(as_tibble(res_negative_h1), starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
    mutate(timepoint_comparison=factor(timepoint_comparison, levels = paste0("h",c(1,3,6,12),"_vs_h0"))) %>%
    filter(pvalue<=1) %>%
    ggplot(aes(x=log2FC, fill=CONDITION, color=CONDITION)) + 
    #geom_histogram(fill="white", alpha=0.5, position="identity", bins = 100) +
    geom_histogram(aes(y=..density..), colour="black", position="identity", alpha=0.3, bins=100)+
    geom_density(alpha=.2) +
    geom_vline(xintercept =0, color='blue', linetype=2) +
    facet_wrap(vars(timepoint_comparison), ncol = 1, nrow = 4, strip.position = "top") + 
    geom_text(data = filter(ttest_pvalue,grepl("ttest", tests)), mapping = aes(x = Inf, y = Inf, label = paste(tests, "p = ",formatC(pvalue, format = "g", digits = 2)), color = NULL,fill= NULL), hjust   = 1.05, vjust   = rep(c(1.5, 3, 4.5, 6),4)) +
    geom_text(data = nDecay, mapping = aes(x = -Inf, y = Inf, label = paste0("Decreasing genes in ",CONDITION, ": n=",n, " (",round(100*percent, 2), "%)"), color = CONDITION), hjust   = -0.05, vjust   = rep(c(1.5, 3),4)) +
    ggtitle(paste0("DEresult.decayDiff.histogram.grouped.", pff ,".res_negative_h1.pdf")) 
  ggsave(file.path(output_dir, paste0("DEresult.decayDiff.histogram.grouped.", pff ,".subset_negative_h1.pdf")), width = 11, height = 10)
  
  # if the decreasing ones have bigger decay rate in SNCA2 vs. SNCA4
  rownames_to_column(res_negative_h1) %>% 
    mutate(SNCA2_minus_SNCA4 = log2FC.SNCA2.h1_vs_h0-log2FC.SNCA4.h1_vs_h0) %>% 
    ggplot(aes(x=SNCA2_minus_SNCA4, fill=ifelse(SNCA2_minus_SNCA4<0, paste0("SNCA2<SNCA4: ",round(100*mean(SNCA2_minus_SNCA4<0),2),"%"), paste0("SNCA2>SNCA4: ",round(100*mean(SNCA2_minus_SNCA4>0),2),"%")))) + 
    geom_histogram(colour="black", position="identity", alpha=0.3, breaks=seq(-3,3,0.1)) + 
    xlab("log2FC.SNCA2.h1_vs_h0 - log2FC.SNCA4.h1_vs_h0") + 
    labs(fill="") + theme(legend.position="top") + 
    ggtitle(paste("Slope difference between SNCA2 vs. SNCA4 for the genes with both negative slopes at h1 (n =", nrow(res_negative_h1),")")) 
  ggsave(file.path(output_dir, paste0("DEresult.decayDiff.histogram.grouped.", pff ,".subset_negative_h1.h1diff.pdf")), width = 10, height = 4)
  
  ## monotonic ones (i.e. keep decreasing)
  # all(x == cummin(x)) # see https://stackoverflow.com/a/13094801
  monotonic_decreasing <- function(x, na.rm = T) all(x == cummin(x))
  df=rownames_to_column(res) %>% rowwise() %>% 
    filter(max(log2FC.SNCA2.h1_vs_h0, log2FC.SNCA2.h3_vs_h0, log2FC.SNCA2.h6_vs_h0, log2FC.SNCA2.h12_vs_h0)<0,
           max(log2FC.SNCA4.h1_vs_h0, log2FC.SNCA4.h3_vs_h0, log2FC.SNCA4.h6_vs_h0, log2FC.SNCA4.h12_vs_h0)<0) %>% 
    mutate(monotonic_SNCA2 = monotonic_decreasing(c(log2FC.SNCA2.h1_vs_h0, log2FC.SNCA2.h3_vs_h0, log2FC.SNCA2.h6_vs_h0, log2FC.SNCA2.h12_vs_h0)),
           monotonic_SNCA4 = monotonic_decreasing(c(log2FC.SNCA4.h1_vs_h0, log2FC.SNCA4.h3_vs_h0, log2FC.SNCA4.h6_vs_h0, log2FC.SNCA4.h12_vs_h0))) %>% 
    filter(monotonic_SNCA2 | monotonic_SNCA4) %>% print(n_extra=17) 
  with(df, table(monotonic_SNCA2, monotonic_SNCA4))
  write.table(df, file=file.path(output_dir, paste0("DEresult.interactionLRT.", pff ,".monotonic.xls")), 
              sep="\t", quote =F, na="", row.names = F)
  
  # if the monotonic ones have bigger decay rate in SNCA2 vs. SNCA4
  df2= df %>% filter(monotonic_SNCA2, monotonic_SNCA4) %>%  #dim()
    #mutate(SNCA2_lt_SNCA4 = log2FC.SNCA2.h1_vs_h0 < log2FC.SNCA4.h1_vs_h0) %>% pull(SNCA2_lt_SNCA4) %>% table()
    #mutate(SNCA2_minus_SNCA4 = log2FC.SNCA2.h1_vs_h0-log2FC.SNCA4.h1_vs_h0) %>% pull(SNCA2_minus_SNCA4) %>% t.test()
    mutate(SNCA2_minus_SNCA4 = log2FC.SNCA2.h1_vs_h0-log2FC.SNCA4.h1_vs_h0) 
  ggplot(df2, aes(x=SNCA2_minus_SNCA4, fill=ifelse(SNCA2_minus_SNCA4<0, paste0("SNCA2<SNCA4: ",round(100*mean(SNCA2_minus_SNCA4<0),2),"%"), paste0("SNCA2>SNCA4: ",round(100*mean(SNCA2_minus_SNCA4>0),2),"%")))) + 
    geom_histogram(colour="black", position="identity", alpha=0.3, breaks=seq(-1,1,0.05)) + xlab("log2FC.SNCA2.h1_vs_h0 - log2FC.SNCA4.h1_vs_h0") + 
    labs(fill="") + theme(legend.position="top") + 
    ggtitle(paste("Slope difference between SNCA2 vs. SNCA4 at h1 for the monotonic decreasing genes (n =", nrow(df %>% filter(monotonic_SNCA2, monotonic_SNCA4)),")")) +
    ggsave(file.path(output_dir, paste0("DEresult.decayDiff.histogram.grouped.", pff ,".monotonic.h1diff.pdf")), width = 10, height = 4)
  
  # # paired difference in slope
  # rownames_to_column(res) %>% 
  #   pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
  #   separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
  #   pivot_wider(names_from = CONDITION, values_from = log2FC) %>%
  #   filter(SNCA2<0, SNCA4<0) %>%
  #   mutate(delta = SNCA2 - SNCA4) %>%
  #   mutate(timepoint_comparison=factor(timepoint_comparison, levels = paste0("h",c(1,3,6,12),"_vs_h0"))) %>%
  #   ggplot(aes(x=delta)) + 
  #   #geom_histogram(fill="white", alpha=0.5, position="identity", bins = 100) +
  #   geom_histogram(aes(y=..density..), colour="black", position="identity", alpha=0.3, bins=200)+
  #   geom_density(alpha=.2) +
  #   facet_wrap(vars(timepoint_comparison), ncol = 1, nrow = 4, strip.position = "top") +
  #   ggtitle(paste0("DEresult.decayDiff.histogram.paired.", pff ,".all.pdf")) + 
  #   ggsave(file.path(output_dir, paste0("DEresult.decayDiff.histogram.paired.", pff ,".all.pdf")))
  
  ## count plot for individual significant genes
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
  pdf(file.path(output_dir, paste0("DEresult.interactionLRT.", pff ,".padj05.plotCounts.pdf")), width=8, height = 10, paper = 'letter')
  for(i in 1:required_n_pages){
    ggplot(dd, aes(x = Timepoint, y = count, color = CONDITION, group = CONDITION)) + 
      geom_point() + stat_summary(fun.y=mean, geom="line") +
      labs(title=paste0("Genes with significant condition-specific changes over time (LRT padj < 0.05)"), caption = paste("page",i,"of",required_n_pages),
           subtitle=file.path(output_dir, paste0("DEresult.interactionLRT.", pff ,".padj05.pdf")), y="Normalized count (log10 scalue)") +
      scale_y_log10() + facet_wrap_paginate(vars(geneSymbol), scales = "free_y", ncol = 3, nrow = 5, page =i, strip.position = "top") +
      theme(strip.background = element_blank(), strip.placement = "outside") -> p
    print(p)
  }
  dev.off()
  
  ## save the normalized count to a file
  df = dd %>% group_by(Timepoint, CONDITION, geneID, geneSymbol) %>% 
    summarise(mean_normalized_count=mean(count)) %>% ungroup() %>%
    pivot_wider(names_from=Timepoint, values_from=mean_normalized_count) %>%
    mutate(log2FC_h1_vs_h0 = log2(h1/h0), 
           log2FC_h3_vs_h1 = log2(h3/h1),
           log2FC_h6_vs_h3 = log2(h6/h3),
           log2FC_h12_vs_h6= log2(h12/h6)) %>% 
    arrange(geneID, CONDITION)
  write.table(df, 
              file=file.path(output_dir, paste0("DEresult.interactionLRT.", pff ,".padj05.consective_log2FC.xls")), 
              sep="\t", quote =F, na="", row.names = F)
  
  write.table(t(cnts),  
              file=file.path(output_dir, paste0("DEresult.interactionLRT.", pff ,".padj05.normalizedCount.xls")), 
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
  annotation_col = data.frame(rowname=names(topDE), CONDITION=sub("\\..*","",names(topDE))) %>% column_to_rownames()
  ann_colors = list(
    #geneType = c(protein_coding = "darkblue", lincRNA = "orange", Mt_rRNA='pink', antisense='yellow', pseudogene='gray', processed_transcript='lightgray',Mt_tRNA='purple',misc_RNA='red',snoRNA='lightblue', polymorphic_pseudogene="gray", processed_pseudogene="gray", TEC="lightred"),
    CONDITION = c(SNCA2 = "#F8766D", SNCA4 = "#00BFC4")
  )
  par(cex=0.5, mar=c(5, 8, 4, 1))
  pheatmap(topDE,
           fontsize = 8,
           main =paste("Heatmap of log2 fold changes for genes with adjusted p < 0.05\n",
                       file.path(output_dir, paste0("DEresult.interactionLRT.", pff ,".all.xls.gz"))),
           width=6, height=0.15*nrow(topDE),
           filename=file.path(output_dir, paste0("DEresult.interactionLRT.", pff ,".topP.heatmap.pdf")),
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
  
  # test for each time point 
  for(hx in levels(covarianceTable$Timepoint)[-1]) {
    
    comparison_name = paste0("CONDITION", levels(covarianceTable$CONDITION)[-1], ".", "Timepoint", hx)
    message(comparison_name)
    res <- results(dds_subset, name = comparison_name,  parallel=TRUE, BPPARAM=MulticoreParam(4), test = "Wald") # now the p-value is for testing each FC
    
    # add annotation
    res$symbol <- genes_annotation$geneSymbol[match(row.names(res), genes_annotation$geneID)]
    res$geneType <- genes_annotation$geneType[match(row.names(res), genes_annotation$geneID)]
    res$geneDescription <- genes_annotation$geneDescription[match(row.names(res), genes_annotation$geneID)]
    
    ## save the data
    res <- res[order(res$padj),]
    head(res); dim(res)
    
    write.table(as.data.frame(res), 
                file=gzfile(file.path(output_dir, paste0("DEresult.interactionWald.", pff, ".", comparison_name,".all.xls.gz"))), 
                sep="\t", quote =F, na="", row.names=T, col.names = NA)
    
    write.table(as.data.frame(subset(res, padj<=0.05)), 
                file=file.path(output_dir, paste0("DEresult.interactionWald.", pff, ".", comparison_name,".padj05.xls")), 
                sep="\t", quote =F, na="", row.names=T, col.names = NA)
    
    
    
    ## MAKING PLOTS
    pdf(file.path(output_dir, paste0("DEresult.interactionWald.", pff,".", comparison_name,".pdf")), paper = 'USr')
    
    DESeq2::plotMA(res, alpha = 0.05, colNonSig = "gray", main=paste0("DEresult.interactionWald.", pff,".", comparison_name))
    
    ## enhanced Vocano plot
    library(EnhancedVolcano)
    p=EnhancedVolcano(res,
                      subtitle = paste0("DEresult.interactionWald.", pff,".", comparison_name),
                      lab = res$symbol,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      drawConnectors = TRUE,
                      arrowheads =F)
    print(p)
    
    dev.off()
    
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
  res_negative_h1 = filter(res, log2FC.PFFno.h1_vs_h0<0, log2FC.PFFyes.h1_vs_h0<0)
  ttest_pvalue = rownames_to_column(res_negative_h1) %>% 
    pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
    pivot_wider(names_from = PFF.status, values_from = log2FC) %>% 
    filter(pvalue<=1) %>%
    group_by(timepoint_comparison) %>% 
    summarise(PFFno=list(PFFno), PFFyes=list(PFFyes)) %>% rowwise() %>% 
    mutate(pvalue=wilcox.test(unlist(PFFno), unlist(PFFyes), paired = T, alternative = "less")$p.value) %>% # alternative = 'less'?
    mutate(timepoint_comparison=factor(timepoint_comparison, levels = paste0("h",c(1,3,6,12),"_vs_h0"))) 
  
  nDecay = rownames_to_column(res_negative_h1) %>% 
    pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
    pivot_wider(names_from = PFF.status, values_from = log2FC) %>% 
    filter(pvalue<=1) %>%
    group_by(timepoint_comparison) %>% 
    summarise(PFFno.percent=mean(PFFno<0), PFFyes.percent=mean(PFFyes<0), PFFno.n=sum(PFFno<0), PFFyes.n=sum(PFFyes<0)) %>% 
    pivot_longer(cols = contains('PFF'), names_to = c('PFF.status',".value"), names_pattern = "(.*)\\.(.*)") %>% 
    mutate(timepoint_comparison=factor(timepoint_comparison, levels = paste0("h",c(1,3,6,12),"_vs_h0"))) 
  
  pivot_longer(as_tibble(res_negative_h1), starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
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
    ggtitle(paste0("DEresult.decayDiff.histogram.grouped.", condition ,".res_negative_h1.pdf")) + 
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
# only subset of h1 and h0
###########################################

for(pff in levels(covarianceTable$PFF.status)){
  message(paste("# Running DEseq on PFF status =", pff, "..."))
  # subsetting
  # debug: pff="PFFno"
  
  dds_subset=dds[,dds$PFF.status==pff & dds$Timepoint %in% c("h0","h1")]; 
  dim(dds_subset);
  
  # re-factorize 
  colData(dds_subset)$Timepoint = factor(colData(dds_subset)$Timepoint, levels = c("h0","h1"))
  
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
  
  # ## Wald tests for the log2 fold changes at individual time points
  # res_h1 = results(dds_subset, name="CONDITIONSNCA4.Timepointh1", test="Wald")
  # res_h1$symbol <- genes_annotation$geneSymbol[match(row.names(res_h1), genes_annotation$geneID)]
  # res_h1$geneType <- genes_annotation$geneType[match(row.names(res_h1), genes_annotation$geneID)]
  # res_h1 <- res_h1[order(res_h1$padj),]
  # head(res_h1)
  
  # extract a matrix of the log2 fold changes for all comparisons using the coef function.
  betas <- coef(dds_subset)
  colnames(betas)  # same as the resultsNames(dds_subset)
  betas = betas[rownames(res), -c(1,2)]
  ## Note that CONDITIONSNCA4.Timepointh1 means difference between SNCA4 vs SNCA2 at h1, controlling for baseline. 
  ## To get the log2 fold change of h1 vs h0 for the SNCA2, Timepoint_h1_vs_h0
  ## To get the log2 fold change of h1 vs h0 for the SNCA4, Timepoint_h1_vs_h0 +  CONDITIONSNCA4.Timepointh1
  betas = as.data.frame(betas) %>% mutate(.keep = "none",
                                          SNCA2.h1_vs_h0  = Timepoint_h1_vs_h0, 
                                          SNCA4.h1_vs_h0  = Timepoint_h1_vs_h0 + CONDITIONSNCA4.Timepointh1)
  
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
              file=gzfile(file.path(output_dir, paste0("DEresult.interactionLRT.", pff ,".h1only.all.xls.gz"))), 
              sep="\t", quote =F, na="", row.names=T, col.names = NA)
  
  write.table(subset(res, padj<=0.05), 
              file=file.path(output_dir, paste0("DEresult.interactionLRT.", pff ,".h1only.padj05.xls")), 
              sep="\t", quote =F, na="", row.names=T, col.names = NA)
  
  ## any globe difference in slope between SNCA2 vs. SNCA4
  ttest_pvalue = rownames_to_column(res) %>% 
    pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
    pivot_wider(names_from = CONDITION, values_from = log2FC) %>% 
    filter(pvalue<=1) %>%
    mutate(timepoint_comparison=factor(timepoint_comparison, levels = paste0("h",c(1,3,6,12),"_vs_h0"))) %>%
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
  
  # pff=no
  # timepoint_comparison SNCA2          SNCA4          pvalue.paired_ttest_less pvalue.paired_ttest pvalue.ttest_less pvalue.ttest pvalue.paired_wilcox_less pvalue.paired_wilcox pvalue.wilcox_less pvalue.wilcox
  # <fct>                <list>         <list>                            <dbl>               <dbl>             <dbl>        <dbl>                     <dbl>                <dbl>              <dbl>         <dbl>
  #   1 h1_vs_h0             <dbl [26,207]> <dbl [26,207]>                    0.942               0.116             0.914        0.172                     0.994               0.0118              0.972        0.0559
  
  # pff=yes
  # timepoint_comparison SNCA2          SNCA4          pvalue.paired_ttest_less pvalue.paired_ttest pvalue.ttest_less pvalue.ttest pvalue.paired_wilcox_less pvalue.paired_wilcox pvalue.wilcox_less pvalue.wilcox
  # <fct>                <list>         <list>                            <dbl>               <dbl>             <dbl>        <dbl>                     <dbl>                <dbl>              <dbl>         <dbl>
  #   1 h1_vs_h0             <dbl [26,210]> <dbl [26,210]>                   0.0926               0.185             0.130        0.259                     0.964               0.0716              0.544         0.911
  
  nDecay = rownames_to_column(res) %>% 
    pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
    pivot_wider(names_from = CONDITION, values_from = log2FC) %>% 
    filter(pvalue<=1) %>%
    group_by(timepoint_comparison) %>% 
    summarise(SNCA2.percent=mean(SNCA2<0), SNCA4.percent=mean(SNCA4<0), SNCA2.n=sum(SNCA2<0), SNCA4.n=sum(SNCA4<0)) %>% 
    pivot_longer(cols = contains('SNCA'), names_to = c('CONDITION',".value"), names_pattern = "(.*)\\.(.*)") 
  
  pivot_longer(as_tibble(res), starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
    filter(pvalue<=1) %>%
    ggplot(aes(x=log2FC, fill=CONDITION, color=CONDITION)) + 
    #geom_histogram(fill="white", alpha=0.5, position="identity", bins = 100) +
    geom_histogram(aes(y=..density..), colour="black", position="identity", alpha=0.3, bins=100)+
    geom_density(alpha=.2) +
    geom_vline(xintercept =0, color='blue', linetype=2) +
    #facet_wrap(vars(timepoint_comparison), ncol = 1, nrow = 4, strip.position = "top") + 
    geom_text(data = ttest_pvalue, mapping = aes(x = -Inf, y = Inf, label = paste("Paired one-sided Wilcox test p = ",formatC(pvalue, format = "g", digits = 2)), color = NULL,fill= NULL), hjust   = -.05, vjust   = 4.5) +
    geom_text(data = nDecay, mapping = aes(x = -Inf, y = Inf, label = paste0("Decreasing genes in ",CONDITION, ": n=",n, " (",round(100*percent, 2), "%)"), color = CONDITION), hjust   = -0.05, vjust   = rep(c(1.5, 3), length(unique(nDecay$timepoint_comparison)))) +
    ggtitle(paste0("DEresult.decayDiff.histogram.grouped.", pff ,".h1only.all.pdf")) + 
    ggsave(file.path(output_dir, paste0("DEresult.decayDiff.histogram.grouped.", pff ,".h1only.all.pdf")), width = 11, height = 3)
  
  # histogram for negative slope genes @ h1
  res_negative_h1 = filter(res, log2FC.SNCA2.h1_vs_h0<0, log2FC.SNCA4.h1_vs_h0<0)
  ttest_pvalue = rownames_to_column(res_negative_h1) %>% 
    pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
    pivot_wider(names_from = CONDITION, values_from = log2FC) %>% 
    filter(pvalue<=1) %>%
    mutate(timepoint_comparison=factor(timepoint_comparison, levels = paste0("h",c(1,3,6,12),"_vs_h0"))) %>%
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
  
  # pff=no
  # timepoint_comparison SNCA2         SNCA4         pvalue.paired_ttest_less pvalue.paired_ttest pvalue.ttest_less pvalue.ttest pvalue.paired_wilcox_less pvalue.paired_wilcox pvalue.wilcox_less pvalue.wilcox
  # <fct>                <list>        <list>                           <dbl>               <dbl>             <dbl>        <dbl>                     <dbl>                <dbl>              <dbl>         <dbl>
  #   1 h1_vs_h0             <dbl [9,561]> <dbl [9,561]>            0.00000000576        0.0000000115         0.0000502     0.000100                  1.58e-50             3.16e-50   0.00000000000123      2.45e-12
  
  # pff=yes
  # timepoint_comparison SNCA2         SNCA4         pvalue.paired_ttest_less pvalue.paired_ttest pvalue.ttest_less pvalue.ttest pvalue.paired_wilcox_less pvalue.paired_wilcox pvalue.wilcox_less pvalue.wilcox
  # <fct>                <list>        <list>                           <dbl>               <dbl>             <dbl>        <dbl>                     <dbl>                <dbl>              <dbl>         <dbl>
  #   1 h1_vs_h0             <dbl [9,878]> <dbl [9,878]>                    0.279               0.559             0.351        0.701                    0.0253               0.0506              0.372         0.744
  
  nDecay = rownames_to_column(res_negative_h1) %>% 
    pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
    pivot_wider(names_from = CONDITION, values_from = log2FC) %>% 
    filter(pvalue<=1) %>%
    group_by(timepoint_comparison) %>% 
    summarise(SNCA2.percent=mean(SNCA2<0), SNCA4.percent=mean(SNCA4<0), SNCA2.n=sum(SNCA2<0), SNCA4.n=sum(SNCA4<0)) %>% 
    pivot_longer(cols = contains('SNCA'), names_to = c('CONDITION',".value"), names_pattern = "(.*)\\.(.*)") 
  
  pivot_longer(as_tibble(res_negative_h1), starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
    filter(pvalue<=1) %>%
    ggplot(aes(x=log2FC, fill=CONDITION, color=CONDITION)) + 
    #geom_histogram(fill="white", alpha=0.5, position="identity", bins = 100) +
    geom_histogram(aes(y=..density..), colour="black", position="identity", alpha=0.3, bins=100)+
    geom_density(alpha=.2) +
    geom_vline(xintercept =0, color='blue', linetype=2) +
    facet_wrap(vars(timepoint_comparison), ncol = 1, nrow = 4, strip.position = "top") + 
    geom_text(data = ttest_pvalue, mapping = aes(x = -Inf, y = Inf, label = paste("Paired one-sided Wilcox test p = ",formatC(pvalue, format = "g", digits = 2)), color = NULL,fill= NULL), hjust   = -0.05, vjust   = 4.5) +
    geom_text(data = nDecay, mapping = aes(x = -Inf, y = Inf, label = paste0("Decreasing genes in ",CONDITION, ": n=",n, " (",round(100*percent, 2), "%)"), color = CONDITION), hjust   = -0.05, vjust   = rep(c(1.5, 3), length(unique(nDecay$timepoint_comparison)))) +
    ggtitle(paste0("DEresult.decayDiff.histogram.grouped.", pff ,".h1only.res_negative_h1.pdf")) + 
    ggsave(file.path(output_dir, paste0("DEresult.decayDiff.histogram.grouped.", pff ,".h1only.subset_negative_h1.pdf")), width = 11, height = 3)
  
  # if the decreasing ones have bigger decay rate in SNCA2 vs. SNCA4
  rownames_to_column(res_negative_h1) %>% 
    mutate(SNCA2_minus_SNCA4 = log2FC.SNCA2.h1_vs_h0-log2FC.SNCA4.h1_vs_h0) %>% 
    ggplot(aes(x=SNCA2_minus_SNCA4, fill=SNCA2_minus_SNCA4<0)) + 
    geom_histogram(colour="black", position="identity", alpha=0.3, breaks=seq(-3,3,0.1)) + 
    xlab("log2FC.SNCA2.h1_vs_h0 - log2FC.SNCA4.h1_vs_h0") + 
    ggtitle(paste("Slope difference between SNCA2 vs. SNCA4 at h1 for the decreasing genes (n =", nrow(res_negative_h1),")")) +
    ggsave(file.path(output_dir, paste0("DEresult.decayDiff.histogram.grouped.", pff ,".h1only.res_negative_h1.h1diff.pdf")), width = 10, height = 4)
  
  # # paired difference in slope
  # rownames_to_column(res) %>% 
  #   pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
  #   separate(name, c("CONDITION","timepoint_comparison"), sep ="\\.") %>% 
  #   pivot_wider(names_from = CONDITION, values_from = log2FC) %>%
  #   filter(SNCA2<0, SNCA4<0) %>%
  #   mutate(delta = SNCA2 - SNCA4) %>%
  #   mutate(timepoint_comparison=factor(timepoint_comparison, levels = paste0("h",c(1,3,6,12),"_vs_h0"))) %>%
  #   ggplot(aes(x=delta)) + 
  #   #geom_histogram(fill="white", alpha=0.5, position="identity", bins = 100) +
  #   geom_histogram(aes(y=..density..), colour="black", position="identity", alpha=0.3, bins=200)+
  #   geom_density(alpha=.2) +
  #   facet_wrap(vars(timepoint_comparison), ncol = 1, nrow = 4, strip.position = "top") +
  #   ggtitle(paste0("DEresult.decayDiff.histogram.paired.", pff ,".all.pdf")) + 
  #   ggsave(file.path(output_dir, paste0("DEresult.decayDiff.histogram.paired.", pff ,".all.pdf")))
  
  ## now we can cluster the top significant genes (in h1 timepoint) by their Log2FC profiles.
  topN=1000
  topDE = filter(res, log2FC.SNCA2.h1_vs_h0<0, padj<=0.05) %>% head(topN) %>% as_tibble() %>% column_to_rownames('symbol') %>%  select(contains("log2FC")) %>% arrange(log2FC.SNCA2.h1_vs_h0)
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
  par(cex=0.5, mar=c(5, 8, 4, 1))
  pheatmap(topDE,
           fontsize = 8,
           main =paste("Heatmap of log2 fold changes for genes with adjusted p < 0.05\n",
                       file.path(output_dir, paste0("DEresult.interactionLRT.", pff ,".h1only.all.xls.gz"))),
           width=6, height=0.15*nrow(topDE),
           filename=file.path(output_dir, paste0("DEresult.interactionLRT.", pff ,".h1only.topP.heatmap.pdf")),
           #border_color = NA,
           color=myColor, breaks=myBreaks,
           annotation_row = annotation_row,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           drop_levels = TRUE,
           scale = "none", border_color = F,
           #clustering_method = 'ward.D', 
           cluster_rows = F,
           #clustering_distance_rows = "correlation",
           cluster_cols = F)
  
}

###########################################
# run DE among PFF, Timepoint jointly
# to test any genes that react in a condition-specific manner over time
# see http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments
# for h1 and h0 only
###########################################

for(condition in levels(covarianceTable$CONDITION)){
  message(paste("# Running DEseq on CONDITION status =", condition, "..."))
  # subsetting
  # debug: condition="SNCA2"
  dds_subset=dds[,dds$CONDITION==condition & dds$Timepoint %in% c("h0","h1")]; 
  dim(dds_subset);
  
  # re-factorize 
  colData(dds_subset)$Timepoint = factor(colData(dds_subset)$Timepoint, levels = c("h0","h1"))
  
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
                                          PFFyes.h1_vs_h0  = Timepoint_h1_vs_h0 + PFF.statusPFFyes.Timepointh1)
  
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
              file=gzfile(file.path(output_dir, paste0("DEresult.interactionLRT.", condition ,".h1only.all.xls.gz"))), 
              sep="\t", quote =F, na="", row.names=T, col.names = NA)
  
  write.table(subset(res, padj<=0.05), 
              file=file.path(output_dir, paste0("DEresult.interactionLRT.", condition ,".h1only.padj05.xls")), 
              sep="\t", quote =F, na="", row.names=T, col.names = NA)
  
  ## any globe difference in slope between PFFno vs. PFFyes
  ttest_pvalue = rownames_to_column(res) %>% 
    pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
    pivot_wider(names_from = PFF.status, values_from = log2FC) %>% 
    filter(pvalue<=1) %>%
    group_by(timepoint_comparison) %>% 
    summarise(PFFno=list(PFFno), PFFyes=list(PFFyes)) %>% rowwise() %>% 
    mutate(pvalue=wilcox.test(unlist(PFFno), unlist(PFFyes), paired = T, alternative = "less")$p.value) 
  
  nDecay = rownames_to_column(res) %>% 
    pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
    pivot_wider(names_from = PFF.status, values_from = log2FC) %>% 
    filter(pvalue<=1) %>%
    group_by(timepoint_comparison) %>% 
    summarise(PFFno.percent=mean(PFFno<0), PFFyes.percent=mean(PFFyes<0), PFFno.n=sum(PFFno<0), PFFyes.n=sum(PFFyes<0)) %>% 
    pivot_longer(cols = contains('PFF'), names_to = c('PFF.status',".value"), names_pattern = "(.*)\\.(.*)") 
  
  pivot_longer(as_tibble(res), starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>%
    filter(pvalue<=1) %>%
    ggplot(aes(x=log2FC, fill=PFF.status, color=PFF.status)) + 
    #geom_histogram(fill="white", alpha=0.5, position="identity", bins = 100) +
    geom_histogram(aes(y=..density..), colour="black", position="identity", alpha=0.3, bins=100)+
    geom_density(alpha=.2) +
    geom_vline(xintercept =0, color='blue', linetype=2) +
    facet_wrap(vars(timepoint_comparison), ncol = 1, nrow = 4, strip.position = "top") + 
    geom_text(data = ttest_pvalue, mapping = aes(x = Inf, y = Inf, label = paste("Paired one-sided Wilcox test p = ",formatC(pvalue, format = "g", digits = 2)), color = NULL,fill= NULL), hjust   = 1.05, vjust   = 1.5) +
    geom_text(data = nDecay, mapping = aes(x = -Inf, y = Inf, label = paste0("Decreasing genes in ",PFF.status, ": n=",n, " (",round(100*percent, 2), "%)"), color = PFF.status), hjust   = -0.05, vjust   = c(1.5, 3)) +
    ggtitle(paste0("DEresult.decayDiff.histogram.grouped.", condition ,".h1only.all.pdf")) + 
    ggsave(file.path(output_dir, paste0("DEresult.decayDiff.histogram.grouped.", condition ,".h1only.all.pdf")), width = 11, height = 3)
  
  # histogram for negative slope genes @ h1
  res_negative_h1 = filter(res, log2FC.PFFno.h1_vs_h0<0, log2FC.PFFyes.h1_vs_h0<0)
  ttest_pvalue = rownames_to_column(res_negative_h1) %>% 
    pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
    pivot_wider(names_from = PFF.status, values_from = log2FC) %>% 
    filter(pvalue<=1) %>%
    group_by(timepoint_comparison) %>% 
    summarise(PFFno=list(PFFno), PFFyes=list(PFFyes)) %>% rowwise() %>% 
    mutate(pvalue=wilcox.test(unlist(PFFno), unlist(PFFyes), paired = T, alternative = "less")$p.value) 
  
  nDecay = rownames_to_column(res_negative_h1) %>% 
    pivot_longer(starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
    pivot_wider(names_from = PFF.status, values_from = log2FC) %>% 
    filter(pvalue<=1) %>%
    group_by(timepoint_comparison) %>% 
    summarise(PFFno.percent=mean(PFFno<0), PFFyes.percent=mean(PFFyes<0), PFFno.n=sum(PFFno<0), PFFyes.n=sum(PFFyes<0)) %>% 
    pivot_longer(cols = contains('PFF'), names_to = c('PFF.status',".value"), names_pattern = "(.*)\\.(.*)") 
  
  pivot_longer(as_tibble(res_negative_h1), starts_with("log2FC"), names_prefix = "log2FC.", values_to='log2FC') %>% 
    separate(name, c("PFF.status","timepoint_comparison"), sep ="\\.") %>% 
    filter(pvalue<=1) %>%
    ggplot(aes(x=log2FC, fill=PFF.status, color=PFF.status)) + 
    #geom_histogram(fill="white", alpha=0.5, position="identity", bins = 100) +
    geom_histogram(aes(y=..density..), colour="black", position="identity", alpha=0.3, bins=100)+
    geom_density(alpha=.2) +
    geom_vline(xintercept =0, color='blue', linetype=2) +
    facet_wrap(vars(timepoint_comparison), ncol = 1, nrow = 4, strip.position = "top") + 
    geom_text(data = ttest_pvalue, mapping = aes(x = -Inf, y = Inf, label = paste("Paired one-sided Wilcox test p = ",formatC(pvalue, format = "g", digits = 2)), color = NULL,fill= NULL), hjust   = -0.05, vjust   = 4.5) +
    geom_text(data = nDecay, mapping = aes(x = -Inf, y = Inf, label = paste0("Decreasing genes in ",PFF.status, ": n=",n, " (",round(100*percent, 2), "%)"), color = PFF.status), hjust   = -0.05, vjust   = c(1.5, 3)) +
    ggtitle(paste0("DEresult.decayDiff.histogram.grouped.", condition ,".h1only.res_negative_h1.pdf")) + 
    ggsave(file.path(output_dir, paste0("DEresult.decayDiff.histogram.grouped.", condition ,".h1only.subset_negative_h1.pdf")), width = 11, height = 3)
  
  # if the decreasing ones have bigger decay rate in PFFno vs. PFFyes
  rownames_to_column(res_negative_h1) %>%
    mutate(PFFno_minus_PFFyes = log2FC.PFFno.h1_vs_h0-log2FC.PFFyes.h1_vs_h0) %>% 
    ggplot(aes(x=PFFno_minus_PFFyes, fill=PFFno_minus_PFFyes<0)) + 
    geom_histogram(colour="black", position="identity", alpha=0.3, breaks=seq(-3,3,0.1)) + 
    xlab("log2FC.PFFno.h1_vs_h0 - log2FC.PFFyes.h1_vs_h0") + 
    ggtitle(paste("Slope difference between PFFno vs. PFFyes at h1 for the decreasing genes (n =", nrow(res_negative_h1),")")) +
    ggsave(file.path(output_dir, paste0("DEresult.decayDiff.histogram.grouped.", condition ,".h1only.res_negative_h1.h1diff.pdf")), width = 10, height = 4)
  
  
  ## now we can cluster the top significant genes (in any timepoint) by their Log2FC profiles.
  topN=1000
  topDE = filter(res, log2FC.PFFno.h1_vs_h0<0, padj<=0.05) %>% head(topN) %>% as_tibble() %>% column_to_rownames('symbol') %>%  select(contains("log2FC"))
  colnames(topDE) = gsub("log2FC.","",colnames(topDE))
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
  par(cex=0.5, mar=c(5, 8, 4, 1))
  if(nrow(topDE)>0) pheatmap(topDE,
                             fontsize = 8,
                             main =paste("Heatmap of log2FC for genes with adj.p < 0.05 and log2FC.PFFno.h1_vs_h0<0\n",
                                         file.path(output_dir, paste0("DEresult.interactionLRT.", condition ,".h1only.all.xls.gz"))),
                             width=6, height=0.3*nrow(topDE)+2,
                             filename=file.path(output_dir, paste0("DEresult.interactionLRT.", condition ,".h1only.topP.heatmap.pdf")),
                             #border_color = NA,
                             color=myColor, breaks=myBreaks,
                             annotation_row = annotation_row,
                             annotation_col = annotation_col,
                             annotation_colors = ann_colors,
                             drop_levels = TRUE,
                             scale = "none", 
                             #clustering_method = 'ward.D', 
                             cluster_rows = F,
                             #clustering_distance_rows = "correlation",
                             cluster_cols = F)
  
}
