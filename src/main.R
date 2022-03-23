library(tidyverse)

################################################
## read the data
################################################
output_dir="~/projects/vik2021/results/DE"
Timepoints_sorted_levels = c(0,1,3,6,12)
for(compare in c('PFFno', 'PFFyes', 'SNCA2', 'SNCA4')){
  # compare="PFFno";
  consecutiveDEgenes=data.frame()
  for(i in 2:length(Timepoints_sorted_levels)){
    time_i = Timepoints_sorted_levels[i];
    time_i_m_1 = Timepoints_sorted_levels[i-1]
    message(paste0("# Running pair-wise comparison between consecutive h",time_i," vs. h",time_i_m_1, " in ", compare, "...")) 
    res= read.table(file=file.path(output_dir, paste0("DEresult.interactionLRT.", compare, ".consecutive_h",time_i,"_vs_h",time_i_m_1,".all.xls.gz")), 
    #res= read.table(file=file.path(output_dir, paste0("DEresult.interactionLRT.", compare, ".consecutive_h",time_i,"_vs_h",time_i_m_1,".padj05.xls")), 
                    sep="\t", header = T)
    res = mutate(res, subset=compare, timeinterval = paste0("h",time_i,"_vs_h",time_i_m_1))
    colnames(res) = gsub(paste0(".h",time_i,"_vs_h",time_i_m_1),"", colnames(res))
    consecutiveDEgenes = rbind(consecutiveDEgenes,res)
    
    # just for each case individually
    x_name=colnames(res)[grep("log2FC",colnames(res))][1];
    y_name=colnames(res)[grep("log2FC",colnames(res))][2]
    filter(res, padj<=.05) %>% 
      mutate(geneType=ifelse(geneType %in% c("lncRNA","protein_coding"), geneType, "others")) %>% 
      mutate(geneType=factor(geneType, levels = c("protein_coding","lncRNA","others"))) %>% 
      ggplot(aes_string(x=x_name, y=y_name)) +
      geom_point(aes(size=-log10(pvalue), alpha=-log10(pvalue))) +
      geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + geom_abline(intercept = 0, slope = 1, linetype=2) + 
      #coord_fixed(ratio = 1) + theme(legend.position="top") +
      coord_fixed(ratio = 1) + theme(legend.position="right", aspect.ratio=1) +
      #facet_wrap(vars(timeinterval), ncol = 4, nrow = 1, strip.position = "top") + 
      #geom_text(aes(x = -Inf, y = Inf, label = timeinterval), hjust   = -0.05, vjust = 1.5, size = 10, color = "#999999") + 
      ggtitle(paste0("LRT significant genes in ", compare)) 
  }
  consecutiveDEgenes$timeinterval=factor(consecutiveDEgenes$timeinterval, levels = unique(consecutiveDEgenes$timeinterval))
  
  
  message("mplot...")
  x_name=colnames(consecutiveDEgenes)[grep("log2FC",colnames(consecutiveDEgenes))][1];
  y_name=colnames(consecutiveDEgenes)[grep("log2FC",colnames(consecutiveDEgenes))][2]
  # scatter plot for significant LRT genes
  df = filter(consecutiveDEgenes, pvalue<=.05) %>% 
    mutate(geneType=ifelse(geneType %in% c("lncRNA","protein_coding"), geneType, "others")) %>% 
    mutate(geneType=factor(geneType, levels = c("protein_coding","lncRNA","others")))
  ggplot(df, aes_string(x=x_name, y=y_name)) +
    geom_point(aes(size=-log10(pvalue), alpha=-log10(pvalue))) +
    geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + geom_abline(intercept = 0, slope = 1, linetype=2) + 
    #coord_fixed(ratio = 1) + theme(legend.position="top") +
    coord_fixed(ratio = 1, xlim=range(select(df, contains("log2FC"))), ylim =range(select(df, contains("log2FC")))) + theme(legend.position="right", aspect.ratio=1) +
    facet_wrap(vars(timeinterval), ncol = 4, nrow = 1, strip.position = "top") + 
    #geom_text(aes(x = -Inf, y = Inf, label = timeinterval), hjust   = -0.05, vjust = 1.5, size = 10, color = "#999999") + 
    ggtitle(paste0("LRT significant genes in ", compare)) 
  ggsave(file.path(output_dir, paste0("DEresult.interactionLRT.", compare, ".consecutive.p05.scatterplot.pdf")), width = 12, height = 4)
}


x_name=colnames(res)[grep("log2FC",colnames(res))][1];
y_name=colnames(res)[grep("log2FC",colnames(res))][2]
# scatter plot for significant LRT genes
filter(res, pvalue<=.05) %>% 
  mutate(geneType=ifelse(geneType %in% c("lncRNA","protein_coding"), geneType, "others")) %>% 
  mutate(geneType=factor(geneType, levels = c("protein_coding","lncRNA","others"))) %>% 
  ggplot(aes_string(x=x_name, y=y_name)) +
  geom_point(aes(size=-log10(pvalue), alpha=-log10(pvalue), col=geneType)) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + geom_abline(intercept = 0, slope = 1, linetype=2) + 
  coord_fixed(ratio = 1) +
  geom_text(x = -Inf, y = Inf, hjust   = -0.05, vjust = 1.5, size = 10, color = "#999999",label = paste0("h",time_i," vs h",time_i_m_1), alpha = .5) + 
  ggtitle(paste0("LRT significant genes in ", X, " between consecutive h",time_i," vs h",time_i_m_1)) 
ggsave(file.path(output_dir, paste0("DEresult.interactionLRT.", X, ".consecutive_h",time_i,"_vs_h",time_i_m_1,".padj05.scatterplot.pdf")), width = 8, height = 7)

consecutiveDEgenes = rbind(mutate(res, subset=X, timeinterval = paste0("h",time_i,"_vs_h",time_i_m_1)), consecutiveDEgenes)


library(plotly)

p<- fert%>%
  ggplot(aes(x=fert, y=life, size = pop, color = continent,frame = year)) +
  labs(x="Fertility Rate", y = "Life expectancy at birth (years)", 
       caption = "(Based on data from Hans Rosling - gapminder.com)", 
       color = 'Continent',size = "Population (millions)") + 
  ylim(30,100) +
  geom_point(aes(text=Country))

ggplotly(p)

