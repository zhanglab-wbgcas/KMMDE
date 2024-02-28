library(reticulate)
library(dplyr)
library(pROC)
require(plotROC)
library(ggplot2)
#source('/home/liukc/project/sdc/02func_algori/revise_dre/timedata_pro/r_compare_proj/FunPat/R/SEL.TS.AREA.R') 
#use_python("/bin/python") #
rm(list=ls())
setwd('E:\\差异贡献\\sdc函数\\修改\\模拟数据\\时序数据\\rna-seq_tcComp\\teststimdata\\new_timedata\\roc_plots')
replicate.range <-c(2)
deg_rang <- c(240,480,960,2400,3600)
for(rep_num in replicate.range){
  result <- read.csv(file=paste("rocplot_tibble_p",rep_num,"res.csv",sep="_"))  # 设置 row.names = FALSE 以避免保存行名
  ####plot####
  
  colnames(result)[2] <- 'Methods'
  head(result);tail(result)
  result <- result %>% filter(!is.na(value))
  
  for(deg_num in deg_rang){
    result$DEG[result$deg_numbers == deg_num] <- 1
    result$DEG[result$deg_numbers == deg_num & result$gene <= deg_num] <- 0
    result$DEG[result$Methods == "KMMDDE" & result$deg_numbers == deg_num] <- 0
    result$DEG[result$Methods == "KMMDDE" & result$deg_numbers == deg_num& result$gene <= deg_num] <- 1
    
  }
  result$Methods <- factor(result$Methods,
                           levels = c("KMMDE", "edgeR","DESeq2","maSigPro","splineTC", "Funpat","ImpulseDE2"))  #"maSigPro","splineTC",
  
  g <- ggplot(result, aes(d = DEG, m = ifelse(Methods == "KMMDE", 1-value, value), colour = Methods)) + scale_color_manual(values =c("red","blue","orange", "darkgreen","grey","cyan","purple")) + theme(legend.text = element_text(size = 14))
  
  g <- g + geom_roc(n.cuts = FALSE, linealpha = 0.9,size = 0.6)
  g <- g + xlab("False positive rate") + ylab("True positive rate") #theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
  # g <- g + theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_grid(deg_numbers ~ timepoints) +theme(axis.text.y = element_text(size = 10))+ theme(axis.text.x = element_text(size = 10,angle = 90, hjust = 1))+theme(strip.text.x = element_text(size = 12))+theme(strip.text.y = element_text(size = 12))+ theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
  g <- g + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_grid(deg_numbers ~ timepoints) +theme(axis.text.y = element_text(size = 12, color = "black"))+ theme(axis.text.x = element_text(size = 12,angle = 90, hjust = 1, color = "black"))+theme(strip.text.x = element_text(size = 13))+theme(strip.text.y = element_text(size = 13))+
    theme(axis.title.x = element_text(size = 17,margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 17, margin = margin(t = 0, r = 10, b = 0, l = 0))) 
  
  g  #
  # ggsave(g, file = "new2rep_rocplot.png" ,width = 8, height =7, dpi = 400)
  ggsave(g, file = "02new_2rep_rocplot.png" ,width =9, height = 8, dpi = 300)
  
  # dev.off()
}
g  #
