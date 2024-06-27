# include generate stimddata and norm data
#2019Comparative analysis of differential gene expression tools for RNA sequencing time course data
#https://github.com/daniel-spies/rna-seq_tcComp
############################################# library loading ###########################################
library(edgeR)
library(dplyr)
############################################# data loading ###########################################

Sys.getlocale()
Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")

stimdataFunction <- function(nRep,nTP,deg){
  cheung <- read.table("/home/liukc/project/sdc/02func_algori/revise_dre/timedata_pro/r_compare_proj/cheung_count_table.txt",header=TRUE,stringsAsFactors=FALSE,row.names=1)
  
  cheung <- cheung[rowSums(cheung)>0,]
  ############################################# estimating dispersions ###########################################
  f <- DGEList(cheung)
  f <- estimateGLMCommonDisp(f) 
  f <- estimateGLMTrendedDisp(f)
  f <- estimateGLMTagwiseDisp(f)
  ############################################# create mu / disp dictionary ###########################################
  pair_f <- data.frame(mu=round(rowMeans(f$counts)),d=f$tagwise.dispersion)
  pair_f <- pair_f[order(pair_f[,1],decreasing=TRUE),]
  pair_f <- aggregate(.~mu,data=pair_f,FUN=mean)
  ############################################# simulation parameters ###########################################
  nSim <- 30000
  # multTP <- c(2,3)  # multiplicator for time points, higher sampling rate, in theory just more replicates but depicted as time points
  multTP <- 1
  ################################### simulating genes by drawing from a neg. binomial distribution ####################################
  
  ################################### simulating patterns ####################################
  
  # multipliers for patterns
  low_up <- c(2,1.2)
  low_down <- c(1/2,1/1.2)
  high_up <- c(6,3.5)
  high_down <- c(1/6,1/3.5)
  early_slow_up <- c(2,6,6,6)
  early_slow_down <- c(1/2,1/6,1/6,1/6)
  mid_slow_up <- c(2,6,6)
  mid_slow_down <- c(1/2,1/6,1/6)
  mid_fast_up <- c(6,6,6)
  mid_fast_down <- c(1/6,1/6,1/6)
  grad_lin_up <- c(2,4,6)
  grad_lin_down <- c(1/2,1/4,1/6)
  grad_cur_up <- c(2,3.5,4.5)
  grad_cur_down <- c(1/2,1/3.5,1/4.5)
  up_down_slow <- c(6,3,1.2)
  up_down_fast <- c(6,1.2)
  
  # column index of DEGs
  col <- list(
    up_early_low = c(1,1),   # selecting transcripts (rows) for DE over 2 timepoints
    up_early_high = c(1,1),
    up_mid_low = c(2,2),
    up_mid_high = c(2,2),
    up_late_low = c(3,3),
    up_late_high = c(3,3),
    down_early_low = c(1,1),
    down_early_high = c(1,1),
    down_mid_low = c(2,2),
    down_mid_high = c(2,2),
    down_late_low = c(3,3),
    down_late_high = c(3,3),
    up_early_slow = c(1,1,1,1),   # selecting transcripts (rows) for DE over 3 timepoints
    up_mid_slow = c(2,2,2),
    down_early_slow = c(1,1,1,1),
    down_mid_slow= c(2,2,2),
    up_mid_fast = c(2,2,2),
    down_mid_fast= c(2,2,2),
    up_grad_lin = c(2,2,2),
    up_grad_cur = c(2,2,2),
    down_grad_lin = c(2,2,2),
    down_grad_cur = c(2,2,2),
    up_down_slow = c(2,2,2),
    up_down_fast = c(2,2))
  
  # pattern modification of each category
  mod <- list(
    up_early_low = low_up,   # selecting transcripts (rows) for DE over 2 timepoints
    up_early_high = high_up,
    up_mid_low = low_up,
    up_mid_high = high_up,
    up_late_low = low_up,
    up_late_high = high_up,
    down_early_low = low_down,
    down_early_high = high_down,
    down_mid_low = low_down,
    down_mid_high = high_down,
    down_late_low = low_down,
    down_late_high = high_down,
    up_early_slow = early_slow_up,   # selecting transcripts (rows) for DE over 3 timepoints
    up_mid_slow = mid_slow_up,
    down_early_slow = early_slow_down,
    down_mid_slow= mid_slow_down,
    up_mid_fast = mid_fast_up,
    down_mid_fast= mid_fast_down,
    up_grad_lin = grad_lin_up,
    up_grad_cur = grad_cur_up,
    down_grad_lin = grad_lin_down,
    down_grad_cur = grad_cur_down,
    up_down_slow = up_down_slow,
    up_down_fast = up_down_fast)
  
  print(nRep);print(nTP);print(deg)
  mu <- rnbinom(n=nSim,mu=mean(rowMeans(f$counts)),size=f$common.dispersion)
  idx <- sapply(mu,function(x) which.min(abs(pair_f[,1]-min(x,max(pair_f[,1])))))
  counts <- data.frame(t(sapply(idx,function(x) replicate(nTP,round(mean(rnbinom(n=min(round(max(pair_f[x,"mu"]/4,1)),500),mu=pair_f[x,"mu"],size=pair_f[x,"d"])))))))
  counts <- counts * 4 # increase library depth
  data <- counts[sample(1:nrow(counts),20000),] # sampling 20k genes
  data <- data[rowSums(data)>0,]
  head(data)
  if(deg==240){
    # row index of DEGs 240
    row <- list(
      up_early_low = seq(1,10),   # selecting transcripts (rows) for DE over 2 timepoints
      up_early_high = seq(11,20),
      up_mid_low = seq(21,30),
      up_mid_high = seq(31,40),
      up_late_low = seq(41,50),
      up_late_high = seq(51,60),
      down_early_low = seq(61,70),
      down_early_high = seq(71,80),
      down_mid_low = seq(81,90),
      down_mid_high = seq(91,100),
      down_late_low = seq(101,110),
      down_late_high = seq(111,120),
      up_early_slow = seq(121,130),   # selecting transcripts (rows) for DE over 3 timepoints
      up_mid_slow = seq(131,140),
      down_early_slow = seq(141,150),
      down_mid_slow= seq(151,160),
      up_mid_fast = seq(161,170),
      down_mid_fast= seq(171,180),
      up_grad_lin = seq(181,190),
      up_grad_cur = seq(191,200),
      down_grad_lin = seq(201,210),
      down_grad_cur = seq(211,220),
      up_down_slow = seq(221,230),
      up_down_fast = seq(231,240))
    
  }
  else if(deg==480){
    # row index of DEGs 480
    row <- list(
      up_early_low = seq(1,20),   # selecting transcripts (rows) for DE over 2 timepoints
      up_early_high = seq(21,40),
      up_mid_low = seq(41,60),
      up_mid_high = seq(61,80),
      up_late_low = seq(81,100),
      up_late_high = seq(101,120),
      down_early_low = seq(121,140),
      down_early_high = seq(141,160),
      down_mid_low = seq(161,180),
      down_mid_high = seq(181,200),
      down_late_low = seq(201,220),
      down_late_high = seq(221,240),
      up_early_slow = seq(241,260),   # selecting transcripts (rows) for DE over 3 timepoints
      up_mid_slow = seq(261,280),
      down_early_slow = seq(281,300),
      down_mid_slow= seq(301,320),
      up_mid_fast = seq(321,340),
      down_mid_fast= seq(341,360),
      up_grad_lin = seq(361,380),
      up_grad_cur = seq(381,400),
      down_grad_lin = seq(401,420),
      down_grad_cur = seq(421,440),
      up_down_slow = seq(441,460),
      up_down_fast = seq(461,480))
  }
  else if(deg==960){
    # row index of DEGs 960
    row <- list(
      up_early_low = seq(1,40),   # selecting transcripts (rows) for DE over 2 timepoints
      up_early_high = seq(41,80),
      up_mid_low = seq(81,120),
      up_mid_high = seq(121,160),
      up_late_low = seq(161,200),
      up_late_high = seq(201,240),
      down_early_low = seq(241,280),
      down_early_high = seq(281,320),
      down_mid_low = seq(321,360),
      down_mid_high = seq(361,400),
      down_late_low = seq(401,440),
      down_late_high = seq(441,480),
      up_early_slow = seq(481,520),   # selecting transcripts (rows) for DE over 3 timepoints
      up_mid_slow = seq(521,560),
      down_early_slow = seq(561,600),
      down_mid_slow= seq(601,640),
      up_mid_fast = seq(641,680),
      down_mid_fast= seq(681,720),
      up_grad_lin = seq(721,760),
      up_grad_cur = seq(761,800),
      down_grad_lin = seq(801,840),
      down_grad_cur = seq(841,880),
      up_down_slow = seq(881,920),
      up_down_fast = seq(921,960))
  }
  else if(deg==1200){
    # row index of DEGs 1200
    row <- list(
      up_early_low = seq(1,50),   # selecting transcripts (rows) for DE over 2 timepoints
      up_early_high = seq(51,100),
      up_mid_low = seq(101,150),
      up_mid_high = seq(151,200),
      up_late_low = seq(201,250),
      up_late_high = seq(251,300),
      down_early_low = seq(301,350),
      down_early_high = seq(351,400),
      down_mid_low = seq(401,450),
      down_mid_high = seq(451,500),
      down_late_low = seq(501,550),
      down_late_high = seq(551,600),
      up_early_slow = seq(601,650),   # selecting transcripts (rows) for DE over 3 timepoints
      up_mid_slow = seq(651,700),
      down_early_slow = seq(701,750),
      down_mid_slow= seq(751,800),
      up_mid_fast = seq(801,850),
      down_mid_fast= seq(851,900),
      up_grad_lin = seq(901,950),
      up_grad_cur = seq(951,1000),
      down_grad_lin = seq(1001,1050),
      down_grad_cur = seq(1051,1100),
      up_down_slow = seq(1101,1150),
      up_down_fast = seq(1151,1200))
  }
  else if(deg==2400){
    # row index of DEGs 2400
    row <- list(
      up_early_low = seq(1,100),   # selecting transcripts (rows) for DE over 2 timepoints
      up_early_high = seq(101,200),
      up_mid_low = seq(201,300),
      up_mid_high = seq(301,400),
      up_late_low = seq(401,500),
      up_late_high = seq(501,600),
      down_early_low = seq(601,700),
      down_early_high = seq(701,800),
      down_mid_low = seq(801,900),
      down_mid_high = seq(901,1000),
      down_late_low = seq(1001,1100),
      down_late_high = seq(1101,1200),
      up_early_slow = seq(1201,1300),   # selecting transcripts (rows) for DE over 3 timepoints
      up_mid_slow = seq(1301,1400),
      down_early_slow = seq(1401,1500),
      down_mid_slow= seq(1501,1600),
      up_mid_fast = seq(1601,1700),
      down_mid_fast= seq(1701,1800),
      up_grad_lin = seq(1801,1900),
      up_grad_cur = seq(1901,2000),
      down_grad_lin = seq(2001,2100),
      down_grad_cur = seq(2101,2200),
      up_down_slow = seq(2201,2300),
      up_down_fast = seq(2301,2400))
  }
  else {
    # row index of DEGs 3600
    row <- list(
      up_early_low = seq(1,150),   # selecting transcripts (rows) for DE over 2 timepoints
      up_early_high = seq(151,300),
      up_mid_low = seq(301,450),
      up_mid_high = seq(451,600),
      up_late_low = seq(601,750),
      up_late_high = seq(751,900),
      down_early_low = seq(901,1050),
      down_early_high = seq(1051,1200),
      down_mid_low = seq(1201,1350),
      down_mid_high = seq(1351,1500),
      down_late_low = seq(1501,1650),
      down_late_high = seq(1651,1800),
      up_early_slow = seq(1801,1950),   # selecting transcripts (rows) for DE over 3 timepoints
      up_mid_slow = seq(1951,2100),
      down_early_slow = seq(2101,2250),
      down_mid_slow= seq(2251,2400),
      up_mid_fast = seq(2401,2550),
      down_mid_fast= seq(2551,2700),
      up_grad_lin = seq(2701,2850),
      up_grad_cur = seq(2851,3000),
      down_grad_lin = seq(3001,3150),
      down_grad_cur = seq(3151,3300),
      up_down_slow = seq(3301,3450),
      up_down_fast = seq(3451,3600))
  }
  
  
  ## create replicates for more time points
  for (TP in 1:length(multTP))
  {
    # simulate replicates for control data set
    idx <- sapply(unlist(t(data)),function(x) which.min(abs(pair_f[,"mu"]-x)))
    outData1 <- matrix(sapply(idx,function(x) replicate(nRep*multTP[TP],round(mean(rnbinom(n=min(round(max(pair_f[x,"mu"]/4,1)),500),mu=pair_f[x,"mu"],size=pair_f[x,"d"]))))),ncol=nRep * nTP * multTP[TP],byrow=TRUE)
    colnames(outData1) <- paste("TP",rep(1:(nTP*multTP[TP]),each=nRep),"_",rep(seq(1:nRep),nTP*multTP[TP]),sep="")
    rownames(outData1) <- paste("gene_",seq(1,nrow(outData1)),sep="")
    # write.table(outData1,file.path(nTP_dir[h],paste("30M",nTP_dir[h],rep[j],deg,"tp_contr.txt",sep="_")))
    print(dim(outData1))
    # simulate pattern of treatment data set
    for (sim in names(row))
    {
      mu <- matrix(round(rep(rowMeans(data[row[[sim]],]),each=length(col[[sim]])) * rep(mod[[sim]],length(row[[sim]]))),ncol=length(col[[sim]]),byrow=T)
      cols <- seq(col[[sim]][1],col[[sim]][1]+ncol(mu)-1)
      data[row[[sim]],cols] <- mu
    }
    # simulate replicates for treatment data set
    idx <- sapply(unlist(t(data)),function(x) which.min(abs(pair_f[,"mu"]-x)))
    outData2 <- matrix(sapply(idx,function(x) replicate(nRep*multTP[TP],round(mean(rnbinom(n=min(round(max(pair_f[x,"mu"]/4,1)),500),mu=pair_f[x,"mu"],size=pair_f[x,"d"]))))),ncol=nRep * nTP * multTP[TP],byrow=TRUE)
    colnames(outData2) <- paste("TP",rep(1:(nTP*multTP[TP]),each=nRep),"_",rep(seq(1:nRep),nTP*multTP[TP]),sep="")
    rownames(outData2) <- paste("gene_",seq(1,nrow(outData2)),sep="")
    # write.table(outData2,file.path(nTP_dir[h],paste("30M",nTP_dir[h],rep[j],deg,"tp_treat.txt",sep="_"))) 
    colnames(outData2) <- gsub("TP", "tp", colnames(outData2))
    outData <- cbind(outData1, outData2)
    outData <-  head(outData, 18000)
    print(dim(outData))
    return(outData)
    
  }
}


normdataFunction <- function(counts){
  # print(Rep);print(TP);print(deg)
  print(head(counts,3))
  # counts_matrix<- as.matrix(counts)
  ############filter######
  # keep <- rowSums(counts>0) >= floor(0.75*ncol(counts))
  # table(keep)
  # filter_count <- counts[keep,]
  # dim(filter_count)
  #normalization#####
  express_cpm <- log2(cpm(counts)+1)
  print(head(express_cpm,3))
  # write.csv(express_cpm,file.path('norm_timedata',paste(norm_egder[2],lib_s,TP,Rep,deg,"both.csv",sep="_"))) #quantile median
  return(express_cpm)
}


