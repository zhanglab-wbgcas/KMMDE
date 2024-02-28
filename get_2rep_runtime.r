library(reticulate)
library(dplyr)
library(BiocParallel)
library(DESeq2)
require(tidyverse)
require(poweRlaw)
require(maSigPro)
require(splineTimeR)
require(limma)
require(limorhyde)
require(data.table)
require(foreach)
require(plotROC)
library(Biobase)
library(lmms)
library(FunPat)
library(edgeR)
#library(statmod)
library(DESeq2)
library(DESeq)

setwd("/home/liukc/project/sdc/02func_algori/revise_dre/timedata_pro/r_compare_proj/ImpulseDE2/R")
all_file <- list.files()
n <- length(all_file)
fp <- paste("/home/liukc/project/sdc/02func_algori/revise_dre/timedata_pro/r_compare_proj/ImpulseDE2/R", all_file, sep = "/")
library(compiler)
for (i in 1:n){
  file = fp[i]
  print(file)
  # print(class(file))
  source(file)
}

source('/home/liukc/project/sdc/02func_algori/revise_dre/timedata_pro/02r_pro/FunPat/R/SEL.TS.AREA.R') 
use_python("/home/liukc/miniconda3/envs/sdcfunc/bin/python3")

setwd('/home/liukc/project/sdc/02func_algori/revise_dre/timedata_pro/run_times')

#修改使得在命令行可以输入参数
args <- commandArgs(trailingOnly = TRUE)

# 确保提供了所有需要的参数
if(length(args) < 4) {
  stop("Not enough arguments. Please provide TP_range, tp_range, and replicate_range.")
}

print(args)
# 解析命令行参数
TP_range <- args[1]
tp.range <- as.numeric(args[2])
replicate.range <- as.numeric(args[3])
deg.range <- as.numeric(args[4])
print(TP_range);print(tp.range);print(replicate.range)
#/home/liukc/local/R_4_2_2/bin/Rscript  get_2rep_runtime.r 4TP 4 2 240

k <- 3
metrics <- ('euclidean')
total_genenum <- 5000
fdr_thres <- 0.05

total <- 1

time_KS<-matrix(NA, nrow=length(tp.range), ncol=length(deg.range))
time_edger<-matrix(NA, nrow=length(tp.range), ncol=length(deg.range))
time_deseq<-matrix(NA, nrow=length(tp.range), ncol=length(deg.range))
time_deseq2<-matrix(NA, nrow=length(tp.range), ncol=length(deg.range))
time_masigpro<-matrix(NA, nrow=length(tp.range), ncol=length(deg.range))
time_splinetc<-matrix(NA, nrow=length(tp.range), ncol=length(deg.range))
time_funpat<-matrix(NA, nrow=length(tp.range), ncol=length(deg.range))
time_impulsede2 <-matrix(NA, nrow=length(tp.range), ncol=length(deg.range))
row_num = 0
#####Figure 2C#########
for(rep in replicate.range){
  for(tp_num in tp.range){
    row_num =row_num +1
    col_num = 1 
    
    for(deg in deg.range){
      print('-----------------')
      print(paste0("deg loop", deg))
      
      TP <- paste(tp_num,'TP',sep="")
      Rep <- paste(rep,'rep',sep="")
      DEG <- paste(deg,'deg',sep="")
      seq <- seq(0, tp_num - 1, length = tp_num)
      
      source('/home/liukc/project/sdc/02func_algori/revise_dre/timedata_pro/02r_pro/newstim-tpdata.R')
      count <- stimdataFunction(rep,tp_num,deg)
      print(dim(count))
      
      #################KD############
      print("running KS")
      start <- proc.time()
      norm_data <- normdataFunction(count)
      
      count <- count[1:total_genenum, ]
      norm_data <- norm_data[1:total_genenum, ]
      #system2("python", args = "/bin/python")
      source_python("/home/liukc/project/sdc/02func_algori/revise_dre/timedata_pro/Kmmd.py")
      kd_res <- kd_func(as.integer(tp_num),as.integer(rep),as.data.frame(norm_data),metrics)  
      
      end <- proc.time()
      KS_runtime <- (end-start)[3][[1]]      
      print(end-start)
      print(KS_runtime)
      # 
      # ############differential other tools ###################
      count.matrix <- as.matrix(count)
      group_list <- c(rep('1',(dim(count)[2])/2),rep('2',(dim(count)[2])/2))
      group <- factor(group_list)
      
      #################edger############
      print('running edgeR')
      start <- proc.time()
      
      edgeR.dgelist = DGEList(counts = count.matrix, group = group)
      edgeR.dgelist <- calcNormFactors(edgeR.dgelist, method = 'TMM') 
      edgeR.disp <- estimateCommonDisp(edgeR.dgelist)
      edgeR.disp <- estimateTagwiseDisp(edgeR.disp)
      edgeR.test <- exactTest(edgeR.disp)
      end <- proc.time()
      edger_runtime <- (end-start)[3][[1]]     
      print(end-start)
      print(edger_runtime)
      
      #################deseq############
      print('running deseq')
      start <- proc.time()
      
      DESeq.cds = newCountDataSet(countData = count.matrix, conditions = group)
      DESeq.cds = estimateSizeFactors(DESeq.cds)
      DESeq.cds = estimateDispersions(DESeq.cds, sharingMode = "maximum", method = "pooled", fitType = "local")
      DESeq.test = nbinomTest(DESeq.cds, "1", "2")
      DESeq.pvalues = DESeq.test$pval
      
      end <- proc.time()
      deseq_runtime <- (end-start)[3][[1]]      
      print(end-start)
      print(deseq_runtime)
      
      #################deseq2############
      print('running deseq2')
      start <- proc.time()
      
      Data <- data.frame(row.names = colnames(count), group = group)
      dds <- DESeqDataSetFromMatrix(countData = count,colData = Data,design = ~ group)
      dds2 <- DESeq(dds)
      tmp <- results(dds2,contrast = c("group","1","2"))
      DEG_DESeq2 <- as.data.frame(tmp)
      end <- proc.time()
      deseq2_runtime <- (end-start)[3][[1]]      
      print(end-start)
      print(deseq2_runtime)
      
      #################masigpro############
      print('running masigpro')
      start <- proc.time()
      
      timepoints <- paste(rep(seq * 2, each = rep), "h", sep = "")
      header <- paste(timepoints, rep(1:rep), sep = "-")
      print(timepoints);print(header)
      # create design matrix
      ctrl <- rep(c(1,0), each = length(header))
      mat <- cbind(Time = as.numeric(sub("h","", timepoints)),
                   Replicate = rep(1:(as.integer(tp_num) * 2), each = as.integer(rep)),
                   Control = ctrl,
                   Treatment = as.numeric(ctrl == 0))
      #print(ctrl);print(mat)
      rownames(mat) <- colnames(count)
      df = 2
      # run differential expression analysis
      NBp <- p.vector(count,design = make.design.matrix(mat, degree = 2), min.obs=0,counts = TRUE, Q = 1)
      NBt <- T.fit(NBp, step.method = "backward")  #forward both seq shuffle
      end <- proc.time()
      masigpro_runtime <- (end-start)[3][[1]]      
      print(end-start)
      print(masigpro_runtime)
      
      #############splinetc############
      print('running splineTC')
      start <- proc.time()
      
      df <-2
      design <- data.frame(row.names = colnames(count),
                           "SampleName" = colnames(count),
                           "Time" = rep(rep(seq, each = as.integer(rep)), 2),
                           "Treatment" = rep(c("control", "case"), each = as.integer(tp_num)* as.integer(rep)),
                           "Replicate" = rep(1:(as.integer(rep)), as.integer(tp_num) * 2))
      phenoData <- new("AnnotatedDataFrame", data = design)
      d <- ExpressionSet(assayData = as.matrix(count), phenoData = phenoData)
      diffExprs <- splineDiffExprs(eSetObject = d, df = df, reference = c("control", "case")[1], intercept = FALSE)
      end <- proc.time()
      splinetc_runtime <- (end-start)[3][[1]]      
      print(end-start)
      print(splinetc_runtime)
      
      #################funpat############
      print('running funpat')
      start <- proc.time()
      
      data_fp <- list("contr"=count[,1:(as.integer(tp_num)*as.integer(rep))],"treat"=count[,(as.integer(tp_num)*as.integer(rep)+1):(dim(count)[2])])
      #print(head(data_fp))
      norm <- lapply(names(data_fp),function(exp){expr=DGEList(data_fp[[exp]]);expr=calcNormFactors(expr);cpm(expr, normalized.lib.sizes=TRUE)})
      names(norm) <- names(data_fp)
      avg <- lapply(names(norm),function(exp)
        sapply(1:as.integer(tp_num),function(T)
          rowMeans(subset(norm[[exp]],select=colnames(norm[[exp]])[((T-1) * rep + 1):((T-1) * rep + rep)]))))
      names(avg) <- names(norm)
      nC <- avg[["contr"]]
      nD <- avg[["treat"]]
      dim(nD);dim(nC)
      ## create replicate combinations for treatment and control for each gene and each time point
      replic <- Reduce(rbind,lapply(1:as.integer(tp_num),function(idx) cbind(nC[,idx],nD[,idx])))
      dim(replic)
      generankout <- 'funpat_res'
      rank.res <- SEL.TS.AREA(replicates=replic,data1=nC,data2=nD,NAcontrol=1,is.interactive=F,takelog=T,allgenes=TRUE,mainDir=generankout)
      # dim(rank.res)
      end <- proc.time()
      funpat_runtime <- (end-start)[3][[1]]     
      print(end-start)
      print(funpat_runtime)
      
      #####impulsede####
      print('running impulseDE2')
      start <- proc.time()
      header <- paste(rep(((1:tp_num)-1)*3,each=rep),"h-",rep(1:rep),sep="")
      time <- ((1:tp_num)-1)*3
      type <- c("control","case")
      colnames(count) <- paste(rep(type,each=length(header)),rep(header,2),sep="_")
      data_matrix<- as.matrix(count[rowSums(count)>0,])
      print(head(data_matrix))
      # specify experimental design
      design <- data.frame("Sample"=colnames(count),
                           "Condition"=rep(type,each=tp_num*rep),
                           "Time"=rep(rep(time,each=rep),2),
                           "Batch"=rep("B_NULL",ncol(count)), 
                           row.names=colnames(count))
      
      # DEG analysis
      impulse_results <- runImpulseDE2(matCountData = data_matrix,
                                       dfAnnotation =design,
                                       boolCaseCtrl = TRUE,
                                       scaQThres = 1,
                                       boolIdentifyTransients = TRUE)
      end <- proc.time()
      
      impulsede2_runtime <- (end-start)[3][[1]]      
      print(end-start)
      print(impulsede2_runtime)

      time_KS[row_num, col_num] <- KS_runtime
      time_edger[row_num, col_num] <- edger_runtime
      time_deseq[row_num, col_num] <- deseq_runtime
      time_deseq2[row_num, col_num] <- deseq2_runtime
      time_masigpro[row_num, col_num] <- masigpro_runtime
      time_splinetc[row_num, col_num] <- splinetc_runtime
      time_funpat[row_num, col_num] <- funpat_runtime
      time_impulsede2[row_num, col_num] <- impulsede2_runtime
      
      
      print(time_KS);print(time_edger);
      print(time_deseq);print(time_deseq2);print(time_masigpro);print(time_splinetc);print(time_funpat);
      print(impulsede2_runtime)
      
      col_num = col_num +1
      
    }
  }
}

folder_path <- "runtime_res_2rep"
write.table(time_KS, file=file.path(folder_path, paste("ks_runtime", tp_num, rep, deg, "res.txt", sep="_")), sep="\t", row.names=FALSE, col.names=TRUE)
write.table(time_edger, file=file.path(folder_path, paste("edger_runtime", tp_num, rep, deg, "res.txt", sep="_")), sep="\t", row.names=FALSE, col.names=TRUE)
write.table(time_deseq, file=file.path(folder_path, paste("deseq_runtime", tp_num, rep, deg, "res.txt", sep="_")), sep="\t", row.names=FALSE, col.names=TRUE)
write.table(time_deseq2, file=file.path(folder_path, paste("deseq2_runtime", tp_num, rep, deg, "res.txt", sep="_")), sep="\t", row.names=FALSE, col.names=TRUE)
write.table(time_masigpro, file=file.path(folder_path, paste("masigpro_runtime", tp_num, rep, deg, "res.txt", sep="_")), sep="\t", row.names=FALSE, col.names=TRUE)
write.table(time_splinetc, file=file.path(folder_path, paste("splinetc_runtime", tp_num, rep, deg, "res.txt", sep="_")), sep="\t", row.names=FALSE, col.names=TRUE)
write.table(time_funpat, file=file.path(folder_path, paste("funpat_runtime", tp_num, rep, deg, "res.txt", sep="_")), sep="\t", row.names=FALSE, col.names=TRUE)


