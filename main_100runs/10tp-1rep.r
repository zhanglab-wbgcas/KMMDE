library(reticulate)
library(dplyr)
library(BiocParallel)
library(DESeq2)

require(tidyverse)
require(poweRlaw)
require(maSigPro)
require(splineTimeR)
require(limma)
require(data.table)
require(foreach)
require(plotROC)
library(Biobase)
library(lmms)
library(FunPat)
library(edgeR)
library(DESeq2)
library(DESeq)

source('/home/liukc/project/sdc/02func_algori/revise_dre/timedata_pro/02r_pro/FunPat/R/SEL.TS.AREA.R') 
use_python("/bin/python") #

setwd('/home/liukc/project/sdc/02func_algori/revise_dre/timedata_pro/100runs/plot_pro/no_deseq/100runs-1rep/')

replicate.range <- c(1)
tp.range <- c(10)
deg.range <- c(240,480,960,2400,3600)

fdr_thres <- 0.05  #0.01,0.1
r_thres<-c(0.5)
r_value = 0.5

metrics <- ('euclidean')
k <- 3
total <- 100
total_genenum <- 5000

#total <- 1  #
print(paste0("total run times", total))

#####Figure 2C#########
for(replicate.num in replicate.range){
  for(tp_num in tp.range){
    totRejection_ks<- Power_ks <- tp_ks<-tn_ks<-fn_ks<-fp_ks<-TPR_ks <-TNR_ks <- FNR_ks <- FDR_ks <-acc_ks<-pre_ks<-recall_ks<-spec_ks<-f1_ks<-matrix(NA, nrow=length(deg.range), ncol=total)
    totRejection_edger<- Power_edger <- tp_edger<-tn_edger<-fn_edger<-fp_edger<-TPR_edger <-TNR_edger <- FNR_edger <- FDR_edger <-acc_edger<-pre_edger<-recall_edger<-spec_edger<-f1_edger<-matrix(NA, nrow=length(deg.range), ncol=total)
    totRejection_deseq2<- Power_deseq2 <- tp_deseq2<-tn_deseq2<-fn_deseq2<-fp_deseq2<-TPR_deseq2 <-TNR_deseq2 <- FNR_deseq2 <- FDR_deseq2 <-acc_deseq2<-pre_deseq2<-recall_deseq2<-spec_deseq2<-f1_deseq2<-matrix(NA, nrow=length(deg.range), ncol=total)
    totRejection_2masigpro<- Power_2masigpro <- tp_2masigpro<-tn_2masigpro<-fn_2masigpro<-fp_2masigpro<-TPR_2masigpro <-TNR_2masigpro <- FNR_2masigpro <- FDR_2masigpro <-acc_2masigpro<-pre_2masigpro<-recall_2masigpro<-spec_2masigpro<-f1_2masigpro<-matrix(NA, nrow=length(deg.range), ncol=total)
    totRejection_splinetc<- Power_splinetc <- tp_splinetc<-tn_splinetc<-fn_splinetc<-fp_splinetc<-TPR_splinetc <-TNR_splinetc <- FNR_splinetc <- FDR_splinetc <-acc_splinetc<-pre_splinetc<-recall_splinetc<-spec_splinetc<-f1_splinetc<-matrix(NA, nrow=length(deg.range), ncol=total)
    totRejection_funpat<- Power_funpat <- tp_funpat<-tn_funpat<-fn_funpat<-fp_funpat<-TPR_funpat <-TNR_funpat <- FNR_funpat <- FDR_funpat <-acc_funpat<-pre_funpat<-recall_funpat<-spec_funpat<-f1_funpat<-matrix(NA, nrow=length(deg.range), ncol=total)
    
    row_num = 0
    for(deg_num in deg.range){
      print('-----------------')
      print(paste0("deg loop", deg_num))
      
      row_num = row_num + 1
      
      num=1
      
      ###run 100times
      while (num <= total){
        timestart<-Sys.time()
        print(paste0("run times", num))
        
        TP <- paste(tp_num,'TP',sep="")
        Rep <- paste(replicate.num,'rep',sep="")
        DEG <- paste(deg_num,'deg',sep="")
        seq <- seq(0, tp_num - 1, length = tp_num)
        
        source('/home/liukc/project/sdc/02func_algori/revise_dre/timedata_pro/02r_pro/newstim-tpdata.R')
        count <- stimdataFunction(replicate.num,tp_num,deg_num)
        norm_data <- normdataFunction(count)
        count <- count[1:total_genenum, ]
        norm_data <- norm_data[1:total_genenum, ]
        print(dim(count)[1])
        
        #################KS############
        print("running KS")
        #system2("python", args = "/bin/python")
        source_python("/home/liukc/project/sdc/02func_algori/revise_dre/timedata_pro/Kmmd.py")
        ks_res <- kd_func(as.integer(tp_num),as.integer(replicate.num),as.data.frame(norm_data),metrics)  #
        
        ks_res$DEG <- c(rep(0, deg_num), rep(1, (total_genenum-deg_num)))
        print(head(ks_res))
        ##
        TP <- dim(ks_res[as.numeric(ks_res$DEG)==0 & ks_res$p_value_list <= fdr_thres, ])[1]
        FN <- dim(ks_res[as.numeric(ks_res$DEG)==0 & ks_res$p_value_list > fdr_thres, ])[1]
        FP <- dim(ks_res[as.numeric(ks_res$DEG)==1 & ks_res$p_value_list <= fdr_thres, ])[1]
        TN <- dim(ks_res[as.numeric(ks_res$DEG)==1 & ks_res$p_value_list > fdr_thres, ])[1]
        
        TPR <- TP / (TP + FN)
        TNR <- TN / (TN + FP)
        FNR <- FN/(TP+FN)
        FDR <- FP/(TP+FP)
        accuracy <- (TP + TN)/(TP + TN + FP + FN)
        recall <- TP/(TP + FN)
        precision <- TP/(TP+FP)
        specificity <- TN/(TN + FP)
        f1 <- 2/((1/precision)+(1/recall))
        totRejection <- sum(ks_res$p_value_list <= fdr_thres, na.rm=TRUE)
        Power <- sum(ks_res$p_value_list[1:deg_num] <= fdr_thres, na.rm=TRUE)/(deg_num)
        print(f1);print(TP);print(FDR);print(FNR);print(totRejection);print(Power);print(accuracy)
        
        
        tp_ks[row_num, num]<-TP
        tn_ks[row_num, num]<-TN
        fn_ks[row_num, num]<-FN
        fp_ks[row_num, num]<-FP
        acc_ks[row_num, num] <- accuracy
        recall_ks[row_num, num] <- recall
        spec_ks[row_num, num] <- specificity
        pre_ks[row_num, num] <- precision
        f1_ks[row_num, num] <- f1
        FNR_ks[row_num, num] <- FNR
        FDR_ks[row_num, num] <- FDR
        TPR_ks[row_num, num] <- TPR
        TNR_ks[row_num, num] <- TNR
        totRejection_ks[row_num, num] <- totRejection
        Power_ks[row_num, num] <-Power
        
        ############differential other tools ###################
        count.matrix <- as.matrix(count)
        group_list <- c(rep('1',(dim(count)[2])/2),rep('2',(dim(count)[2])/2))
        group <- factor(group_list)
        
        print('running edgeR')
        edgeR.dgelist = DGEList(counts = count.matrix, group = group)
        # keep.exprs <- filterByExpr(dge) #
        # dge <- dge[keep.exprs,,keep.lib.sizes=FALSE]
        edgeR.dgelist <- calcNormFactors(edgeR.dgelist, method = 'TMM') # normalizes the library sizes
        edgeR.disp <- estimateCommonDisp(edgeR.dgelist)#, design, robust=T)
        edgeR.disp <- estimateTagwiseDisp(edgeR.disp)
        edgeR.test <- exactTest(edgeR.disp)
        edgeR.pvalues <- edgeR.test$table$PValue
        edgeR.adjp <- p.adjust(edgeR.pvalues, "BH")  #ֵ
        edger_res <- data.frame(
          gene = rownames(count),
          adjp = edgeR.adjp,
          pvalue = edgeR.pvalues)
        edger_res$DEG <- c(rep(0, (deg_num)), rep(1, (total_genenum-(deg_num))))
        print(head(edger_res))
        
        TP <- dim(edger_res[as.numeric(edger_res$DEG)==0 & edger_res$adjp <= fdr_thres, ])[1]
        FN <- dim(edger_res[as.numeric(edger_res$DEG)==0 & edger_res$adjp > fdr_thres, ])[1]
        FP <- dim(edger_res[as.numeric(edger_res$DEG)==1 & edger_res$adjp <= fdr_thres, ])[1]
        TN <- dim(edger_res[as.numeric(edger_res$DEG)==1 & edger_res$adjp > fdr_thres, ])[1]
        TPR <- TP / (TP + FN)
        TNR <- TN / (TN + FP)
        FNR <- FN/(TP+FN)
        FDR <- FP/(TP+FP)
        accuracy <- (TP + TN)/(TP + TN + FP + FN)
        recall <- TP/(TP + FN)
        precision <- TP/(TP+FP)
        specificity <- TN/(TN + FP)
        f1 <- 2/((1/precision)+(1/recall))
        totRejection <- sum(edger_res$adjp <= fdr_thres, na.rm=TRUE)
        Power <- sum(edger_res$adjp[1:(deg_num)] <= fdr_thres, na.rm=TRUE)/(deg_num)
        totRejection_edger[row_num, num] <- totRejection
        Power_edger[row_num, num] <-Power
        print(f1);print(TP);print(FDR);print(FNR);print(totRejection);print(Power);print(accuracy)
        
        tp_edger[row_num, num]<-TP
        tn_edger[row_num, num]<-TN
        fn_edger[row_num, num]<-FN
        fp_edger[row_num, num]<-FP
        acc_edger[row_num, num] <- accuracy
        recall_edger[row_num, num] <- recall
        spec_edger[row_num, num] <- specificity
        pre_edger[row_num, num] <- precision
        f1_edger[row_num, num] <- f1
        FNR_edger[row_num, num] <- FNR
        FDR_edger[row_num, num] <- FDR
        TPR_edger[row_num, num] <- TPR
        TNR_edger[row_num, num] <- TNR
        
        #################deseq2############
        print('running deseq2')
        Data <- data.frame(row.names = colnames(count), group = group)
        dds <- DESeqDataSetFromMatrix(countData = count,colData = Data,design = ~ group)
        dds2 <- DESeq(dds)
        tmp <- results(dds2,contrast = c("group","1","2"))
        deseq2_res <- as.data.frame(tmp)
        deseq2_res$DEG <-c(rep(0, (deg_num)), rep(1, (total_genenum-(deg_num))))
        print(head(deseq2_res))
        
        TP <- dim(deseq2_res[as.numeric(deseq2_res$DEG)==0 & deseq2_res$padj <= fdr_thres, ])[1]
        FN <- dim(deseq2_res[as.numeric(deseq2_res$DEG)==0 & deseq2_res$padj > fdr_thres, ])[1]
        FP <- dim(deseq2_res[as.numeric(deseq2_res$DEG)==1 & deseq2_res$padj <= fdr_thres, ])[1]
        TN <- dim(deseq2_res[as.numeric(deseq2_res$DEG)==1 & deseq2_res$padj > fdr_thres, ])[1]
        TPR <- TP / (TP + FN)
        TNR <- TN / (TN + FP)
        FNR <- FN/(TP+FN)
        FDR <- FP/(TP+FP)
        accuracy <- (TP + TN)/(TP + TN + FP + FN)
        recall <- TP/(TP + FN)
        precision <- TP/(TP+FP)
        specificity <- TN/(TN + FP)
        f1 <- 2/((1/precision)+(1/recall))
        totRejection <- sum(deseq2_res$padj <= fdr_thres, na.rm=TRUE)
        Power <- sum(deseq2_res$padj[1:(deg_num)] <= fdr_thres, na.rm=TRUE)/(deg_num)
        totRejection_deseq2[row_num, num] <- totRejection
        Power_deseq2[row_num, num] <-Power
        print(f1);print(TP);print(FDR);print(FNR);print(totRejection);print(Power);print(accuracy)
        
        
        tp_deseq2[row_num, num]<-TP
        tn_deseq2[row_num, num]<-TN
        fn_deseq2[row_num, num]<-FN
        fp_deseq2[row_num, num]<-FP
        acc_deseq2[row_num, num] <- accuracy
        recall_deseq2[row_num, num] <- recall
        spec_deseq2[row_num, num] <- specificity
        pre_deseq2[row_num, num] <- precision
        f1_deseq2[row_num, num] <- f1
        FNR_deseq2[row_num, num] <- FNR
        FDR_deseq2[row_num, num] <- FDR
        TPR_deseq2[row_num, num] <- TPR
        TNR_deseq2[row_num, num] <- TNR
        
        #################masigpro############
        print('running masigpro')
        
        timepoints <- paste(rep(seq * 2, each = replicate.num), "h", sep = "")
        header <- paste(timepoints, rep(1:replicate.num), sep = "-")
        print(timepoints);print(header)
        # create design matrix
        ctrl <- rep(c(1,0), each = length(header))
        mat <- cbind(Time = as.numeric(sub("h","", timepoints)),
                     Replicate = rep(1:(as.integer(tp_num) * 2), each = as.integer(replicate.num)),
                     Control = ctrl,
                     Treatment = as.numeric(ctrl == 0))
        #print(ctrl);print(mat)
        rownames(mat) <- colnames(count)
        df = 2
        # run differential expression analysis
        NBp <- p.vector(count,design = make.design.matrix(mat, degree = 2), min.obs=0,counts = TRUE, Q = 1)
        NBt <- T.fit(NBp, step.method = "backward")  #forward both seq shuffle
        head(NBt$sol)
        # p-value R-squared p.valor_TreatmentvsControl
        # masigpro.df <- NBt$sol[,c(1, 2,4)]  #p-value R-squared
        masigpro.df <- NBt$sol[,c(1, 2)]  #p-value  p.valor_TreatmentvsControl
        tail(masigpro.df);head(masigpro.df)
        dim(masigpro.df)
        ms_gene <- rownames(masigpro.df)
        Genes <- rownames(count)
        del_gene <- Genes[!Genes %in% ms_gene]  #"gene_928"  "gene_4414"
        length(del_gene)
        new_genedf <- data.frame('p-value' = rep(1, length(del_gene)) , 'R-squared' = rep(0, length(del_gene)) ,check.names=FALSE)
        rownames(new_genedf) <- del_gene
        masigpro.df <- rbind(masigpro.df, new_genedf)
        masigpro.df <- masigpro.df[order(match(rownames(masigpro.df), Genes)), ]
        masigpro.df$DEG <- c(rep(1, (deg_num)), rep(0, (total_genenum-(deg_num))))
        print(tail(masigpro.df))
        
        TP <- dim(masigpro.df[as.numeric(masigpro.df$DEG)==1 & masigpro.df[,2] > r_value, ])[1]
        FN <- dim(masigpro.df[as.numeric(masigpro.df$DEG)==1 & masigpro.df[,2] < r_value, ])[1]
        FP <- dim(masigpro.df[as.numeric(masigpro.df$DEG)==0 & masigpro.df[,2] > r_value, ])[1]
        TN <- dim(masigpro.df[as.numeric(masigpro.df$DEG)==0 & masigpro.df[,2] < r_value, ])[1]
        TPR <- TP / (TP + FN)
        TNR <- TN / (TN + FP)
        FNR <- FN/(TP+FN)
        FDR <- FP/(TP+FP)
        accuracy <- (TP + TN)/(TP + TN + FP + FN)
        recall <- TP/(TP + FN)
        precision <- TP/(TP+FP)
        specificity <- TN/(TN + FP)
        f1 <- 2/((1/precision)+(1/recall))
        totRejection <- sum(masigpro.df[,2] > r_thres, na.rm=TRUE)
        Power <- sum((masigpro.df[,2])[1:(deg_num)] > r_thres, na.rm=TRUE)/(deg_num)
        print(f1);print(TP);print(FDR);print(FNR);print(totRejection);print(Power);print(accuracy)
        
        totRejection_2masigpro[row_num, num] <- totRejection
        Power_2masigpro[row_num, num] <- Power
        tp_2masigpro[row_num, num]<-TP
        tn_2masigpro[row_num, num]<-TN
        fn_2masigpro[row_num, num]<-FN
        fp_2masigpro[row_num, num]<-FP
        acc_2masigpro[row_num, num] <- accuracy
        recall_2masigpro[row_num, num] <- recall
        spec_2masigpro[row_num, num] <- specificity
        pre_2masigpro[row_num, num] <- precision
        f1_2masigpro[row_num, num] <- f1
        FNR_2masigpro[row_num, num] <- FNR
        FDR_2masigpro[row_num, num] <- FDR
        TPR_2masigpro[row_num, num] <- TPR
        TNR_2masigpro[row_num, num] <- TNR
        
        #############splinetc############
        print('running splineTC')
        df <-2
        design <- data.frame(row.names = colnames(count),
                             "SampleName" = colnames(count),
                             "Time" = rep(rep(seq, each = as.integer(replicate.num)), 2),
                             "Treatment" = rep(c("control", "case"), each = as.integer(tp_num)* as.integer(replicate.num)),
                             "Replicate" = rep(1:(as.integer(replicate.num)), as.integer(tp_num) * 2))
        phenoData <- new("AnnotatedDataFrame", data = design)
        d <- ExpressionSet(assayData = as.matrix(count), phenoData = phenoData)
        diffExprs <- splineDiffExprs(eSetObject = d, df = df, reference = c("control", "case")[1], intercept = FALSE)
        splineTC.df <- diffExprs[,c("P.Value", "adj.P.Val")]
        splineTC.df$DEG <- c(rep(0, (deg_num)), rep(1, (total_genenum-(deg_num))))
        print(head(splineTC.df))
        
        TP <- dim(splineTC.df[as.numeric(splineTC.df$DEG)==0 & splineTC.df$adj.P.Val <= fdr_thres, ])[1]
        FN <- dim(splineTC.df[as.numeric(splineTC.df$DEG)==0 & splineTC.df$adj.P.Val > fdr_thres, ])[1]
        FP <- dim(splineTC.df[as.numeric(splineTC.df$DEG)==1 & splineTC.df$adj.P.Val <= fdr_thres, ])[1]
        TN <- dim(splineTC.df[as.numeric(splineTC.df$DEG)==1 & splineTC.df$adj.P.Val > fdr_thres, ])[1]
        TPR <- TP / (TP + FN)
        TNR <- TN / (TN + FP)
        FNR <- FN/(TP+FN)
        FDR <- FP/(TP+FP)
        accuracy <- (TP + TN)/(TP + TN + FP + FN)
        recall <- TP/(TP + FN)
        precision <- TP/(TP+FP)
        specificity <- TN/(TN + FP)
        f1 <- 2/((1/precision)+(1/recall))
        totRejection <- sum(splineTC.df$adj.P.Val <= fdr_thres, na.rm=TRUE)
        Power <- sum(splineTC.df$adj.P.Val[1:(deg_num)] <= fdr_thres, na.rm=TRUE)/(deg_num)
        totRejection_splinetc[row_num, num] <- totRejection
        Power_splinetc[row_num, num] <-Power
        print(f1);print(TP);print(FDR);print(FNR);print(totRejection);print(Power);print(accuracy)
        
        tp_splinetc[row_num, num]<-TP
        tn_splinetc[row_num, num]<-TN
        fn_splinetc[row_num, num]<-FN
        fp_splinetc[row_num, num]<-FP
        acc_splinetc[row_num, num] <- accuracy
        recall_splinetc[row_num, num] <- recall
        spec_splinetc[row_num, num] <- specificity
        pre_splinetc[row_num, num] <- precision
        f1_splinetc[row_num, num] <- f1
        FNR_splinetc[row_num, num] <- FNR
        FDR_splinetc[row_num, num] <- FDR
        TPR_splinetc[row_num, num] <- TPR
        TNR_splinetc[row_num, num] <- TNR
        
        #################funpat############
        print('running funpat')
        data_fp <- list("contr"=count[,1:(as.integer(tp_num)*as.integer(replicate.num))],"treat"=count[,(as.integer(tp_num)*as.integer(replicate.num)+1):(dim(count)[2])])
        #print(head(data_fp))
        norm <- lapply(names(data_fp),function(exp){expr=DGEList(data_fp[[exp]]);expr=calcNormFactors(expr);cpm(expr, normalized.lib.sizes=TRUE)})
        names(norm) <- names(data_fp)
        avg <- lapply(names(norm),function(exp)
          sapply(1:as.integer(tp_num),function(T)
            rowMeans(subset(norm[[exp]],select=colnames(norm[[exp]])[((T-1) * replicate.num + 1):((T-1) * replicate.num + replicate.num)]))))
        names(avg) <- names(norm)
        nC <- avg[["contr"]]
        nD <- avg[["treat"]]
        dim(nD);dim(nC)
        ## create replicate combinations for treatment and control for each gene and each time point
        replic <- Reduce(rbind,lapply(1:as.integer(tp_num),function(idx) cbind(nC[,idx],nD[,idx])))
        dim(replic)
        generankout <- 'funpat_res'
        rank.res <- SEL.TS.AREA(replicates=replic,data1=nC,data2=nD,NAcontrol=1,is.interactive=F,takelog=T,allgenes=TRUE,mainDir=generankout)
        dim(rank.res)
        funpat.df <- cbind(rank.res$element_ID,rank.res$adjusted_p_value)
        funpat.df <- as.data.frame(funpat.df)
        dim(funpat.df)
        colnames(funpat.df) <- c('gene','adjusted_p_value')  #
        head(funpat.df);tail(funpat.df);class(funpat.df)
        
        ms_gene <- funpat.df$gene
        Genes <- rownames(count)
        del_gene <- Genes[!Genes %in% ms_gene]
        length(del_gene)
        new_genedf <- data.frame('gene'=del_gene, 'adjusted_p_value' = rep(1, length(del_gene)) ,check.names=FALSE)
        funpat.df <- rbind(funpat.df, new_genedf)
        dim(funpat.df)
        funpat.df <- funpat.df[order(match(funpat.df$gene, Genes)), ]
        funpat.df$DEG <- c(rep(0, (deg_num)), rep(1, (total_genenum-(deg_num))))
        print(tail(funpat.df))
        
        TP <- dim(funpat.df[as.numeric(funpat.df$DEG)==0 & funpat.df$adjusted_p_value <= fdr_thres, ])[1]
        FN <- dim(funpat.df[as.numeric(funpat.df$DEG)==0 & funpat.df$adjusted_p_value > fdr_thres, ])[1]
        FP <- dim(funpat.df[as.numeric(funpat.df$DEG)==1 & funpat.df$adjusted_p_value <= fdr_thres, ])[1]
        TN <- dim(funpat.df[as.numeric(funpat.df$DEG)==1 & funpat.df$adjusted_p_value > fdr_thres, ])[1]
        TPR <- TP / (TP + FN)
        TNR <- TN / (TN + FP)
        FNR <- FN/(TP+FN)
        FDR <- FP/(TP+FP)
        accuracy <- (TP + TN)/(TP + TN + FP + FN)
        recall <- TP/(TP + FN)
        precision <- TP/(TP+FP)
        specificity <- TN/(TN + FP)
        f1 <- 2/((1/precision)+(1/recall))
        totRejection <- sum(funpat.df$adjusted_p_value <= fdr_thres, na.rm=TRUE)
        Power <- sum(funpat.df$adjusted_p_value[1:(deg_num)] <= fdr_thres, na.rm=TRUE)/(deg_num)
        totRejection_funpat[row_num, num] <- totRejection
        Power_funpat[row_num, num] <-Power
        print(f1);print(TP);print(FDR);print(FNR);print(totRejection);print(Power);print(accuracy)
        
        tp_funpat[row_num, num]<-TP
        tn_funpat[row_num, num]<-TN
        fn_funpat[row_num, num]<-FN
        fp_funpat[row_num, num]<-FP
        acc_funpat[row_num, num] <- accuracy
        recall_funpat[row_num, num] <- recall
        spec_funpat[row_num, num] <- specificity
        pre_funpat[row_num, num] <- precision
        f1_funpat[row_num, num] <- f1
        FNR_funpat[row_num, num] <- FNR
        FDR_funpat[row_num, num] <- FDR
        TPR_funpat[row_num, num] <- TPR
        TNR_funpat[row_num, num] <- TNR
        
        num <- num+1
        
        timeend<-Sys.time()
        runningtime<-timeend-timestart
        print(runningtime)
      }
      write.csv(FNR_ks, file=paste("FNR_25kmmde",replicate.num,tp_num,"res.csv",sep="_"))
      write.csv(FDR_ks, file=paste("FDR_25kmmde",replicate.num,tp_num,"res.csv",sep="_"))
      write.csv(f1_ks, file=paste("f1_25kmmde",replicate.num,tp_num,"res.csv",sep="_"))
      write.csv(totRejection_ks,file=paste("totRejection_25kmmde",replicate.num,tp_num,"res.csv",sep="_"))
      write.csv(Power_ks,file=paste("Power_25kmmde",replicate.num,tp_num,"res.csv",sep="_"))
      
      write.csv(FNR_edger, file=paste("FNR_edger",replicate.num,tp_num,"res.csv",sep="_"))
      write.csv(FDR_edger, file=paste("FDR_edger",replicate.num,tp_num,"res.csv",sep="_"))
      write.csv(f1_edger, file=paste("f1_edger",replicate.num,tp_num,"res.csv",sep="_"))
      write.csv(totRejection_edger,file=paste("totRejection_edger",replicate.num,tp_num,"res.csv",sep="_"))
      write.csv(Power_edger,file=paste("Power_edger",replicate.num,tp_num,"res.csv",sep="_"))
      
      write.csv(FNR_deseq2, file=paste("FNR_deseq2",replicate.num,tp_num,"res.csv",sep="_"))
      write.csv(FDR_deseq2, file=paste("FDR_deseq2",replicate.num,tp_num,"res.csv",sep="_"))
      write.csv(f1_deseq2, file=paste("f1_deseq2",replicate.num,tp_num,"res.csv",sep="_"))
      write.csv(totRejection_deseq2,file=paste("totRejection_deseq2",replicate.num,tp_num,"res.csv",sep="_"))
      write.csv(Power_deseq2,file=paste("Power_deseq2",replicate.num,tp_num,"res.csv",sep="_"))
      
      write.csv(FNR_2masigpro, file=paste("FNR_2masigpro",replicate.num ,tp_num,"res.csv",sep="_"))
      write.csv(FDR_2masigpro, file=paste("FDR_2masigpro",replicate.num ,tp_num,"res.csv",sep="_"))
      write.csv(f1_2masigpro, file=paste("f1_2masigpro",replicate.num ,tp_num,"res.csv",sep="_"))
      write.csv(totRejection_2masigpro,file=paste("totRejection_2masigpro",replicate.num ,tp_num,"res.csv",sep="_"))
      write.csv(Power_2masigpro,file=paste("Power_2masigpro",replicate.num ,tp_num,"res.csv",sep="_"))
      
      write.csv(FNR_splinetc, file=paste("FNR_splinetc",replicate.num,tp_num,"res.csv",sep="_"))
      write.csv(FDR_splinetc, file=paste("FDR_splinetc",replicate.num,tp_num,"res.csv",sep="_"))
      write.csv(f1_splinetc, file=paste("f1_splinetc",replicate.num,tp_num,"res.csv",sep="_"))
      write.csv(totRejection_splinetc,file=paste("totRejection_splinetc",replicate.num,tp_num,"res.csv",sep="_"))
      write.csv(Power_splinetc,file=paste("Power_splinetc",replicate.num,tp_num,"res.csv",sep="_"))
      # 
      write.csv(FNR_funpat, file=paste("FNR_funpat",replicate.num,tp_num,"res.csv",sep="_"))
      write.csv(FDR_funpat, file=paste("FDR_funpat",replicate.num,tp_num,"res.csv",sep="_"))
      write.csv(f1_funpat, file=paste("f1_funpat",replicate.num,tp_num,"res.csv",sep="_"))
      write.csv(totRejection_funpat,file=paste("totRejection_funpat",replicate.num,tp_num,"res.csv",sep="_"))
      write.csv(Power_funpat,file=paste("Power_funpat",replicate.num,tp_num,"res.csv",sep="_"))
    }
    
    
    ############ֵ##########
    vartotRejection_ks <- meantotRejection_ks <- meanPower_ks <- meantp_ks <- meantn_ks <-meanfn_ks <-meanfp_ks <- meanTPR_ks <- meanTNR_ks <-meanFNR_ks <- meanFDR_ks <- meanacc_ks <-  meanpre_ks <- meanrecall_ks <- meanspec_ks <- meanf1_ks <-rep(NA, length(deg.range))
    vartotRejection_edger <- meantotRejection_edger <- meanPower_edger <- meantp_edger <- meantn_edger <-meanfn_edger <-meanfp_edger <- meanTPR_edger <- meanTNR_edger <-meanFNR_edger <- meanFDR_edger <- meanacc_edger <- meanpre_edger <- meanrecall_edger <- meanspec_edger <- meanf1_edger <-rep(NA, length(deg.range))
    vartotRejection_deseq2 <- meantotRejection_deseq2 <- meanPower_deseq2 <- meantp_deseq2 <- meantn_deseq2 <-meanfn_deseq2 <-meanfp_deseq2 <- meanTPR_deseq2 <- meanTNR_deseq2 <-meanFNR_deseq2 <- meanFDR_deseq2 <- meanacc_deseq2 <- meanpre_deseq2 <- meanrecall_deseq2 <- meanspec_deseq2 <- meanf1_deseq2 <- rep(NA, length(deg.range))
    vartotRejection_2masigpro <- meantotRejection_2masigpro <- meanPower_2masigpro <- meantp_2masigpro <- meantn_2masigpro <-meanfn_2masigpro <-meanfp_2masigpro <- meanTPR_2masigpro <- meanTNR_2masigpro <-meanFNR_2masigpro <- meanFDR_2masigpro <- meanacc_2masigpro  <- meanpre_2masigpro <- meanrecall_2masigpro <- meanspec_2masigpro <- meanf1_2masigpro <-rep(NA, length(deg.range))
    vartotRejection_splinetc <- meantotRejection_splinetc <- meanPower_splinetc <- meantp_splinetc <- meantn_splinetc <-meanfn_splinetc <-meanfp_splinetc <- meanTPR_splinetc <- meanTNR_splinetc <-meanFNR_splinetc <- meanFDR_splinetc <- meanacc_splinetc  <- meanpre_splinetc <- meanrecall_splinetc <- meanspec_splinetc <- meanf1_splinetc <-rep(NA, length(deg.range))
    vartotRejection_funpat <- meantotRejection_funpat <- meanPower_funpat <- meantp_funpat <- meantn_funpat <-meanfn_funpat <-meanfp_funpat <- meanTPR_funpat <- meanTNR_funpat <-meanFNR_funpat <- meanFDR_funpat <- meanacc_funpat  <- meanpre_funpat <- meanrecall_funpat <- meanspec_funpat <- meanf1_funpat <-rep(NA, length(deg.range))
    
    
    for (i in 1:length(deg.range))
    {
      meantotRejection_ks[i] <-mean(totRejection_ks[i, ], na.rm=TRUE)
      vartotRejection_ks[i]<-sqrt(var(totRejection_ks[i, ], na.rm=TRUE))
      meanPower_ks[i]<-mean(Power_ks[i, ], na.rm=TRUE)
      meanacc_ks[i]<-mean(acc_ks[i, ], na.rm=TRUE)
      meanpre_ks[i]<-mean(pre_ks[i, ], na.rm=TRUE)
      meanrecall_ks[i]<-mean(recall_ks[i, ], na.rm=TRUE)
      meanspec_ks[i]<-mean(spec_ks[i, ], na.rm=TRUE)
      meanf1_ks[i]<-mean(f1_ks[i, ], na.rm=TRUE)
      meanFNR_ks[i]<-mean(FNR_ks[i, ], na.rm=TRUE)
      meanFDR_ks[i]<-mean(FDR_ks[i, ], na.rm=TRUE)
      meantp_ks[i]<-mean(tp_ks[i, ], na.rm=TRUE)
      meantn_ks[i]<-mean(tn_ks[i, ], na.rm=TRUE)
      meanfn_ks[i]<-mean(fn_ks[i, ], na.rm=TRUE)
      meanfp_ks[i]<-mean(fp_ks[i, ], na.rm=TRUE)
      meanTPR_ks[i]<-mean(TPR_ks[i, ], na.rm=TRUE)
      meanTNR_ks[i]<-mean(TNR_ks[i, ], na.rm=TRUE)
      
      meantotRejection_edger[i] <-mean(totRejection_edger[i, ], na.rm=TRUE)
      vartotRejection_edger[i]<-sqrt(var(totRejection_edger[i, ], na.rm=TRUE))
      meanPower_edger[i]<-mean(Power_edger[i, ], na.rm=TRUE)
      meanacc_edger[i]<-mean(acc_edger[i, ], na.rm=TRUE)
      meanpre_edger[i]<-mean(pre_edger[i, ], na.rm=TRUE)
      meanrecall_edger[i]<-mean(recall_edger[i, ], na.rm=TRUE)
      meanspec_edger[i]<-mean(spec_edger[i, ], na.rm=TRUE)
      meanf1_edger[i]<-mean(f1_edger[i, ], na.rm=TRUE)
      meanFNR_edger[i]<-mean(FNR_edger[i, ], na.rm=TRUE)
      meanFDR_edger[i]<-mean(FDR_edger[i, ], na.rm=TRUE)
      meantp_edger[i]<-mean(tp_edger[i, ], na.rm=TRUE)
      meantn_edger[i]<-mean(tn_edger[i, ], na.rm=TRUE)
      meanfn_edger[i]<-mean(fn_edger[i, ], na.rm=TRUE)
      meanfp_edger[i]<-mean(fp_edger[i, ], na.rm=TRUE)
      meanTPR_edger[i]<-mean(TPR_edger[i, ], na.rm=TRUE)
      meanTNR_edger[i]<-mean(TNR_edger[i, ], na.rm=TRUE)
      
      meantotRejection_deseq2[i] <-mean(totRejection_deseq2[i, ], na.rm=TRUE)
      vartotRejection_deseq2[i]<-sqrt(var(totRejection_deseq2[i, ], na.rm=TRUE))
      meanPower_deseq2[i]<-mean(Power_deseq2[i, ], na.rm=TRUE)
      meanacc_deseq2[i]<-mean(acc_deseq2[i, ], na.rm=TRUE)
      meanpre_deseq2[i]<-mean(pre_deseq2[i, ], na.rm=TRUE)
      meanrecall_deseq2[i]<-mean(recall_deseq2[i, ], na.rm=TRUE)
      meanspec_deseq2[i]<-mean(spec_deseq2[i, ], na.rm=TRUE)
      meanf1_deseq2[i]<-mean(f1_deseq2[i, ], na.rm=TRUE)
      meanFNR_deseq2[i]<-mean(FNR_deseq2[i, ], na.rm=TRUE)
      meanFDR_deseq2[i]<-mean(FDR_deseq2[i, ], na.rm=TRUE)
      meantp_deseq2[i]<-mean(tp_deseq2[i, ], na.rm=TRUE)
      meantn_deseq2[i]<-mean(tn_deseq2[i, ], na.rm=TRUE)
      meanfn_deseq2[i]<-mean(fn_deseq2[i, ], na.rm=TRUE)
      meanfp_deseq2[i]<-mean(fp_deseq2[i, ], na.rm=TRUE)
      meanTPR_deseq2[i]<-mean(TPR_deseq2[i, ], na.rm=TRUE)
      meanTNR_deseq2[i]<-mean(TNR_deseq2[i, ], na.rm=TRUE)
      
      meantotRejection_2masigpro[i] <-mean(totRejection_2masigpro[i, ], na.rm=TRUE)
      vartotRejection_2masigpro[i]<-sqrt(var(totRejection_2masigpro[i, ], na.rm=TRUE))
      meanPower_2masigpro[i]<-mean(Power_2masigpro[i, ], na.rm=TRUE)
      meanacc_2masigpro[i]<-mean(acc_2masigpro[i, ], na.rm=TRUE)
      meanpre_2masigpro[i]<-mean(pre_2masigpro[i, ], na.rm=TRUE)
      meanrecall_2masigpro[i]<-mean(recall_2masigpro[i, ], na.rm=TRUE)
      meanspec_2masigpro[i]<-mean(spec_2masigpro[i, ], na.rm=TRUE)
      meanf1_2masigpro[i]<-mean(f1_2masigpro[i, ], na.rm=TRUE)
      meanFNR_2masigpro[i]<-mean(FNR_2masigpro[i, ], na.rm=TRUE)
      meanFDR_2masigpro[i]<-mean(FDR_2masigpro[i, ], na.rm=TRUE)
      meantp_2masigpro[i]<-mean(tp_2masigpro[i, ], na.rm=TRUE)
      meantn_2masigpro[i]<-mean(tn_2masigpro[i, ], na.rm=TRUE)
      meanfn_2masigpro[i]<-mean(fn_2masigpro[i, ], na.rm=TRUE)
      meanfp_2masigpro[i]<-mean(fp_2masigpro[i, ], na.rm=TRUE)
      meanTPR_2masigpro[i]<-mean(TPR_2masigpro[i, ], na.rm=TRUE)
      meanTNR_2masigpro[i]<-mean(TNR_2masigpro[i, ], na.rm=TRUE)
      
      meantotRejection_splinetc[i] <-mean(totRejection_splinetc[i, ], na.rm=TRUE)
      vartotRejection_splinetc[i]<-sqrt(var(totRejection_splinetc[i, ], na.rm=TRUE))
      meanPower_splinetc[i]<-mean(Power_splinetc[i, ], na.rm=TRUE)
      meanacc_splinetc[i]<-mean(acc_splinetc[i, ], na.rm=TRUE)
      meanpre_splinetc[i]<-mean(pre_splinetc[i, ], na.rm=TRUE)
      meanrecall_splinetc[i]<-mean(recall_splinetc[i, ], na.rm=TRUE)
      meanspec_splinetc[i]<-mean(spec_splinetc[i, ], na.rm=TRUE)
      meanf1_splinetc[i]<-mean(f1_splinetc[i, ], na.rm=TRUE)
      meanFNR_splinetc[i]<-mean(FNR_splinetc[i, ], na.rm=TRUE)
      meanFDR_splinetc[i]<-mean(FDR_splinetc[i, ], na.rm=TRUE)
      meantp_splinetc[i]<-mean(tp_splinetc[i, ], na.rm=TRUE)
      meantn_splinetc[i]<-mean(tn_splinetc[i, ], na.rm=TRUE)
      meanfn_splinetc[i]<-mean(fn_splinetc[i, ], na.rm=TRUE)
      meanfp_splinetc[i]<-mean(fp_splinetc[i, ], na.rm=TRUE)
      meanTPR_splinetc[i]<-mean(TPR_splinetc[i, ], na.rm=TRUE)
      meanTNR_splinetc[i]<-mean(TNR_splinetc[i, ], na.rm=TRUE)
      
      meantotRejection_funpat[i] <-mean(totRejection_funpat[i, ], na.rm=TRUE)
      vartotRejection_funpat[i]<-sqrt(var(totRejection_funpat[i, ], na.rm=TRUE))
      meanPower_funpat[i]<-mean(Power_funpat[i, ], na.rm=TRUE)
      meanacc_funpat[i]<-mean(acc_funpat[i, ], na.rm=TRUE)
      meanpre_funpat[i]<-mean(pre_funpat[i, ], na.rm=TRUE)
      meanrecall_funpat[i]<-mean(recall_funpat[i, ], na.rm=TRUE)
      meanspec_funpat[i]<-mean(spec_funpat[i, ], na.rm=TRUE)
      meanf1_funpat[i]<-mean(f1_funpat[i, ], na.rm=TRUE)
      meanFNR_funpat[i]<-mean(FNR_funpat[i, ], na.rm=TRUE)
      meanFDR_funpat[i]<-mean(FDR_funpat[i, ], na.rm=TRUE)
      meantp_funpat[i]<-mean(tp_funpat[i, ], na.rm=TRUE)
      meantn_funpat[i]<-mean(tn_funpat[i, ], na.rm=TRUE)
      meanfn_funpat[i]<-mean(fn_funpat[i, ], na.rm=TRUE)
      meanfp_funpat[i]<-mean(fp_funpat[i, ], na.rm=TRUE)
      meanTPR_funpat[i]<-mean(TPR_funpat[i, ], na.rm=TRUE)
      meanTNR_funpat[i]<-mean(TNR_funpat[i, ], na.rm=TRUE)
      
    }
    
    #cbind
    print(cbind(meanacc_ks, meanacc_edger,  meanacc_deseq2,meanacc_2masigpro,meanacc_splinetc, meanacc_funpat))#, meanacc_lmms))
    print(cbind(meanpre_ks, meanpre_edger,  meanpre_deseq2,meanpre_2masigpro,meanpre_splinetc, meanpre_funpat))#, meanpre_lmms))
    print(cbind(meanrecall_ks, meanrecall_edger,  meanrecall_deseq2,meanrecall_2masigpro, meanrecall_splinetc, meanrecall_funpat))#, meanrecall_lmms))
    print(cbind(meanspec_ks, meanspec_edger,  meanspec_deseq2,meanspec_2masigpro,meanspec_splinetc, meanspec_funpat))#, meanspec_lmms))
    print(cbind(meanf1_ks, meanf1_edger,  meanf1_deseq2,meanf1_2masigpro,meanf1_splinetc, meanf1_funpat))#, meanf1_lmms))
    print(cbind(meanFNR_ks, meanFNR_edger, meanFNR_deseq2,meanFNR_2masigpro,meanFNR_splinetc, meanFNR_funpat))#, meanf1_lmms))
    print(cbind(meanFDR_ks, meanFDR_edger,  meanFDR_deseq2,meanFDR_2masigpro,meanFDR_splinetc, meanFDR_funpat))#, meanf1_lmms))
    print(cbind(meantp_ks, meantp_edger, meantp_deseq2,meantp_2masigpro,meantp_splinetc, meantp_funpat))#, meanf1_lmms))
    print(cbind(meantn_ks, meantn_edger, meantn_deseq2,meantn_2masigpro, meantn_splinetc, meantn_funpat))#, meanf1_lmms))
    print(cbind(meanfn_ks, meanfn_edger,  meanfn_deseq2,meanfn_2masigpro,meanfn_splinetc, meanfn_funpat))#, meanf1_lmms))
    print(cbind(meanfp_ks, meanfp_edger,  meanfp_deseq2,meanfp_2masigpro, meanfp_splinetc, meanfp_funpat))#, meanf1_lmms))
    print(cbind(meanTPR_ks, meanTPR_edger,  meanTPR_deseq2,meanTPR_2masigpro, meanTPR_splinetc, meanTPR_funpat))#, meanf1_lmms))
    print(cbind(meanTNR_ks, meanTNR_edger,  meanTNR_deseq2,meanTNR_2masigpro, meanTNR_splinetc, meanTNR_funpat))#, meanf1_lmms))
    print(cbind(meantotRejection_ks, meantotRejection_edger,  meantotRejection_deseq2,meantotRejection_2masigpro, meantotRejection_splinetc, meantotRejection_funpat))#, meanacc_lmms))
    print(cbind(meanPower_ks, meanPower_edger, meanPower_deseq2,meanPower_2masigpro, meanPower_splinetc, meanPower_funpat))#, meanpre_lmms))
    print(cbind(vartotRejection_ks, vartotRejection_edger, vartotRejection_deseq2,vartotRejection_2masigpro, vartotRejection_splinetc, vartotRejection_funpat))#, meanpre_lmms))
    
    meanacc<- cbind(meanacc_ks, meanacc_edger, meanacc_deseq2,meanacc_2masigpro,meanacc_splinetc, meanacc_funpat)#, meanacc_lmms)
    meanpre<-cbind(meanpre_ks, meanpre_edger,  meanpre_deseq2,meanpre_2masigpro, meanpre_splinetc, meanpre_funpat)#, meanpre_lmms)
    meanrecall<-cbind(meanrecall_ks, meanrecall_edger,  meanrecall_deseq2,meanrecall_2masigpro,meanrecall_splinetc, meanrecall_funpat)#, meanrecall_lmms)
    meanspec <-cbind(meanspec_ks, meanspec_edger, meanspec_deseq2,meanspec_2masigpro,meanspec_splinetc, meanspec_funpat)#, meanspec_lmms)
    meanf1 <-cbind(meanf1_ks, meanf1_edger, meanf1_deseq2, meanf1_2masigpro,  meanf1_splinetc, meanf1_funpat)#, meanf1_lmms)
    meanFNR <- cbind(meanFNR_ks, meanFNR_edger, meanFNR_deseq2,meanFNR_2masigpro,meanFNR_splinetc, meanFNR_funpat)#, meanf1_lmms))
    meanFDR <- cbind(meanFDR_ks, meanFDR_edger,  meanFDR_deseq2,meanFDR_2masigpro, meanFDR_splinetc, meanFDR_funpat)#, meanf1_lmms))
    meantp <- cbind(meantp_ks, meantp_edger, meantp_deseq2,meantp_2masigpro, meantp_splinetc, meantp_funpat)#, meanf1_lmms))
    meantn <- cbind(meantn_ks, meantn_edger,  meantn_deseq2,meantn_2masigpro,meantn_splinetc, meantn_funpat)#, meanf1_lmms))
    meanfn <- cbind(meanfn_ks, meanfn_edger,  meanfn_deseq2,meanfn_2masigpro, meanfn_splinetc, meanfn_funpat)#, meanf1_lmms))
    meanfp <- cbind(meanfp_ks, meanfp_edger,  meanfp_deseq2,meanfp_2masigpro, meanfp_splinetc, meanfp_funpat)#, meanf1_lmms))
    meanTPR <- cbind(meanTPR_ks, meanTPR_edger,  meanTPR_deseq2,meanTPR_2masigpro, meanTPR_splinetc, meanTPR_funpat)#, meanf1_lmms))
    meanTNR <- cbind(meanTNR_ks, meanTNR_edger,  meanTNR_deseq2,meanTNR_2masigpro, meanTNR_splinetc, meanTNR_funpat)#, meanf1_lmms))
    meantotRejection <- cbind(meantotRejection_ks, meantotRejection_edger,  meantotRejection_deseq2,meantotRejection_2masigpro,  meantotRejection_splinetc, meantotRejection_funpat)#, meanacc_lmms))
    meanPower <- cbind(meanPower_ks, meanPower_edger,  meanPower_deseq2,meanPower_2masigpro, meanPower_splinetc, meanPower_funpat)#, meanpre_lmms)
    vartotRejection <- cbind(vartotRejection_ks, vartotRejection_edger,  vartotRejection_deseq2,vartotRejection_2masigpro,vartotRejection_splinetc, vartotRejection_funpat)#, meanpre_lmms))
    
    #save
    write.csv(meanacc, file=paste("meanacc",replicate.num ,tp_num,"res.csv",sep="_"))
    write.csv(meanpre, file=paste("meanpre",replicate.num ,tp_num,"res.csv",sep="_"))
    write.csv(meanrecall, file=paste("meanrecall",replicate.num ,tp_num,"res.csv",sep="_"))
    write.csv(meanspec, file=paste("meanspec",replicate.num ,tp_num,"res.csv",sep="_"))
    write.csv(meanf1, file=paste("meanf1",replicate.num ,tp_num,"res.csv",sep="_"))
    write.csv(meanFNR, file=paste("meanFNR",replicate.num ,tp_num,"res.csv",sep="_"))
    write.csv(meanFDR, file=paste("meanFDR",replicate.num ,tp_num,"res.csv",sep="_"))
    write.csv(meantp, file=paste("meantp",replicate.num ,tp_num,"res.csv",sep="_"))
    write.csv(meantn,file=paste("meantn",replicate.num ,tp_num,"res.csv",sep="_"))
    write.csv(meanfn,file=paste("meanfn",replicate.num ,tp_num,"res.csv",sep="_"))    
    write.csv(meanfp,file=paste("meanfp",replicate.num ,tp_num,"res.csv",sep="_"))    
    write.csv(meanTPR,file=paste("meanTPR",replicate.num ,tp_num,"res.csv",sep="_"))    
    write.csv(meanTNR,file=paste("meanTNR",replicate.num ,tp_num,"res.csv",sep="_"))
    write.csv(meantotRejection,file=paste("meantotRejection",replicate.num ,tp_num,"res.csv",sep="_")) 
    write.csv(meanPower,file=paste("meanPower",replicate.num ,tp_num,"res.csv",sep="_"))
    write.csv(vartotRejection,file=paste("vartotRejection",replicate.num ,tp_num,"res.csv",sep="_"))
    
    colnames(meanacc) <- colnames(meanpre) <-  colnames(meanspec) <- colnames(meanf1) <- colnames(meanFNR) <- colnames(meanFDR) <- colnames(meantp) <- colnames(meantn) <- colnames(meanfn) <- colnames(meanfp) <- colnames(meanTPR) <- colnames(meanTNR) <-colnames(meantotRejection) <- colnames(meanPower) <-  colnames(vartotRejection) <- c("KMMDE", "edgeR",  "DESeq2", "maSigPro","splineTC",'FunPat')#,'lmms')
    
    pdf(file=paste("RNA_Seq_comparison",replicate.num ,tp_num,"boxp.pdf",sep="_"), width=4, height=4)
    par(mfrow=c(2, 2), mar=c(5,4,1,1)+0.2, mgp=c(3, 1, 0), cex=0.7)
    boxplot(meanf1, col=c("red","blue","orange", "darkgreen","grey","cyan"), xlab=NULL, ylab="F1 score", las=2)
    boxplot(meanPower,col=c("red","blue","orange", "darkgreen","grey","cyan"), xlab=NULL, ylab="Power", las=2)
    boxplot(meanFDR, col=c("red","blue","orange", "darkgreen","grey","cyan"), xlab=NULL, ylab="False Discovery Rate", las=2)
    boxplot(vartotRejection,col=c("red","blue","orange", "darkgreen","grey","cyan"), xlab=NULL, ylab="SD of total discoveries", las=2)
    
    dev.off()
    
    
  }
}









