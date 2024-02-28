library(reticulate)
library(dplyr)
library(BiocParallel)
library(DESeq2)
require(tidyverse)
require(poweRlaw)
require(maSigPro)
require(splineTimeR)
require(limma)
#require(limorhyde)
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
library(pROC)

source('/home/liukc/project/sdc/02func_algori/revise_dre/timedata_pro/r_compare_proj/FunPat/R/SEL.TS.AREA.R') 
use_python("/bin/python") #
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

setwd('/home/liukc/project/sdc/02func_algori/revise_dre/timedata_pro/r_compare_proj/new_plots')

min_max_norm <- function(x) {
  return((x - min(x))/(max(x) - min(x)))
}
min_max_norm_df <- function(df) {
  df_norm <- as.data.frame(lapply(df,min_max_norm))
  return(df_norm)
}

get_auc <- function(x, y) {
  roc_obj <- roc(x, y)
  auc_val <- auc(roc_obj)
  return(auc_val)
}

replicate.range <-c(2)
tp.range <- c(4,6,8,10)
TP.range <- c('4TP','6TP','8TP','10TP')
deg.range <- c(240,480,960,2400,3600)
metrics <- c('tp','fp','fn','tn','accuracy','recall','precision','specificity','f1')
# total <- 100
total <- 1  #
total_genenum <- 5000  #10000

metrics <- ('euclidean')
k <- 3

#####Figure 2C#########
for(rep_num in replicate.range){
  result <- tibble(gene = NA, method = NA, timepoints = NA, deg_numbers=NA, value = NA)
  
  for(tp_num in tp.range){
    for(deg_num in deg.range){
      print('deg looping');print(deg_num)
      num <-1
      while (num<=total){
        timestart<-Sys.time()
        
        #stimdata & norm data
        # run_res_col =1;run_res_row =+ 1

        print(tp_num);print(deg_num)
        # deg <- deg_num
        # t <- tp 
        TP <- paste(tp_num,'TP',sep="")
        Rep <- paste(rep_num,'rep',sep="")
        DEG <- paste(deg_num,'deg',sep="")
        seq <- seq(0, tp_num - 1, length = tp_num)

        source('/home/liukc/project/sdc/02func_algori/revise_dre/timedata_pro/02r_pro/newstim-tpdata.R')
        count <- stimdataFunction(rep_num,tp_num,deg_num)
        norm_data <- normdataFunction(count)
  
        count <- count[1:total_genenum, ]
        norm_data <- norm_data[1:total_genenum, ]
        gene.num <- (dim(count))[1]
 
        #################KMMDE############

        print("running KMMDE")
        source_python("/home/liukc/project/sdc/02func_algori/revise_dre/timedata_pro/Kmmd.py")
        kd_res <- kd_func(as.integer(tp_num),as.integer(rep_num),as.data.frame(norm_data),metrics)  #
        
        rownames(kd_res) <- kd_res$gene
        kd_res$gene <- NULL
        kd_res <- kd_res[, 1, drop = FALSE]
        print(head(kd_res,3))
        # sorted_kd_res= kd_res[order(kd_res$position),]  #
        kd_res_norm <- min_max_norm_df(kd_res)
        print(head(kd_res,3));dim(kd_res)
        
        result <- result %>% add_row(gene = 1:(gene.num),
                                     method = "KMMDE",
                                     timepoints = tp_num,
                                     deg_numbers = deg_num,
                                     value = kd_res_norm[,1])
        
        
        
        ############differential other tools ###################
        count.matrix <- as.matrix(count)
        group_list <- c(rep('1',(dim(count)[2])/2),rep('2',(dim(count)[2])/2))
        group <- factor(group_list)
        
        #################edger############

        print('running edgeR')
        edgeR.dgelist = DGEList(counts = count.matrix, group = group)
        # keep.exprs <- filterByExpr(dge) #
        # dge <- dge[keep.exprs,,keep.lib.sizes=FALSE]
        edgeR.dgelist <- calcNormFactors(edgeR.dgelist, method = 'TMM') #normalizes the library sizes
        edgeR.disp <- estimateCommonDisp(edgeR.dgelist)#, design, robust=T)
        edgeR.disp <- estimateTagwiseDisp(edgeR.disp)
        edgeR.test <- exactTest(edgeR.disp)
        edgeR.pvalues <- edgeR.test$table$PValue
        print(head(edgeR.test$table))
        result <- result %>% add_row(gene = 1:(gene.num),
                                     method = "edgeR",
                                     timepoints = tp_num,
                                     deg_numbers = deg_num,
                                     value = edgeR.pvalues)
      
        ################deseq2############

        print('running deseq2')
        Data <- data.frame(row.names = colnames(count), group = group)
        dds <- DESeqDataSetFromMatrix(countData = count,colData = Data,design = ~ group)
        dds2 <- DESeq(dds)
        tmp <- results(dds2,contrast = c("group","1","2"))
        deseq2_res <- as.data.frame(tmp)
        print(head(deseq2_res));print(dim(deseq2_res))
        
        result <- result %>% add_row(gene = 1:(gene.num),
                                     method = "DESeq2",
                                     timepoints = tp_num,
                                     deg_numbers = deg_num,
                                     value = deseq2_res$pvalue)
        
        ################masigpro############
        print('running masigpro')
        timepoints <- paste(rep(seq * 3, each = rep), "h", sep = "")
        header <- paste(timepoints, rep(1:rep_num), sep = "-")
        
        # create design matrix
        ctrl <- rep(c(1,0), each = length(header))
        mat <- cbind(Time = as.numeric(sub("h","", timepoints)),
                     Replicate = rep(1:(as.integer(tp_num) * 2), each = as.integer(rep_num)),
                     Control = ctrl,
                     Treatment = as.numeric(ctrl == 0))
        dim(mat)
        rownames(mat) <- colnames(count)
        df = 2
        # run differential expression analysis
        NBp <- p.vector(count,design = make.design.matrix(mat, degree = df), min.obs=0,counts = TRUE, Q = 1)
        NBt <- T.fit(NBp, step.method = "backward")  #forward both seq shuffle
        head(NBt$sol)
        # p-value R-squared p.valor_TreatmentvsControl
        masigpro.df <- NBt$sol[,c(1, 2,4)]  #p-value R-squared
        # masigpro.df <- NBt$sol[,c(1, 4)]  #p-value  p.valor_TreatmentvsControl
        tail(masigpro.df);head(masigpro.df)
        dim(masigpro.df)
        ms_gene <- rownames(masigpro.df)
        Genes <- rownames(count)
        del_gene <- Genes[!Genes %in% ms_gene]  #"gene_928"  "gene_4414"
        length(del_gene)
        new_genedf <- data.frame('p-value' = rep(1, length(del_gene)) , 'R-squared' = rep(0, length(del_gene)),'p.valor_TreatmentvsControl'=rep(1, length(del_gene)) ,check.names=FALSE)
        rownames(new_genedf) <- del_gene
        masigpro.df <- rbind(masigpro.df, new_genedf)
        print(head(masigpro.df));print(dim(masigpro.df))
        
        result <- result %>% add_row(gene = 1:(gene.num),
                                     method = "maSigPro",
                                     timepoints = tp_num,
                                     deg_numbers = deg_num,
                                     value = masigpro.df[,1])
        
        # # #############splinetc############
        print('running splineTC')
        df <-2
        design <- data.frame(row.names = colnames(count),
                             "SampleName" = colnames(count),
                             "Time" = rep(rep(seq, each = as.integer(rep_num)), 2),
                             "Treatment" = rep(c("control", "case"), each = as.integer(tp_num)* as.integer(rep_num)),
                             "Replicate" = rep(1:(as.integer(rep_num)), as.integer(tp_num) * 2))
        phenoData <- new("AnnotatedDataFrame", data = design)
        d <- ExpressionSet(assayData = as.matrix(count), phenoData = phenoData)
        diffExprs <- splineDiffExprs(eSetObject = d, df = df, reference = c("control", "case")[1], intercept = FALSE)
        splineTC.df <- diffExprs[,c("P.Value", "adj.P.Val")]
        dim(splineTC.df)
        
        # sorted_genes <- (splineTC.df$gene)[order(gsub("gene_", "", splineTC.df$gene, fixed = TRUE) %>% as.integer)]
        # splineTC.df <- splineTC.df[order(match(splineTC.df$gene, sorted_genes)), ]
        # rownames(splineTC.df) <- splineTC.df$gene
        # splineTC.df$gene <- NULL
        print(class(splineTC.df$P.Value))
        result <- result %>% add_row(gene = 1:(gene.num),
                                     method = "splineTC",
                                     timepoints = tp_num,
                                     deg_numbers = deg_num,
                                     value = splineTC.df$P.Value)
        print(result)
        
        #################funpat############
        print('running funpat')
        dim(count)
        data_fp <- list("contr"=count[,1:(tp_num*rep_num)],"treat"=count[,(tp_num*rep_num):(dim(count)[2])])
        norm <- lapply(names(data_fp),function(exp){expr=DGEList(data_fp[[exp]]);expr=calcNormFactors(expr);cpm(expr, normalized.lib.sizes=TRUE)})
        names(norm) <- names(data_fp)
        avg <- lapply(names(norm),function(exp)
          sapply(1:tp_num,function(T)
            rowMeans(subset(norm[[exp]],select=colnames(norm[[exp]])[((T-1) * rep_num + 1):((T-1) * rep_num + rep_num)]))))
        names(avg) <- names(norm)
        nC <- avg[["contr"]]
        nD <- avg[["treat"]]
        dim(nD);dim(nC)
        ## create replicate combinations for treatment and control for each gene and each time point
        replic <- Reduce(rbind,lapply(1:tp_num,function(idx) cbind(nC[,idx],nD[,idx])))
        dim(replic)
        generankout <- 'funpat_res'
        rank.res <- SEL.TS.AREA(replicates=replic,data1=nC,data2=nD,NAcontrol=1,is.interactive=F,takelog=T,allgenes=TRUE,mainDir=generankout)
        dim(rank.res)
        funpat.df <- cbind(rank.res$element_ID,rank.res$p_value)
        funpat.df <- as.data.frame(funpat.df)
        # dim(funpat.df)
        colnames(funpat.df) <- c('gene','p_value')  #
        # head(funpat.df);tail(funpat.df);class(funpat.df)
        
        ms_gene <- funpat.df$gene
        Genes <- rownames(count)
        del_gene <- Genes[!Genes %in% ms_gene]
        length(del_gene)
        new_genedf <- data.frame('gene'=del_gene, 'p_value' = rep(1, length(del_gene)) ,check.names=FALSE)
        funpat.df <- rbind(funpat.df, new_genedf)

        sorted_genes <- (funpat.df$gene)[order(gsub("gene_", "", funpat.df$gene, fixed = TRUE) %>% as.integer)]
        
        funpat.df <- funpat.df[order(match(funpat.df$gene, sorted_genes)), ]
        
        rownames(funpat.df) <- funpat.df$gene
        
        funpat.df$gene <- NULL
        

        print(head(funpat.df));print(dim(funpat.df))
        
        print(class(gene.num));print(class(tp_num));print(class(deg_num));print(class(funpat.df$p_value))
       
        result <- result %>% add_row(gene = 1:(gene.num),
                                     method = "Funpat",
                                     timepoints = tp_num,
                                     deg_numbers = deg_num,
                                     value =  as.numeric(funpat.df$p_value))
        ####running impulseDE2####
	print('running impulseDE2')
        header <- paste(rep(((1:tp_num)-1)*3,each=as.integer(rep_num)),"h-",rep(1:as.integer(rep_num)),sep="")
        time <- ((1:tp_num)-1)*3
        type <- c("control","case")
        colnames(count) <- paste(rep(type,each=length(header)),rep(header,2),sep="_")
        data_matrix<- as.matrix(count[rowSums(count)>0,])
        print(head(data_matrix))
        # specify experimental design
        design <- data.frame("Sample"=colnames(count),
                             "Condition"=rep(type,each=tp_num*rep_num),
                             "Time"=rep(rep(time,each=rep_num),2),
                             "Batch"=rep("B_NULL",ncol(count)), 
                             row.names=colnames(count))
        
        # DEG analysis
        impulse_results <- runImpulseDE2(matCountData = data_matrix,
                                         dfAnnotation =design,
                                         boolCaseCtrl = TRUE,
                                         scaQThres = 1,
                                         boolIdentifyTransients = TRUE)
        
        impulsede2_res <- impulse_results$dfImpulseDE2Results[,c(1,2)]
        
        ms_gene <- impulsede2_res$Gene
        Genes <- rownames(count)
        del_gene <- Genes[!Genes %in% ms_gene]
        length(del_gene)
        new_genedf <- data.frame('Gene'=del_gene, 'p' = rep(1, length(del_gene)) ,check.names=FALSE)
        impulsede2_res <- rbind(impulsede2_res, new_genedf)
        print(dim(impulsede2_res))
        impulsede2_res <- impulsede2_res[order(match(impulsede2_res$Gene, Genes)), ]
        impulsede2_res$DEG <- c(rep(0, (deg_num)), rep(1, (total_genenum-(deg_num))))  
        result <- result %>% add_row(gene = 1:(gene.num),
                                     method = "ImpulseDE2",
                                     timepoints = tp_num,
                                     deg_numbers = deg_num,
                                     value =  as.numeric(impulsede2_res$p))                           
                                                                  
        print('One loop complete')
        num =num +1
        print(tail(result))
      }
      
    print('One DEG complete')
    }
  }
  result %>% filter(!is.na(value))
  write.csv(result, file=paste("newrocplot_tibble_p",rep_num,"res.csv",sep="_"), row.names = FALSE)  
      

  
}

