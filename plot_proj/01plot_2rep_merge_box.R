rm(list=ls())
library(reticulate)
library(dplyr)
library(pROC)
require(plotROC)
library(ggplot2)
setwd('E:\\差异贡献\\sdc函数\\修改\\模拟数据\\时序数据\\rna-seq_tcComp\\teststimdata\\new_timedata\\roc_plots\\box_plots\\newboxplot\\2rep_dir')

draw_and_save_plots <- function(replicate.num, tp_num) {
  meanPower <- read.csv(paste("meanPower", replicate.num, tp_num, "res.csv", sep = "_"), row.names = 1)
  vartotRejection <- read.csv(paste("vartotRejection", replicate.num, tp_num, "res.csv", sep = "_"), row.names = 1)
  meanFDR <- read.csv(paste("meanFDR", replicate.num, tp_num, "res.csv", sep = "_"), row.names = 1)
  meanf1 <- read.csv(paste("meanf1", replicate.num, tp_num, "res.csv", sep = "_"), row.names = 1)
  colnames(meanPower) <- colnames(vartotRejection) <- colnames(meanFDR) <- colnames(meanf1) <- c("KMMDE", "edgeR", "DESeq2", "maSigPro", "splineTC", 'FunPat','ImpulseDE2')
  print(head(meanPower));head(vartotRejection)
  
  png(file = paste("01RNA_Seq_test", replicate.num, tp_num, "boxp.png", sep = "_"),
      width = 1500, height =1800, res = 300)   
  
  par(mfrow = c(2, 2), mar = c(6, 4, 3, 2) + 0.1, mgp = c(3, 1, 0), cex = 0.8) #family = "Arial",
  
  boxplot(meanf1, col = c("red", "blue", "orange", "darkgreen", "grey", "cyan","purple"), xlab = NULL, ylab = "F1 score", las = 2, cex.lab = 1.1, cex.axis = 1.1)
  boxplot(meanPower, col = c("red", "blue", "orange", "darkgreen", "grey", "cyan","purple"), xlab = NULL, ylab = "Power", las = 2, cex.lab = 1.1)
  boxplot(meanFDR, col = c("red", "blue", "orange", "darkgreen", "grey", "cyan","purple"), xlab = NULL, ylab = "False Discovery Rate", las = 2, cex.lab = 1.1)
  boxplot(vartotRejection, col = c("red", "blue", "orange", "darkgreen", "grey", "cyan","purple"), xlab = NULL, ylab = "SD of total discoveries", las = 2, cex.lab = 1.1)
  
  mtext(paste("Time point =", tp_num), side = 3, line = -1.5, outer = TRUE, cex = 1.2) # cex = 1.2 可以增大标题
  
  dev.off()
}

# 调用函数
draw_and_save_plots(2, 4)
draw_and_save_plots(2, 6)
draw_and_save_plots(2, 8)
draw_and_save_plots(2, 10)  


####合并####
install.packages("magick")
install.packages("gridExtra")
install.packages("png")

library(magick)
library(gridExtra)


# 读取图片文件
image1 <- image_read("01RNA_Seq_test_2_4_boxp.png")

image2 <- image_read("01RNA_Seq_test_2_6_boxp.png")
image3 <- image_read("01RNA_Seq_test_2_8_boxp.png")
image4 <- image_read("01RNA_Seq_test_2_10_boxp.png")

# 水平拼接
top_row <- image_append(c(image1, image2), stack = FALSE)
bottom_row <- image_append(c(image3, image4), stack = FALSE)

# 垂直拼接
combined_image <- image_append(c(top_row, bottom_row), stack = TRUE)

# 保存为新的图片文件
image_write(combined_image, "combined_2rep_box_image.png",format = "png")


