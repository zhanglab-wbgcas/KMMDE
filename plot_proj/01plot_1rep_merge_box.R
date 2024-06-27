rm(list=ls())
library(reticulate)
library(dplyr)
library(pROC)
require(plotROC)
library(ggplot2)
setwd('.\\newboxplot\\1rep_dir')

draw_and_save_plots <- function(replicate.num, tp_num) {
  meanPower <- read.csv(paste("meanPower", replicate.num, tp_num, "res.csv", sep = "_"), row.names = 1)
  vartotRejection <- read.csv(paste("vartotRejection", replicate.num, tp_num, "res.csv", sep = "_"), row.names = 1)
  meanFDR <- read.csv(paste("meanFDR", replicate.num, tp_num, "res.csv", sep = "_"), row.names = 1)
  meanf1 <- read.csv(paste("meanf1", replicate.num, tp_num, "res.csv", sep = "_"), row.names = 1)
  colnames(meanPower) <- colnames(vartotRejection) <- colnames(meanFDR) <- colnames(meanf1) <- c("KMMDE", "edgeR", "DESeq2", "maSigPro", "splineTC", 'FunPat')

  png(file = paste("01RNA_Seq_test", replicate.num, tp_num, "boxp.png", sep = "_"),
      width = 1500, height =1800, res = 300)   #02 width = 1200, height = 1000   01 width = 900, height = 1200

  par(mfrow = c(2, 2), mar = c(5, 4, 3, 2) + 0.1, mgp = c(3, 1, 0), cex = 0.8) #family = "Arial",
  
  boxplot(meanf1, col = c("red", "blue", "orange", "darkgreen", "grey", "cyan"), xlab = NULL, ylab = "F1 score", las = 2, cex.lab = 1.1, cex.axis = 1.1)
  boxplot(meanPower, col = c("red", "blue", "orange", "darkgreen", "grey", "cyan"), xlab = NULL, ylab = "Power", las = 2, cex.lab = 1.1)
  boxplot(meanFDR, col = c("red", "blue", "orange", "darkgreen", "grey", "cyan"), xlab = NULL, ylab = "False Discovery Rate", las = 2, cex.lab = 1.1)
  boxplot(vartotRejection, col = c("red", "blue", "orange", "darkgreen", "grey", "cyan"), xlab = NULL, ylab = "SD of total discoveries", las = 2, cex.lab = 1.1)
  
  mtext(paste("Time point =", tp_num), side = 3, line = -1.5, outer = TRUE, cex = 1.2) # cex = 1.2
  
  dev.off()
}

draw_and_save_plots(1, 4)
draw_and_save_plots(1, 6)
draw_and_save_plots(1, 8)
draw_and_save_plots(1, 10)

####merge####
install.packages("magick")
install.packages("gridExtra")
install.packages("png")

library(magick)
library(gridExtra)

image1 <- image_read("01RNA_Seq_test_1_4_boxp.png")
image2 <- image_read("01RNA_Seq_test_1_6_boxp.png")
image3 <- image_read("01RNA_Seq_test_1_8_boxp.png")
image4 <- image_read("01RNA_Seq_test_1_10_boxp.png")

top_row <- image_append(c(image1, image2), stack = FALSE)
bottom_row <- image_append(c(image3, image4), stack = FALSE)

combined_image <- image_append(c(top_row, bottom_row), stack = TRUE)

#save
image_write(combined_image, "combined_1repbox_image.png",format = "png")

