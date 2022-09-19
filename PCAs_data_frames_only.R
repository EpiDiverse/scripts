#!/usr/bin/env R
#Date: September 2022
#Author: Bárbara Díez Rodríguez
#Description: Create csv files with PCA coordinates and variance explained, using the Epidiverse multisample bed files. Requires samples.tsv file with ID and group columns for each individual.
#Usage: Rscript PCAs_data_frames_only.R input_file samples.tsv context

library(data.table)
args <- commandArgs(trailingOnly=TRUE)
data <- fread(args[1], header=TRUE, check.names=FALSE, data.table = FALSE)
samples <- read.table(args[2], header=FALSE)
context <- args[3]

data <- data[c(-1,-2,-3)]
pca <- prcomp(t(na.omit(data)))
dims <- as.data.frame(pca$x)
groups <- character()
origin <- character()


for (a in rownames(dims)){
    groups <- c(groups, as.character(samples[samples$V1==a,]$V2))
}

pcavar <- round((pca$sdev^2)/sum((pca$sdev^2)),3)*100
dims$groups <- groups
write.csv(dims, paste("",context,"PCA.csv",sep=""))
write.csv(pcavar, paste("",context,"_PCA_var.csv",sep=""))
