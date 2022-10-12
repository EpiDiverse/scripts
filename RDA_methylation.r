# Author: Iris Sammarco
# Date: 19/09/2022
# Description: perform redundancy analysis (RDA) with methylated positions, using country and growing condition (field and garden) as predictors, and test if the interaction is significant (anova). The output files are txt files with the anova significance levels (no RDA plots are plotted here). The significance level can be useful to assess whether clusters found in a PCA analysis show statistically significant differences in your data (in this case, my data clustered according to the country of the populations)

library(vegan) # Used to run RDA
library(data.table) # Used to run fread

#Read files as dataframes:
meth = fread("CG.bed", header=T, check.names=F, sep="\t", dec=".", data.table = FALSE) # define the methylation input file here. This is a multisample bed file containing methylated positions on the row and samples on the columns. It looks like:
# chrom start end sample1 sample2
# 1 23  24  4 0
predictor = read.table("predictor.txt", header=T, sep="\t") # define the predictor input file here
#the predictor file looks like this:
#ID	Country	Condition
#FV_CZ_01_01_P0_WC0_M1_1	Czechia Field
#FV_CZ_01_03_P0_WC0_M1_1	Czechia Garden
# I have samples (ID), and Country and Condition are my predictors

## Methylation
# Remove the first 3 col (chr, start, end) and transpose the table
meth = meth[,c(4:171)] #replace 171 with the number of columns you have (number of samples + 3)
meth = t(meth)
#meth = t(na.omit(meth)) #rda can't handle NAs, remove them if you have them

## Predictor
# Make names characters (not factors)
predictor$ID = as.character(predictor$ID)
predictor$Country = as.character(predictor$Country)
predictor$Condition = as.character(predictor$Condition)
#if you had numbers, they had to be factors (as.factor)

# Confirm that methylation data and predictors are in the same order
identical(rownames(meth), predictor[,1])

# Perform Hellinger transformation on the methylation data:
meth = decostand(meth, method = "hellinger", MARGIN=1) #MARGIN=1 means rows, on rows I have samples now because I transposed the matrix

## Run RDA
#Test the country using growing condition as a Condition () (covariate):
meth_rda = rda(meth ~ Country + Condition(Condition), data = predictor, scale = F) #I remove scale since I performed Hellinger transformation
summary(meth_rda) #Constrained Proportion: variance of Y explained by X, Unconstrained Proportion: unexplained variance in Y
sink("anova_rda.txt") #with this you can save the output of the anova in a txt file, define the output directory
print(anova(meth_rda, test="Chi")) #run Anova on the result to see the significance level
sink()

#Test the interaction of your predictors:
meth_rda3 = rda(meth ~ Country * Condition + Condition(Country + Condition), data = predictor, scale = F)
summary(meth_rda3)
sink("anova_rda_interaction.txt")
print(anova(meth_rda3, test="Chi"))
sink()
