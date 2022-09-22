#Author: Iris Sammarco
#Date: 20/02/2021
#Description: performs Gene Ontology (GO) enrichment analysis using custom annotation files and plot the results

library(clusterProfiler)
library(ggplot2)

setwd("") #set working directory

#Import the custom GO files
#The files look like this:
#GO:0055114	FvH4_1g00020.t10	Biological Process:oxidation-reduction process
#GO:0055114	FvH4_1g00020.t11	Biological Process:oxidation-reduction process
#GO:0055114	FvH4_1g00020.t12	Biological Process:oxidation-reduction process
Fragaria_BP <- read.table("BP_genes2go.txt", sep = "\t", header = FALSE) #this is your custom annotation file containing all the GO terms for your species in a 3 tabs separated file: 1 column with the GO terms, 2nd with the Gene ID and 3rd with the GO term. I would reccommend to make 3 files, 1 for each GO category (BP, CC, MF), since it seems that with custom annotation it can't distinguish them and the results may be biased
Fragaria_CC <- read.table("CC_genes2go.txt", sep = "\t", header = FALSE) #same but for cellular component
Fragaria_MF <- read.table("MF_genes2go.txt", sep = "\t", header = FALSE) #molecular function

#Extract the TERM2GENE and TERM2NAME info (you'll need it later, for the enricher() function)
term2gene_BP <- Fragaria_BP[,c(1,2)] #extract first 2 columns (GO term and gene ID)
term2name_BP <- Fragaria_BP[,c(1,3)] #extract Go term and function
term2gene_CC <- Fragaria_CC[,c(1,2)] #extract first 2 columns (GO term and gene ID)
term2name_CC <- Fragaria_CC[,c(1,3)] #extract Go term and function
term2gene_MF <- Fragaria_MF[,c(1,2)] #extract first 2 columns (GO term and gene ID)
term2name_MF <- Fragaria_MF[,c(1,3)] #extract Go term and function

#Import the DMR-related genes: extract the gene ID from your annotated unionbed files. The file looks like this:
#FvH4_1g00010.t1
#FvH4_1g00010.t1
#FvH4_1g00010.t1
#FvH4_1g00010.t1
#I have 2 conditions (FF is field and GG is garden):
CpG_FF_genes <- read.table("CpG_FF_gene_t.txt", header = FALSE)
CpG_FF_genes <- as.vector(CpG_FF_genes[,1])
CpG_GG_genes <- read.table("CpG_GG_gene_t.txt", header = FALSE)
CpG_GG_genes <- as.vector(CpG_GG_genes[,1])

intersect(term2gene_BP$V2, CpG_FF_genes) #try to see if your dmr-related genes actually intersect with the universe of your species genes 

#### CpG ####
#BP
CpG_FF_genes_GO_BP <- enricher(gene= CpG_FF_genes, pvalueCutoff = 0.05, pAdjustMethod = "fdr", TERM2GENE= term2gene_BP, TERM2NAME =  term2name_BP) #run the enrichment function
head(CpG_FF_genes_GO_BP)
summary(head(CpG_FF_genes_GO_BP))
#dotplot(CpG_FF_genes_GO_BP, showCategory=30)
#barplot(CpG_FF_genes_GO_BP)
df_CpG_FF_genes_GO_BP <- CpG_FF_genes_GO_BP[,c(2,3,6,9)] # extract the columns of interest (e.g. function, adjusted p-value, gene count, gene ratio ..)
write.csv(df_CpG_FF_genes_GO_BP, 'CpG_FF_genes_GO_BP.csv', col.names = TRUE) #write them down to a file; I then checked them manually and selected for example only the top 10 categories (ordered by gene ratio)

#same for garden:
CpG_GG_genes_GO_BP <- enricher(gene= CpG_GG_genes, pvalueCutoff = 0.05, pAdjustMethod = "fdr", TERM2GENE= term2gene_BP, TERM2NAME =  term2name_BP)
head(CpG_GG_genes_GO_BP)
df_CpG_GG_genes_GO_BP <- CpG_GG_genes_GO_BP[,c(2,3,6,9)] #function, p-value, count
write.csv(df_CpG_GG_genes_GO_BP, 'CpG_GG_genes_GO_BP.csv', col.names = TRUE)

#CC
CpG_FF_genes_GO_CC <- enricher(gene= CpG_FF_genes, pvalueCutoff = 0.05, pAdjustMethod = "fdr", TERM2GENE= term2gene_CC, TERM2NAME =  term2name_CC)
head(CpG_FF_genes_GO_CC)
df_CpG_FF_genes_GO_CC <- CpG_FF_genes_GO_CC[,c(2,3,6,9)] #function, p-value, count
write.csv(df_CpG_FF_genes_GO_CC, 'CpG_FF_genes_GO_CC.csv', col.names = TRUE)

CpG_GG_genes_GO_CC <- enricher(gene= CpG_GG_genes, pvalueCutoff = 0.05, pAdjustMethod = "fdr", TERM2GENE= term2gene_CC, TERM2NAME =  term2name_CC)
head(CpG_GG_genes_GO_CC)
df_CpG_GG_genes_GO_CC <- CpG_GG_genes_GO_CC[,c(2,3,6,9)] #function, p-value, count
write.csv(df_CpG_GG_genes_GO_CC, 'CpG_GG_genes_GO_CC.csv', col.names = TRUE)

#MF
CpG_FF_genes_GO_MF <- enricher(gene= CpG_FF_genes, pvalueCutoff = 0.05, pAdjustMethod = "fdr", TERM2GENE= term2gene_MF, TERM2NAME =  term2name_MF)
head(CpG_FF_genes_GO_MF)
summary(head(CpG_FF_genes_GO_MF))
df_CpG_FF_genes_GO_MF <- CpG_FF_genes_GO_MF[,c(2,3,6,9)] #function, p-value, count
write.csv(df_CpG_FF_genes_GO_MF, 'CpG_FF_genes_GO_MF.csv', col.names = TRUE)

CpG_GG_genes_GO_MF <- enricher(gene= CpG_GG_genes, pvalueCutoff = 0.05, pAdjustMethod = "fdr", TERM2GENE= term2gene_MF, TERM2NAME =  term2name_MF)
head(CpG_GG_genes_GO_MF)
df_CpG_GG_genes_GO_MF <- CpG_GG_genes_GO_MF[,c(2,3,6,9)] #function, p-value, count
write.csv(df_CpG_GG_genes_GO_MF, 'CpG_GG_genes_GO_MF.csv', col.names = TRUE)

#plot: barplot with gene ratio (or gene count or whatever you want) and GO terms; you can plot the BP, CC and MF together if you want, or separate them in different plots. Since I had many enriched terms, I decided to plot only the BP terms. You can add other columns of course, here at the end I'm using only GeneRatio% (it's the GeneRatio*100, Condition (Field or Garden) and if you put together the 3 GO catgories you'll use the GO_category probably):
#Read the file, it looks like this:
#	Description	pvalue	Count	Condition	GeneRatio	GeneRatio%	GO_category
#GO:0007165	signal transduction	8.79E-04	148	Field	0.029892951	2.9892951	Biological Process
#GO:0006281	DNA repair	2.09E-05	80	Field	0.016158352	1.6158352	Biological Process
#GO:0007165	signal transduction	7.99E-06	150	Garden	0.031152648	3.1152648	Biological Process
#GO:0006281	DNA repair	2.78E-06	80	Garden	0.016614746	1.6614746	Biological Process

CpG_top10 = read.csv2('CpG_GO_top10_GeneRatio.csv', sep= ',', header = TRUE)
CpG_top10$Description <- as.factor(CpG_top10$Description)
CpG_top10$Condition <- factor(CpG_top10$Condition, levels = c('Field', 'Garden'), ordered = TRUE)
CpG_top10$GO_category <- as.character(CpG_top10$GO_category)
CpG_top10$GeneRatio. <- as.numeric(CpG_top10$GeneRatio.)

darkcyan_t <- adjustcolor( "darkcyan", alpha.f = 0.4)
png("CpG_GO_top10.png", width = 600, height = 600)
ggplot(data=CpG_top10, aes(x = Description)) +
  geom_bar(stat="identity", aes(y = GeneRatio., fill = Condition)) +
  xlab("") +
  ylab("Gene ratio") +
  scale_fill_manual("Condition", values = c("darkcyan", "#FF8C00")) +
  ggtitle("CpG, GO enrichment") + 
  coord_flip() + 
  theme(axis.title.x = element_text(size = 17), axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 16)) + 
  theme(legend.title = element_blank(), legend.text = element_text(size = 16)) + 
  theme(plot.title = element_text(size=20))
dev.off()