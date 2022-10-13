#Author: Iris Sammarco
#Date: 04/04/2022
#Description: calculate average methylation over regions of interest (eg a list of DMRs, of genes etc), from several files. The example here is for 2 files, but you can add as many as you want. For each samples, it outputs a bedgraph file with average methylation (0-1) per region of interest that looks like:
#chr	start	end	Meth
#Fvb1	733185	733248	0.0228136882129278
#Fvb1	3135712	3136036	0.831408775981524

library(genomation)
library(GenomicRanges)
library(methylKit)
library(plyr)

genomic_regions = gffToGRanges("DMRs.bed") #change here the input file, it can be a bed file, gff file (see the documentation of gffToGRanges)
#The genomic_regions input file looks like this (but remove the header):
#chrom	start	end
#Fvb1	12913261	12913431
#Fvb1	15132390	15132587

# Read bedgraphs
files_meth <- as.list(unlist(list.files(path="", pattern="*.bedGraph", full.names=TRUE, recursive=TRUE))) #insert path and pattern of the list of bedgraph files you want to extract average methylation from.
#The bedgraph input files look like this (but remove the header):
#chrom	start	end	Meth(%)	Methylated	Unmethylated	Coverage
#Fvb1	74	75	0	0	2	2
#Fvb1	150	151	66	2	1	3

meth <- methRead(location= files_meth, sample.id=list("sample1", "sample2"), pipeline= list(fraction=FALSE, chr.col=1, start.col=2, end.col=3, freqC.col=4, coverage.col=7, strand.col= 6), header = FALSE, sep = "\t", context = "CpG", resolution = "base", mincov = 5, treatment=c(0, 0), assembly="") #change names to sample.id, sequence context, treatment and genome assembly
meth=regionCounts(meth, genomic_regions, cov.bases=0, strand.aware = FALSE)
#create a dataframe for each sample and call them CpG1, CpG2 etc
for (i in 1:length(meth)){
nam <- paste("CpG", i, sep = "")
assign(nam, as.data.frame(meth[[i]], header = TRUE))
}
#calculate methylation for each sample and create a df with only chr start end Meth
res_CpG = lapply(list(CpG1, CpG2), function(x) cbind(chr=x$chr, start=x$start, end=x$end, Meth=x$numCs/x$coverage))
#res[[1]]
#create a dataframe for each sample
for (i in 1:length(res_CpG)){
nam <- paste("CpG_Meth", i, sep = "")
assign(nam, as.data.frame(res_CpG[[i]], header = TRUE))
}
# write files (edit path and name of the output files):
write.table(CpG_Meth1,"DMR1.bedgraph", row.names=FALSE, sep = '\t', dec= ".")
write.table(CpG_Meth2,"DMR2.bedgraph", row.names=FALSE, sep = '\t', dec= ".")