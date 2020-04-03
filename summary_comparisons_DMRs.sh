#!/bin/sh

#Title: summary_comparisons_DMRs.sh
#Date: 20200403
#Author: Cristian Peña
#Description
#	This script counts hyper and hypo DMRs for all comaparisons according to the given context 

#Procedure
#	The script collects information from DMR results, specifally from "visual" directory, each comparison, and *0.05.txt files

#Usage
#	input directory is the directory with all DMR results
#		bash summary_comparisons_DMRs.sh [absolute_path_to_input directory] [output_directory] [context]
#		bash summary_comparisons_DMRs.sh /mnt/nfs/bioinfdata/home/NIOO/cristianp/DMRs CpG_short_summary CpG
#
#	it creates a new tsv file with headers:
#context	group1	group2	hyper	hypo
#

inputDir=$1
outputDir=$2
context=$3

mkdir ${outputDir}
CURRENT=${outputDir}/${context}
mkdir ${CURRENT}
cp ${inputDir}/${context}/visual/*/*.txt ${CURRENT}/.

FILES=$(ls ${CURRENT} | tr "\n" " ")
array=($FILES)
echo -e context"\t"group1"\t"group2"\t"hyper"\t"hypo > ${CURRENT}/${context}_comparisons_summary.tsv
for i in "${array[@]}"
do
group1=$(echo ${i} | cut -f1 -d . | awk -F '_vs_' '{print $1}' )
group2=$(echo ${i} | cut -f1 -d . | awk -F '_vs_' '{print $2}' )
hyper=$(grep hypermethylated ${CURRENT}/${i} | wc -l )
hypo=$(grep hypomethylated ${CURRENT}/${i} | wc -l )
echo -e $context"\t"$group1"\t"$group2"\t"$hyper"\t"$hypo >> ${CURRENT}/${context}_comparisons_summary.tsv
done
rm ${CURRENT}/*.txt
