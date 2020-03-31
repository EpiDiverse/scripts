#!/bin/bash

#Title: merge_DMRs.sh
#Date: 20203103
#Author: Cristian Peña
#Description
#	This script merge all bed files (included in an inputDirectory) in one
#	unique bed file. In this way each context can be summarized.

#Procedure
#	The script needs an inputDirectory with all the bed files to merge, 
#	this directory can be obtained with the following code e.g. for CpG:
#		cd path_to_DMR_results_folder
#		mkdir CpG/beds
#		cp CpG/metilene/*/*0.05.bed CpG/beds/.
#
#	so "CpG/beds" is the inputDirectory

#Usage
#	first change to the directory with all DMR results
#	i.e. cd path_to_DMR_results_folder
#	bash merge_DMRs.sh inputDirectory outputDirectory
#	bash merge_DMRs.sh CpG/beds CpG_summary
#
#	it creates 3 new bed files
#		DMRs_${output}_merged.bed						contains all DMRs coordinates (just coordinates)
#		raw_intersect_all_${output}_DMRs.bed			contains all DMRs coordinates + specific coordinates per comparison
#		comp_intersect_all_${output}_DMRs.bed			contains all specific DMRs coordinates per comparison

#this script uses bedtools

inputDir=$1
output=$2
mkdir ${output}
FILES=$(ls ${inputDir} | awk -v prefix="${inputDir}/" '{print prefix $0}' | tr "\n" " ")
NAMES=$(ls ${inputDir} | cut -f1 -d . | tr "\n" " ")

cat ${inputDir}/*.bed > ${output}/DMRs.bed
bedtools sort -i ${output}/DMRs.bed > ${output}/DMRs_sorted.bed
bedtools merge -i ${output}/DMRs_sorted.bed > ${output}/DMRs_${output}_merged.bed
bedtools intersect -a ${output}/DMRs_${output}_merged.bed -b ${FILES} -wa -wb -names ${NAMES} > ${output}/raw_intersect_all_${output}_DMRs.bed
awk '{print $1 "\t" $6 "\t" $7 "\t" $4 "\t" $9}' ${output}/raw_intersect_all_${output}_DMRs.bed > ${output}/comp_intersect_all_${output}_DMRs.bed
rm ${output}/DMRs.bed
rm ${output}/DMRs_sorted.bed
