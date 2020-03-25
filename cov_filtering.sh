#!/bin/sh

#This script works only for one sample
#Takes the bedGraphs from the output of WGBS pipeline (main folder) as input and creates a new folder with same structure but bedGraphs filtered by coverage
#in this way the output can directly be used by VIEWBS or DMR EpiDiverse pipelines 

#USAGE:
#bash cov_filtering.sh /pathToWGBSOutput sample_name filtered_folder_name coverage
#EXAMPLES:
#bash cov_filtering.sh /scr/episan/RP13/stress_G1/WGBS sample_1 WGBS_filtered 5
#bash cov_filtering.sh /mnt/nfs/bioinfdata/home/NIOO/cristianp/WGBS PN_CZ_01_51_R1_HC0_M1_1 WGBS_filtered 5

#Original WGBS output structure:
#Folder:    WGBS
#Subfolders:  sample_1
#               bam
#               bedGraph
#                 sample_1_CpG.bedGraph
#                 sample_1_CHG.bedGraph
#                 sample_1_CHH.bedGraph
#               stats
#               sample_1.bam
#             sample_2
#             sample_3
#             sample_n

#Filtered WGBS output structure:
#Folder:    WGBS_filtered
#Subfolders:  sample_1
#               bedGraph
#                 sample_1_CpG.bedGraph
#                 sample_1_CHG.bedGraph
#                 sample_1_CHH.bedGraph

inputpath=$1
#e.g. /scr/episan/RP13/stress_G1/WGBS
samples=$2
#e.g. PN_CZ_01_51_R1_HC0_M1_1
outputfolder=$3
#e.g. WGBS_filtered
cov=$4
#e.g. 5

mkdir ${outputfolder}

sample=${samples}
INPUT=${inputpath}/${sample}/bedGraph
mkdir ${outputfolder}/${sample}
CURRENT=${outputfolder}/${sample}/bedGraph
mkdir ${CURRENT}

head -1 ${INPUT}/${sample}_CpG.bedGraph > ${CURRENT}/header.${sample}_CpG.bedGraph
head -1 ${INPUT}/${sample}_CHG.bedGraph > ${CURRENT}/header.${sample}_CHG.bedGraph
head -1 ${INPUT}/${sample}_CHH.bedGraph > ${CURRENT}/header.${sample}_CHH.bedGraph

awk '{ if ($5+$6 > '${cov}') { print } }' ${INPUT}/${sample}_CpG.bedGraph > ${CURRENT}/tmp.${sample}_CpG.bedGraph
cat ${CURRENT}/header.${sample}_CpG.bedGraph ${CURRENT}/tmp.${sample}_CpG.bedGraph > ${CURRENT}/${sample}_CpG.bedGraph
awk '{ if ($5+$6 > '${cov}') { print } }' ${INPUT}/${sample}_CHG.bedGraph > ${CURRENT}/tmp.${sample}_CHG.bedGraph
cat ${CURRENT}/header.${sample}_CHG.bedGraph ${CURRENT}/tmp.${sample}_CHG.bedGraph > ${CURRENT}/${sample}_CHG.bedGraph
awk '{ if ($5+$6 > '${cov}') { print } }' ${INPUT}/${sample}_CHH.bedGraph > ${CURRENT}/tmp.${sample}_CHH.bedGraph
cat ${CURRENT}/header.${sample}_CHH.bedGraph ${CURRENT}/tmp.${sample}_CHH.bedGraph > ${CURRENT}/${sample}_CHH.bedGraph

rm ${CURRENT}/tmp.*
rm ${CURRENT}/header.*
