#!/bin/sh
#author: Cristian Peña
#Date: 19/06/2020
#Description: This script looks into the WGBS output folder (recomendation: first filter it by coverage) and uses a list of samples (from a text file, first column)
#             to extract the correspondent bedGraphs, then uses bedtools unionbedg to merge all samples (only the methylation column -4th column),
#             and finally, it filters every position based in a fraction of allowed NAs (eg. 0.2 = maximum 20% of NAs per position)
#             All contexts are analyzed in a for loop. If you want to do it in separate process (screens) use unionbedg_with_NAs_context.sh

#Usage: copy this script in a new empty directory where you want to create the new output. Also copy the sample text file in this directory.  Create a new screen and run the script.
#bash unionbedg_with_NAs.sh {absolute_path_to_WGBS_output_folder} {samples.txt} {output_directory_name} {allowed_NAs_fraction}

#examples:
#bash unionbedg_with_NAs.sh /mnt/nfs/bioinfdata/home/NIOO/cristianp/WGBS/filtered samplesTEST.txt output_with_NAs 0.2
#bash unionbedg_with_NAs.sh /mnt/nfs/bioinfdata/home/NIOO/cristianp/WGBS/filtered samplesDMR.tsv output_with_NAs 0.2

inputDir=$1
samples=($(cut -f1 ${2} | tr "\n" " "))
outputDir=$3
percent=$4
contexts=("CpG" "CHG" "CHH")
number_samples=${#samples[@]}
allowed_NAs=$(printf "%.0f" $(echo ${number_samples}*${percent} | bc -l))
mkdir ${outputDir}

for c in "${contexts[@]}"
do
 context=$c
 echo "processing "$number_samples" samples in "$context
 echo "NAs allowed: "$allowed_NAs" ("$4")"
 echo "process: extracting and reformatting samples"
 for i in "${samples[@]}"
 do
  id=${inputDir}/${i}/bedGraph/${i}_${context}.bedGraph
  echo ${i}_${context}.txt >> samples_${context}.txt
  cut -f-4 ${id} > ${i}_${context}.txt
 done
 input_bedtools=$(cat samples_${context}.txt | tr "\n" " ")
 echo "process: bedtools unionbedg"
 bedtools unionbedg -i ${input_bedtools} -filler NA -header -names ${samples[@]} > ${outputDir}/nonfiltered_${context}.bed
 echo "process: filtering positions with more than "$allowed_NAs" NAs"
 cat ${outputDir}/nonfiltered_${context}.bed | awk -v var="$allowed_NAs" '{count=0;for (i=4;i<=NF;i++){ if ($i=="NA"){count=count+1}}}{if (count<=var) print $0}' > ${outputDir}/unsorted_${context}.bed
 echo "process: sorting bed file"
 sort -k1,1 -k2,2n ${outputDir}/unsorted_${context}.bed >> ${outputDir}/${context}.bed
 echo "process: removing temporary files"
 rm *_${context}.txt
 rm ${outputDir}/nonfiltered_${context}.bed
 rm ${outputDir}/unsorted_${context}.bed
done


