#!/bin/bash

#Author: Iris Sammarco
#Date: 16/03/21
#Description: create a tsv file containing the number of DMRs overlapping with different genomic regions. The input file is a DMR bed file annotated with gene and TE annotations (you can annotate it with bedtools intersect, https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html). It works with Fragaria vesca's annotation v4.0.a2 (https://www.rosaceae.org/species/fragaria_vesca/genome_v4.0.a2) (Jung et al. 2019), you may need to change the name of the genomic features if you use other annotations.
#Usage: bash annotation_DMR_summary.bash gene_TE_file context comparisons

gene_TE_file=$1 #DMR bed file annotated with genes and TEs
context=$2 #CpG, CHG, CHH
comparisons=$3 #ID of the comparison (e.g. Control_Drought)

tot_dmrs=$(cat $1 | wc -l)

echo -e "Context""\t"$context"\n""tot_annotation_matches = "$tot_dmrs"\n""Region""\t""Number_DMRs" > ${context}_intersect_summary_${comparisons}.tsv

gene=$(cut -f1,8 "$1" | awk 'BEGIN{FS="\t"} {if($2 ~ /gene/) print $1, $2'} | wc -l)
promoter=$(cut -f1,8 "$1" | awk 'BEGIN{FS="\t"} {if($2 ~ /promoter/) print $1, $2'} | wc -l)
CDS=$(cut -f1,8 "$1" | awk 'BEGIN{FS="\t"} {if($2 ~ /CDS/) print $1, $2'} | wc -l)
mRNA=$(cut -f1,8 "$1" | awk 'BEGIN{FS="\t"} {if($2 ~ /mRNA/) print $1, $2'} | wc -l)
exon=$(cut -f1,8 "$1" | awk 'BEGIN{FS="\t"} {if($2 ~ /exon/) print $1, $2'} | wc -l)
three_prime_UTR=$(cut -f1,8 "$1" | awk 'BEGIN{FS="\t"} {if($2 ~ /three_prime_UTR/) print $1, $2'} | wc -l)
five_prime_UTR=$(cut -f1,8 "$1" | awk 'BEGIN{FS="\t"} {if($2 ~ /five_prime_UTR/) print $1, $2'} | wc -l)
TE=$(cut -f1,8 "$1" | awk 'BEGIN{FS="\t"} {if($2 ~ /transposable_element/) print $1, $2'} | wc -l)

echo -e "gene""\t"$gene"\n""promoter""\t"$promoter"\n""CDS""\t"$CDS"\n""mRNA""\t"$mRNA"\n""exon""\t"$exon"\n""three_prime_UTR""\t"$three_prime_UTR"\n""five_prime_UTR""\t"$five_prime_UTR"\n""TE""\t"$TE"\n" >> ${context}_intersect_summary_${comparisons}.tsv