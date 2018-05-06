#!/bin/bash

## KEGG pathway enrichment analysis for DEGs or DMGs.
## Usage:
# usage="./kegg_enrich.sh deg.txt eg.txt gene_ko_map.txt kegg_pathway.txt out_prefix out_folder"
# deg.txt: differentially expressed genes, one gene per line, sorted file.
# eg.txt: expressed genes (expresion value >= ?), one gene per line, sorted file.
## gene_ko.txt: relation of genes and ko..., from blast2KO2Map-batch.pl. format:
## 			SpisGene25151	nr	154416763	tva:TVAG_110290	ko:K12472
## 			SpisGene17082	nr	156368546	nve:NEMVE_v1g245946	ko:K10284
## 			SpisGene12113	nr	156372874	nve:NEMVE_v1g245128	ko:K10274
## gene_map.txt: relation of genes and maps..., from blast2KO2Map-batch.pl. format:
## 			SpisGene25151	nr	154416763	tva:TVAG_110290	ko:K12472	map04144
## 			SpisGene3690	nr	156385000	nve:NEMVE_v1g242831	ko:K15307	ko05166
## 			SpisGene3690	nr	156385000	nve:NEMVE_v1g242831	ko:K15307	map05166
# gene_ko_map.txt: relation of genes, ko and maps. format:
# 			SpisGene10001	K05695	map04514
# 			SpisGene10001	K05695	map04520
# 			SpisGene10001	K05695	map04910
# kegg_pathway.txt: relation of map numbers and map names, from KEGG database. format:
# 			map00010	Glycolysis / Gluconeogenesis
# 			map00020	Citrate cycle (TCA cycle)
# 			map00030	Pentose phosphate pathway
# 

## variables
deg_file=$1
eg_file=$2
gene_ko_map=$3
pathway=$4
o_prefix=$5
o_folder=$6

## make out folder.
if [ ! -d $o_folder ]; then
	mkdir $o_folder
fi
cd $o_folder

# EG mapping.
eg_mapped_file=${o_prefix}".eg.ko_map.txt"
grep -w -F -f ../$eg_file ../$gene_ko_map > $eg_mapped_file
# how many EGs.
#wc -l ../$eg_file
# how many EGs were found in pathways.
eg_mapped_cnt=$(cut -f 1 $eg_mapped_file | sort | uniq | wc -l)
# 5704 of 22060 EGs mapped to pathways.

# DEG mapping.
deg_mapped_file=${o_prefix}".deg.ko_map.txt"
grep -w -F -f ../$deg_file ../$gene_ko_map > $deg_mapped_file
## how many DEGs.
#wc -l ../$deg_file
# how many DEGs were found in pathways.
deg_mapped_cnt=$(cut -f 1 $deg_mapped_file | sort | uniq | wc -l)
# 271 of 842 DEGs mapped to pathways.

# how many genes in every map.
genome_geneNum_per_map=$o_prefix".all.geneNum.txt"
cut -f 1,3 ../$gene_ko_map | sort -k2,2 -k 1,1 | uniq | cut -f 2 | uniq -c | sed -e 's/^\s\{1,\}//' -e 's/\s\{1,\}/\t/' > $genome_geneNum_per_map
# how many EGs in every map.
eg_geneNum_per_map=$o_prefix".eg.geneNum.txt"
cut -f 1,3 $eg_mapped_file | sort -k2,2 -k 1,1 | uniq | cut -f 2 | uniq -c | sed -e 's/^\s\{1,\}//' -e 's/\s\{1,\}/\t/' > $eg_geneNum_per_map
# how many DEGs in every map.
deg_geneNum_per_map=$o_prefix".deg.geneNum.txt"
cut -f 1,3 $deg_mapped_file | sort -k2,2 -k 1,1 | uniq | cut -f 2 | uniq -c | sed -e 's/^\s\{1,\}//' -e 's/\s\{1,\}/\t/' > $deg_geneNum_per_map

# merge the two geneNum files to one.
deg_geneNum_merged=$o_prefix".merged.geneNum.txt"
join -j 2 -a 1 -e 0 -t $'\t' -o 0,1.1,2.1 $eg_geneNum_per_map $deg_geneNum_per_map > $deg_geneNum_merged
# add the pathway names.
deg_fisher=$o_prefix".deg.fisher.txt"
join -j 1 -a 1 -t $'\t' -o 0,2.2,1.2,1.3 $deg_geneNum_merged ../$pathway > $deg_fisher
# adjust the format.
sed -i "s/$/\t$deg_mapped_cnt\t$eg_mapped_cnt/" $deg_fisher

## fisher's exact test.
#		Path	~Path
#DEG		n11*	n12		n1p*
#~DEG	n21		n22		n2p
#		np1*	np2		npp*
## npp=5704, n1p=271, 

## R - fisher and FDR.
fdr_file=$o_prefix".fisher.fdr.txt"
r_file=$o_prefix".fisher.fdr.RData"
Rscript ../deg_fisher.R $deg_fisher $fdr_file $r_file


