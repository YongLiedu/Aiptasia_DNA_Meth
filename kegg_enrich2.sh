#!/bin/bash

## KEGG pathway enrichment analysis for DEGs and DMGs.
## Usage:
# usage="./kegg_enrich.sh deg.txt eg.txt dmg.txt mg.txt gene_ko_map.txt kegg_pathway.txt out_prefix out_folder"
# deg.txt: differentially expressed genes, one gene per line, sorted file.
# eg.txt: expressed genes (expresion value >= ?), one gene per line, sorted file.
# dmg.txt: differentially methylated genes, one gene per line, sorted file.
# mg.txt: methylated genes, one gene per line, sorted file.
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
dmg_file=$3
mg_file=$4
gene_ko_map=$5
pathway=$6
o_prefix=$7
o_folder=$8

## make out folder.
if [ ! -d $o_folder ]; then
	mkdir $o_folder
fi
cd $o_folder

# merge EG and MG.
emg_file=$o_prefix".em.txt"
cat ../$eg_file ../$mg_file | sort | uniq > $emg_file
emg_mapped_file=$o_prefix".emg.ko_map.txt"
grep -w -F -f $emg_file ../$gene_ko_map > $emg_mapped_file
## how many EMGs.
#wc -l spis_EMG.sorted.txt
# how many EMGs were found in pathways.
emg_mapped_cnt=$(cut -f 1 $emg_mapped_file | sort | uniq | wc -l)

# merge DEG and DMG.
demg_file=$o_prefix".dem.txt"
cat ../$deg_file ../$dmg_file | sort | uniq > $demg_file
demg_mapped_file=$o_prefix".demg.ko_map.txt"
grep -w -F -f $demg_file ../$gene_ko_map > $demg_mapped_file
## how many DEMGs.
#wc -l spis_DEMG-0.05.txt
# how many DEMGs were found in pathways.
demg_mapped_cnt=$(cut -f 1 $demg_mapped_file | sort | uniq | wc -l)

# how many genes in every map.
genome_geneNum_per_map=$o_prefix".all.geneNum.txt"
cut -f 1,3 ../$gene_ko_map | sort -k2,2 -k 1,1 | uniq | cut -f 2 | uniq -c | sed -e 's/^\s\{1,\}//' -e 's/\s\{1,\}/\t/' > $genome_geneNum_per_map
# how many EMGs in every map.
emg_geneNum_per_map=$o_prefix".emg.geneNum.txt"
cut -f 1,3 $emg_mapped_file | sort -k2,2 -k 1,1 | uniq | cut -f 2 | uniq -c | sed -e 's/^\s\{1,\}//' -e 's/\s\{1,\}/\t/' > $emg_geneNum_per_map
# how many DEMGs in every map.
demg_geneNum_per_map=$o_prefix".demg.geneNum.txt"
cut -f 1,3 $demg_mapped_file | sort -k2,2 -k 1,1 | uniq | cut -f 2 | uniq -c | sed -e 's/^\s\{1,\}//' -e 's/\s\{1,\}/\t/' > $demg_geneNum_per_map

# merge the two geneNum files to one.
demg_geneNum_merged=$o_prefix".merged.geneNum.txt"
join -j 2 -a 1 -e 0 -t $'\t' -o 0,1.1,2.1 $emg_geneNum_per_map $demg_geneNum_per_map > $demg_geneNum_merged
# add the pathway names.
demg_fisher=$o_prefix".demg.fisher.txt"
join -j 1 -a 1 -t $'\t' -o 0,2.2,1.2,1.3 $demg_geneNum_merged ../$pathway > $demg_fisher
# adjust the format.
sed -i "s/$/\t$demg_mapped_cnt\t$emg_mapped_cnt/" $demg_fisher

## fisher's exact test.
#		Path	~Path
#DEG		n11*	n12		n1p*
#~DEG	n21		n22		n2p
#		np1*	np2		npp*
## npp=5704, n1p=271, 

## R - fisher and FDR.
fdr_file=$o_prefix".fisher.fdr.txt"
r_file=$o_prefix".fisher.fdr.RData"
Rscript ../deg_fisher.R $demg_fisher $fdr_file $r_file

