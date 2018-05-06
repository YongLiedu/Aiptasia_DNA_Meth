#!/bin/bash

## plot DMGs and DEGs on maps together, with up/down-regulated values.
## Usage:
# usage="./kegg_pathview_map.sh deg_logFC.txt dmg_ratio.txt mg.txt gene_ko.txt gene_ko_map.txt out_prefix out_folder"
# deg_logFC.txt: differentially expressed genes and their logFoldChange values, one gene per line, two columns, sorted file.
# dmg_ratio.txt: differentially methylated genes and their ratio values, one gene per line, two columns, sorted file.
# mg.txt: methylated genes, one gene per line, one column, sorted file.
# gene_ko.txt: all genes in genome and their ko, one gene per line, two columns, sorted file.
# gene_ko_map.txt: relation of genes, ko and maps. format:
# 			SpisGene10001	K05695	map04514
# 			SpisGene10001	K05695	map04520
# 			SpisGene10001	K05695	map04910
# 

## variables
deg_log=$1
dmg_ratio=$2
mg_gene=$3
gene_ko=$4
gene_ko_map=$5
o_prefix=$6
o_folder=$7

## make output folder.
if [ ! -d $o_folder ]; then
	mkdir $o_folder
fi
cd $o_folder

## plot DMGs and DEGs on maps together, with up/down-regulated values.
# merge DEGs and DMGs.
demg_log_ratio=$o_prefix".demg.log_ratio.txt"
join -j 1 -a 1 -a 2 -e 0 -t $'\t' -o 0,1.2,2.2 ../$deg_log ../$dmg_ratio > $demg_log_ratio
sed -i -e 's/\t0$/\t1/' $demg_log_ratio
demg_mapped=$o_prefix".demg.mapped.txt"
join -j 1 -t $'\t' ../$gene_ko_map $demg_log_ratio | cut -f 2,3,4,5 | sort | uniq > $demg_mapped
# may different genes have different expression and methylation state,
# but they match same KO, then the changed one was used, if some up some down then use 0.4.
# edit the file to add two column which show up/down expressed/methylated genes.
ko_color=$o_prefix".ko.colors.txt"
awk 'START{FS="\t";OFS="\t"}{if($3==0){expr=0}else if($3<0){expr=-1}else{expr=1} if($4==1){meth=0}else if($4>1){meth=1}else{meth=-1} if($1==lastko){if(expr==0){expr=lastexpr}else if(lastexpr!=0 && expr!=lastexpr){expr=0.4} if(meth==0){meth=lastmeth}else if(lastmeth!=0 && meth!=lastmeth){meth=0.4} }else{print lastline"\t"lastexpr"\t"lastmeth} lastline=$0;lastko=$1;lastexpr=expr;lastmeth=meth}END{print lastline"\t"lastexpr"\t"lastmeth}' $demg_mapped | sed '1d' | cut -f 1,3,4,5,6 > $ko_color
# get the corresponding map list.
maps=${o_prefix}".maps.txt"
cut -f 2 $demg_mapped | sed 's/map//' | sort | uniq > $maps

# add genes present in genome and methylated genes.
# genes present in genome
gene_present=$o_prefix".ko.present.txt"
cut -f 2 ../$gene_ko | sort | uniq | sed 's/$/\t1/' > $gene_present
# methylated genes
gene_meth=$o_prefix".ko.meth.txt"
grep -w -F -f ../$mg_gene ../$gene_ko | cut -f 2 | sort | uniq | sed 's/$/\t1/' > $gene_meth
# merge DEGs, genes present in genome, methylated genes and DMGs together.
gene_present_meth=$o_prefix".ko.present.meth.txt"
join -j 1 -a 1 -e 0 -t $'\t' -o 0,1.2,2.2 $gene_present $gene_meth > $gene_present_meth
mapplot=$o_prefix".plotmap.txt"
join -j 1 -a 2 -e 0 -t $'\t' -o 0,1.2,1.3,1.4,2.2,2.3,1.5 $ko_color $gene_present_meth > $mapplot

## R - pathview.
map_suffix=$o_prefix
r_file=$o_prefix".pathview.RData"
Rscript ../deg_pathview2.R $mapplot $maps $map_suffix $r_file


