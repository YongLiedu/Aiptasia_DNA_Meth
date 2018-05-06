#!/bin/bash

## plot DEGs on maps, with up/down-regulated values.
## Usage:
# usage="./kegg_pathview_map.sh deg_logFC.txt gene_ko_map.txt out_prefix out_folder"
# deg_logFC.txt: differentially expressed genes and their logFoldChange values, one gene per line, two columns, sorted file.
# gene_ko_map.txt: relation of genes, ko and maps. format:
# 			SpisGene10001	K05695	map04514
# 			SpisGene10001	K05695	map04520
# 			SpisGene10001	K05695	map04910
# 

## variables
deg_log=$1
gene_ko_map=$2
o_prefix=$3
o_folder=$4

## make out folder.
if [ ! -d $o_folder ]; then
	mkdir $o_folder
fi
cd $o_folder

## plot DEGs on maps, with up/down-regulated values.
# match gene, ko and map.
ko_match=${o_prefix}".ko.matched.txt"
join -j 1 -t $'\t' ../$gene_ko_map ../$deg_log | cut -f 2,3,4 | sort | uniq > $ko_match
# may different genes have different expression and methylation state,
# but they match same KO, then the changed ones are used,
# two columns are used to show colors, if all corresponding genes are up- or down- regulated,
# then both used the same color, or else use different colors.
ko_color=${o_prefix}".ko.colors.txt"
awk 'START{FS="\t";OFS="\t";expr=0}{if($3==0){expr=0}else if($3<0){expr=-1}else{expr=1} if($1==lastko){if(expr==0){expr=lastexpr}else if(lastexpr!=0 && expr!=lastexpr){expr=0.4}  } else {if(lastexpr==1){rcol=lcol}else if(lastexpr==-1){lcol=rcol}else if(lastexpr==0){lcol=0;rcol=0} print lastko"\t"lcol"\t"rcol} lastko=$1;lastexpr=expr;if($3<0){rcol=$3}else{lcol=$3} } END{if(lastexpr==1){rcol=lcol}else if(lastexpr==-1){lcol=rcol}else if(lastexpr==0){lcol=0;rcol=0} print lastko"\t"lcol"\t"rcol}' $ko_match | sed '1d' > $ko_color
# get the corresponding map list.
maps=${o_prefix}".maps.txt"
cut -f 2 $ko_match | sed 's/map//' | sort | uniq > $maps

## R - pathview.
map_suffix=$o_prefix
r_file=$o_prefix".pathview.RData"
Rscript ../deg_pathview.R $ko_color $maps $map_suffix $r_file

