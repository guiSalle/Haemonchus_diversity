#!/bin/bash
####----- Take coordinates and retrieve associated gene names
###- Uses haemonchus_contortus.PRJEB506.WBPS9.canonical_geneset.gtf.gz annotation (Laing et al. 2013 Genome Biology paper)

blast=$1

echo $blast
head -1 $blast |cut -f2,9,10|awk '{if($3<$2){print $1"\t"$3"\t"$2}else{print $1"\t"$2"\t"$3}}'> $blast.bed
geneid=`intersectBed -a ./$blast.bed -b haemonchus_contortus.PRJEB506.WBPS9.canonical_geneset.gtf.gz -wb |head -1|cut -f12|sed 's/;/\t/g'|cut -f1`
echo $blast $geneid > $blast.txt