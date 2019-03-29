#!/bin/bash
#####-----------------------
## Reference v3 had no annotation
## To retrieve gene id overlapped by selection signals, v3 sequence was blasted against Laing et al. 2013 reference (haemonchus_contortus.PRJEB506.WBPS8.genomic.fa)
####-------------------------

hitfile=$1
blastdir=$2
reffa=reference_assembly_v3 # available at: ftp://ngs.sanger.ac.uk/production/pathogens/Haemonchus_contortus

rm ./top.fasta

###--- Fetch sequence corresponding to every interval
cut -f1,2,3 $hitfile > tmp
while read chr up dw
do
    if [ $chr == '5' ]
    then
	samtools faidx $reffa $chr"_Celeg_TT_axel":${up}-${dw} >> top.fasta
   else
	samtools faidx $reffa $chr"_Celeg_TT":${up}-${dw} >> top.fasta
    fi
done<tmp
rm tmp

mkdir $blastdir
cd $blastdir
python ~/scripts/split_ref_by_contigs.py ../top.fasta 

## blastall to retrieve coordinates
n=0
for i in `ls *.fasta`
do
    let "n+=1"
    echo  "blastall -p blastn -e 1e-50 -d haemonchus_contortus.PRJEB506.WBPS8.genomic.fa -m 8 -i $i -o $i.blastn" > blastn$n.sh
    echo "~/scripts/retrieve_gene_gtf_from_blast.sh $i.blastn" >> ./blastn$n.sh
done
chmod +x *sh
## Then submit each job
