# ! /bin/bash

aln1=$1
aln2=$2
outpath=$3

sourcedir=`pwd`


sampling.py $aln1 sample1.FASTA 50
sampling.py $aln2 sample2.FASTA 50

opal --in sample1.FASTA --in2 sample2.FASTA --out sample_merged.FASTA

$sourcedir/merge_in.py $aln1 sample_merged.FASTA aln1_merged.FASTA
$sourcedir/merge_in.py $aln2 aln1_merged.FASTA $outpath


