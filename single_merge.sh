# ! /bin/bash

aln1=$1
aln2=$2
outpath=$3

sourcedir=`pwd`

s1=`mktemp`
s2=`mktemp`
smerged=`mktemp`
alnMerged=`mktemp`

sampling.py $aln1 $s1 50
sampling.py $aln2 $s2 50

opal --in $s1 --in2 $s2 --out $smerged

$sourcedir/merge_in.py $aln1 $smerged $alnMerged
$sourcedir/merge_in.py $aln2 $alnMerged $outpath

rm $s1 $s2 $smerged $alnMerged
