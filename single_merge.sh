# ! /bin/bash

aln1=$1
aln2=$2
outpath=$3

sourcedir=`pwd`

r1=0.1
r2=0.1
maxN=2000

s1=`mktemp`
s2=`mktemp`
smerged=`mktemp`
alnMerged=`mktemp`

len1=`grep ">" $aln1 | wc -l`
len2=`grep ">" $aln2 | wc -l`

n1=`echo $len1*$r1|bc|numlist -int`
n2=`echo $len2*$r2|bc|numlist -int`

echo $n1
echo $n2

if (( n1 > maxN )); then
	n1=$maxN
fi

if (( $n2 > $maxN )); then
	n2=$maxN
fi

echo $n1
echo $n2

sampling.py $aln1 $s1 $n1
sampling.py $aln2 $s2 $n2

java -Xmx256m -jar opal_2.1.3/Opal.jar --in $s1 --in2 $s2 --out $smerged

$sourcedir/merge_in.py $aln1 $smerged $alnMerged
$sourcedir/merge_in.py $aln2 $alnMerged $outpath

rm $s1 $s2 $smerged $alnMerged
