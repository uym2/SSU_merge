# ! /bin/bash

dom1=$1
dom2=$2
outpath=$3

outdir=`dirname $outpath`
sourcedir=`pwd`

mkdir $outdir
cd $outdir

sampling.py $dom1 sample1.FASTA 2000
sampling.py $dom2 sample2.FASTA 2000

opal --in sample1.FASTA --in2 sample2.FASTA --out sample_merged.FASTA

$sourcedir/merge_in.py $dom1 sample_merged.FASTA dom1_merged.FASTA
$sourcedir/merge_in.py $dom2 dom1_merged.FASTA $outpath


