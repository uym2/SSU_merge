#!//bin/bash

directory=$1
tree=$2


grouping=`mktemp`
nleaves=`mktemp`
distances=`mktemp`
pasta_nleaves=`mktemp`

grep ">" $directory/r*d*.fasta | sed -e "s/^.*r/r/g" -e "s/.fasta:>/ /g" > $grouping
awk '{print $1;}' $grouping | sort | uniq -c | awk '{print $2, $1;}' | sort > $pasta_nleaves

placing_and_compute_distances.py $tree $grouping 200 $nleaves $distances 

diff $nleaves $pasta_nleaves

tree_merging.py $nleaves $distances $directory
