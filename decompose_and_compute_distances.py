from dendropy import Tree
from tree import PhylogeneticTree
from decompose_lib import decompose_by_diameter, compute_group_distance_matrix
import sys
import os

path = sys.argv[1]
max_size = int(sys.argv[2])
tree_file = sys.argv[3]
distance_file = sys.argv[4]
#min_size = int(sys.argv[4])
#max_diam = float(sys.argv[5])

t = Tree.get_from_path(path,'newick')
#T = PhylogeneticTree(t)

print('computing treeMap ... ')
#treeMap = T.decompose_tree(max_size,'centroid',decomp_strategy="brlen")#,minSize=min_size,maxDiam=max_diam)
treeMap = decompose_by_diameter(t,'centroid',max_size=max_size)
#treeNum = len(treeMap.keys())
#print(treeMap)
#print('writing tree ... ')



D = compute_group_distance_matrix(t,treeMap)

with open(distance_file,'w') as f:
    for A,B in D:
        f.write(A + " " + B + " " + str(D[(A,B)]) + "\n")

with open(tree_file,'w') as f:
    for name in treeMap:
        node = treeMap[name]
        f.write(name + " " + str(node.nleaf) + "\n")
