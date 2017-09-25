#! /usr/bin/env python

from min_spanning_tree import min_spanning_tree
from sys import argv
import subprocess
import numpy as np


nleaves=argv[1]
distances=argv[2]
directory=argv[3]

with open (nleaves,'r') as f:
    leaves=f.read()

leaves=leaves.split(' ')
leave=[]
for item in leaves:
    leave.extend(item.split('\n'))
leave_set=[]
for x in range(len(leave)):
    if x%2==0:
        leave_set.append(leave[x])
del leave_set[-1]
mapping={}
idx=0
for c in leave_set:
    mapping[c]=idx
    idx=idx+1



with open(distances,'r') as f:
    distance=f.read()

distance=distance.split('\n')

distance_vect=[]
for item in distance:
    distance_vect.append(item.split(' '))
del distance_vect[-1]

distance_matrix=[[0 for c in range(len(leave_set))] for d in range(len(leave_set))]
for [c1,c2,c3] in distance_vect:
    distance_matrix[mapping[c1]][mapping[c2]]=float(c3)
    distance_matrix[mapping[c2]][mapping[c1]]=float(c3)
distance_vector=[]
for c in distance_matrix:
    distance_vector.extend(c)
#print(distance_vector)


edge_set,num=min_spanning_tree(leave_set,distance_vector)
print('finish min spanning tree')
'''
with open(edgefile,'r') as f:
	edges=f.read()
f.close()

edges=edges.split('\r\n')
edge_set=[]
for item in edges:
	edge_set.extend([item.split(' ')])
del edge_set[-1]
'''
print(edge_set)

index=1
mylist=[]
for item in edge_set:
	label1=directory+'/'+item[0]+'.fasta'
	label2=directory+'/'+item[1]+'.fasta'
	aln_merge=directory+'/'+str(index)+'.fasta'
	subprocess.check_call(['java','-jar','/calab_data/mirarab/home/ziyang96/Tools/opal_2.1.3/Opal.jar','--in',label1,'--in2',label2,'--out',aln_merge])
	#print(label1,label2,aln_merge)
	mylist.append(({item[0],item[1]},str(index)+'.fasta'))
	index=index+1

print('finish opal')

def mergeList(L,idxName,directory):
    while len(L) > 1:
        (S,S_name) = L.pop()
        found = False
        for (S1,S1_name) in L:
            if (len(S & S1)):
                found = True
                L.remove((S1,S1_name))
                name = str(idxName) + ".fasta"
                subprocess.check_call(["python", "merge_in.py", directory+'/'+ S_name, directory+'/'+S1_name, directory+'/'+name])
                S = (S | S1)
                L.append((S,name))
                idxName = idxName + 1
                print(idxName)
                break
        if not found:
                print("Could not find the pair for this set ")
                print(S)
                print(S_name)

    return L






merge=mergeList(mylist, index, directory)

