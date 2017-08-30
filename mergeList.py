#! /usr/bin/env python


import subprocess

myList = [({'A','B','C'},'1.fasta'),({'A','D','E'},'2.fasta'),({'D','F','G'},'3.fasta')]


def mergeList(L,idxName,directory):
    while len(L) > 1:
        (S,S_name) = L.pop()
        for (S1,S1_name) in L:
            if (len(S & S1)):
                L.remove((S1,S1_name))
                name = str(idxName) + ".fasta"
                subprocess.check_call(["python", "merge_in.py", directory+'/'+ S_name, directory+'/'+S1_name, directory+'/'+name])
                S = (S | S1)
                L.append((S,name))
                idxName = idxName + 1
                break
    return L[0]


print(mergeList(myList,4))
