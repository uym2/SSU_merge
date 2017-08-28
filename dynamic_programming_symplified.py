#! /usr/bin/env python

from sequence_lib import gap_propagate, read_fasta, write_fasta
from tri_merge_matrix import tri_matrix
from sys import argv
import numpy as np
import subprocess

aln1 = argv[1]
aln2 = argv[2]
aln3 = argv[3]
newAln_file = argv[4]


struct1=tri_matrix(aln1,aln2,aln3)
print('finish constructing matrix, start dunamic programming')


def new_map(mapping,input):
    try:
        mapping[input]
    except:
        return 0
    else:
        return mapping[input]

def max_of_cube(mapping,submap,input):
    L=[0 for x in range(7)]
    L[0]=new_map(submap,(input[0]-1,input[1]-1,input[2]-1))+new_map(mapping,(input[0],input[1],input[2]))
    L[1]=new_map(submap,(input[0],input[1]-1,input[2]-1))+new_map(mapping,(-1,input[1],input[2]))
    L[2]=new_map(submap,(input[0]-1,input[1],input[2]-1))+new_map(mapping,(input[0],-1,input[2]))
    L[3]=new_map(submap,(input[0]-1,input[1]-1,input[2]))+new_map(mapping,(input[0],input[1],-1))
    L[4]=new_map(submap,(input[0],input[1],input[2]-1))+new_map(mapping,(-1,-1,input[2]))
    L[5]=new_map(submap,(input[0],input[1]-1,input[2]))+new_map(mapping,(-1,input[1],-1))
    L[6]=new_map(submap,(input[0]-1,input[1],input[2]))+new_map(mapping,(input[0],-1,-1))
    if max(L)==L[0]:
        return max(L),[input[0]-1,input[1]-1,input[2]-1]
    elif max(L)==L[1]:
        return max(L),[input[0],input[1]-1,input[2]-1]
    elif max(L)==L[2]:
        return max(L),[input[0]-1,input[1],input[2]-1]
    elif max(L)==L[3]:
        return max(L),[input[0]-1,input[1]-1,input[2]]
    elif max(L)==L[4]:
        return max(L),[input[0],input[1],input[2]-1]
    elif max(L)==L[5]:
        return max(L),[input[0],input[1]-1,input[2]]
    elif max(L)==L[6]:
        return max(L),[input[0]-1,input[1],input[2]]


def dynamic_program(start, end, mapping):
    submap={}
    submap[(start[0],start[1],start[2])]=0
    trace_back={}
    '''
    for i  in range(1,end[0]-start[0]+1):
        submap[(start[0]+i,start[1],start[2])]=submap[(start[0]+i-1,start[1],start[2])]+new_map(mapping,(start[0]+i,0,0))
    for j  in range(1,end[1]-start[1]+1):
        submap[(start[0],start[1]+j,start[2])]=submap[(start[0],start[1]+j-1,start[2])]+new_map(mapping,(0,start[1]+j,0))
    for k  in range(1,end[2]-start[2]+1):
        submap[(start[0],start[1],start[2]+k)]=submap[(start[0],start[1],start[2]+k-1)]+new_map(mapping,(0,0,start[2]+k))
    '''

    for i in range(0,end[0]-start[0]+1):
        for j in range(0,end[1]-start[1]+1):
            for k in range(0,end[2]-start[2]+1):
                if i+j+k>=0:
                    maxim,come_from = max_of_cube(mapping,submap,[start[0]+i,start[1]+j,start[2]+k])
                    submap[(start[0]+i,start[1]+j,start[2]+k)]=maxim
                    trace_back[(start[0]+i,start[1]+j,start[2]+k)]=come_from
    subpath=[]
    subpath.extend([end])
    
    while subpath[0]!=start:
        subpath.insert(0,trace_back[(subpath[0][0],subpath[0][1],subpath[0][2])])        
    return submap,subpath




def dynamic_merge(mapping):
    '''
    dat=[]
   
    for x in struct:
        dat.append([x[0][0],x[0][1],x[0][2],x[1]])
        
    
    dat_2=[]
    dat_1=[]
    dat_3=[]
    for seq in dat:
        if seq[3]==3:
            dat_3.append(seq[:-1])
        if seq[3]==1:
            dat_1.append(seq[:-1])
        if seq[3]==2:
            dat_2.append(seq[:-1])
    
    dat_3_no_0=[]
    for seq in dat_3:
        if seq[0]*seq[1]*seq[2]!=0:
            dat_3_no_0.append(seq)
    '''

    path=[(-1,-1,-1)]
    mapping[(-1,-1,-1)]=3
    #path.append([0,0,0])
    #dat_3_no_0.insert(0,[0,0,0])
    #mapping = {}
    #mapping[(0,0,0)]=3
    #for item in dat:
        #mapping[(item[0],item[1],item[2])] = np.power(item[3],2/3.0)
    matching=[(-1,-1,-1)]
    for key in sorted(mapping.keys()):
        if mapping[key]==3 and key[0]!=-1 and key[1]!=-1 and key[2]!=-1:
            matching.append(key)

    for i in range(len(matching)-1):
        if (matching[i+1][0]-matching[i][0])*(matching[i+1][1]-matching[i][1])*(matching[i+1][2]-matching[i][2])==1:
            path.append(matching[i])
	    path.append(matching[i+1])
        else:            
            submap,subpath=dynamic_program(matching[i],matching[i+1],mapping)
            path.extend(subpath)
    new_path=[]
    [new_path.append(i) for i in path if not i in new_path]
    match=[]
    for i in range(1,len(new_path)):
        if new_path[i][0]-new_path[i-1][0]==1 and new_path[i][1]-new_path[i-1][1]==1 and new_path[i][2]-new_path[i-1][2]==1:
            match.append(new_path[i])
        elif new_path[i][0]-new_path[i-1][0]==1 and new_path[i][1]-new_path[i-1][1]==1 and new_path[i][2]-new_path[i-1][2]==0:
            match.append([new_path[i][0],new_path[i][1],'-'])
        elif new_path[i][0]-new_path[i-1][0]==1 and new_path[i][1]-new_path[i-1][1]==0 and new_path[i][2]-new_path[i-1][2]==1:
            match.append([new_path[i][0],'-',new_path[i][2]])
        elif new_path[i][0]-new_path[i-1][0]==0 and new_path[i][1]-new_path[i-1][1]==1 and new_path[i][2]-new_path[i-1][2]==1:
            match.append(['-',new_path[i][1],new_path[i][2]])
        elif new_path[i][0]-new_path[i-1][0]==1 and new_path[i][1]-new_path[i-1][1]==0 and new_path[i][2]-new_path[i-1][2]==0:
            match.append([new_path[i][0],'-','-'])
        elif new_path[i][0]-new_path[i-1][0]==0 and new_path[i][1]-new_path[i-1][1]==1 and new_path[i][2]-new_path[i-1][2]==0:
            match.append(['-',new_path[i][1],'-'])
        elif new_path[i][0]-new_path[i-1][0]==0 and new_path[i][1]-new_path[i-1][1]==0 and new_path[i][2]-new_path[i-1][2]==1:
            match.append(['-','-',new_path[i][2]])            
    seq1=[c[0] for c in match]
    seq2=[c[1] for c in match]
    seq3=[c[2] for c in match]
    return seq1,seq2,seq3

sq1,sq2,sq3=dynamic_merge(struct1)


name1,align1=read_fasta(aln1)
name2,align2=read_fasta(aln2)
name3,align3=read_fasta(aln3)

new_aln=[]
for seq in align1:
        new_seq = gap_propagate(sq1,seq)
        new_aln.append(new_seq)
for seq in align2:
        new_seq = gap_propagate(sq2,seq)
        new_aln.append(new_seq)
for seq in align3:
        new_seq = gap_propagate(sq3,seq)
        new_aln.append(new_seq)

name=[]
for seq in name1:
        name.append(seq)
for seq in name2:
        name.append(seq)
for seq in name3:
        name.append(seq)

write_fasta(newAln_file,name,new_aln)
print('finish dynamic merging')


