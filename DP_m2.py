#! /usr/bin/env python

from sequence_lib import gap_propagate, read_fasta, write_fasta
#from tri_merge_matrix import tri_matrix
from matching_lib import match_AB_BC, match_AB
from sys import argv
import numpy as np
import subprocess
from os.path import splitext, basename


fileA = argv[1]
fileB = argv[2]
fileC = argv[3]
newAln_file = argv[4]


#fileA_name,ext = splitext(basename(fileA))
#fileB_name,ext = splitext(basename(fileB))
#fileC_name,ext = splitext(basename(fileC))

#fileAB = 'temp1'
fileAB = subprocess.check_output(['mktemp'])
try:
#        subprocess.check_call(["java","-jar","/calab_data/mirarab/home/ziyang96/Tools/opal_2.1.3/Opal.jar","--in",fileA,"--in2",fileB,"--out",fileAB])
        subprocess.check_call(["java","-jar","/Users/uym2/Packages_N_Libraries/opal_2.1.3/Opal.jar","--in",fileA,"--in2",fileB,"--out",fileAB])
except:
        subprocess.check_call(['python','sample_merge.py',fileA,fileB,fileAB])
#fileBC = 'temp2'
fileBC = subprocess.check_output(['mktemp'])
try:
#        subprocess.check_call(["java","-jar","/calab_data/mirarab/home/ziyang96/Tools/opal_2.1.3/Opal.jar","--in",fileB,"--in2",fileC,"--out",fileBC])
        subprocess.check_call(["java","-jar","/Users/uym2/Packages_N_Libraries/opal_2.1.3/Opal.jar","--in",fileB,"--in2",fileC,"--out",fileBC])
except:
        subprocess.check_call(['python','sample_merge.py',fileB,fileC,fileBC])

fileCA = subprocess.check_output(['mktemp'])
#fileCA = 'temp3'
try:
#        subprocess.check_call(["java","-jar","/calab_data/mirarab/home/ziyang96/Tools/opal_2.1.3/Opal.jar","--in",fileC,"--in2",fileA,"--out",fileCA])
        subprocess.check_call(["java","-jar","/Users/uym2/Packages_N_Libraries/opal_2.1.3/Opal.jar","--in",fileC,"--in2",fileA,"--out",fileCA])
except:
        subprocess.check_call(['python','sample_merge.py',fileC,fileA,fileCA])

#fileAB_fix='temp4'
fileAB_fix = subprocess.check_output(['mktemp'])
#fileBC_fix='temp5'
fileBC_fix = subprocess.check_output(['mktemp'])
#fileCA_fix='temp6'
fileCA_fix = subprocess.check_output(['mktemp'])

subprocess.check_call(['stdFAS.py', fileAB, fileAB_fix])
subprocess.check_call(['stdFAS.py', fileCA, fileCA_fix])
subprocess.check_call(['stdFAS.py', fileBC, fileBC_fix])

nameA, alnA = read_fasta(fileA)
nameB, alnB = read_fasta(fileB)
nameC, alnC = read_fasta(fileC)
nameAB, alnAB = read_fasta(fileAB_fix)
nameBC, alnBC = read_fasta(fileBC_fix)
nameCA, alnCA = read_fasta(fileCA_fix)

M_AB = match_AB(alnA,alnB,alnAB)
M_BC = match_AB(alnB,alnC,alnBC)
M_CA = match_AB(alnC,alnA,alnCA)

#print(M_AB)
#print(M_BC)
#print(M_CA)

M_ABC = {}

match_AB_BC(M_AB,M_BC,M_ABC)
#printMap(M_ABC)
#print M_ABC[(1948,-1,-1)]
match_AB_BC(M_BC,M_CA,M_ABC,idx_type="bca")
#print M_ABC[(1948,-1,-1)]
#printMap(M_ABC)
match_AB_BC(M_CA,M_AB,M_ABC,idx_type="cab")
#print M_ABC[(1948,-1,-1)]
#printMap(M_ABC)

print('finish constructing matrix, start dunamic programming')
'''
for key in M_ABC.keys():
        #print(key)
        if (key[0]+1)*(key[1]+1)+(key[1]+1)*(key[2]+1)+(key[2]+1)*(key[0]+1) == 0 and M_ABC[key]!=3:
            M_ABC[key]=1.5
        elif (key[0]+1)*(key[1]+1)*(key[2]+1)==0:
            M_ABC[key]=0
	elif M_ABC[key]!=3:
	    M_ABC[key]=0
	elif (key[0]+1)*(key[1]+1)*(key[2]+1)==0 and (key[0]+1)*(key[1]+1)+(key[1]+1)*(key[2]+1)+(key[2]+1)*(key[0]+1) != 0:
	    
print(M_ABC)
'''



def new_map(mapping,i,j,k,start):
    i1 = i + start[0] if i>=0 else i
    j1 = j + start[1] if j>=0 else j
    k1 = k + start[2] if k>=0 else k
    try:
        return mapping[(i1,j1,k1)]
    except:
        return 0
#    else:
#        return mapping[(input[0]+start[0],input[1]+start[1],input[2]+start[2])]
        #return mapping[input]

def new_map_list(submap,i,j,k):
    if i < 0 or j < 0 or k < 0:
        return 0
    return submap[i][j][k]

def max_of_cube(mapping,submap,start,i,j,k):
#    print(submap)
    L=[0 for x in range(7)]
#    L[0]=new_map(submap,(input[0]-1,input[1]-1,input[2]-1))+new_map(mapping,(input[0],input[1],input[2]))
#    L[1]=new_map(submap,(input[0],input[1]-1,input[2]-1))+new_map(mapping,(-1,input[1],input[2]))
#    L[2]=new_map(submap,(input[0]-1,input[1],input[2]-1))+new_map(mapping,(input[0],-1,input[2]))
#    L[3]=new_map(submap,(input[0]-1,input[1]-1,input[2]))+new_map(mapping,(input[0],input[1],-1))
#    L[4]=new_map(submap,(input[0],input[1],input[2]-1))+new_map(mapping,(-1,-1,input[2]))
#    L[5]=new_map(submap,(input[0],input[1]-1,input[2]))+new_map(mapping,(-1,input[1],-1))
#    L[6]=new_map(submap,(input[0]-1,input[1],input[2]))+new_map(mapping,(input[0],-1,-1))
    L[0]=new_map_list(submap,i-1,j-1,k-1) + new_map(mapping,i,j,k,start)
    L[1]=new_map_list(submap,i,j-1,k-1) + new_map(mapping,-1,j,k,start)
    L[2]=new_map_list(submap,i-1,j,k-1) + new_map(mapping,i,-1,k,start)
    L[3]=new_map_list(submap,i-1,j-1,k) + new_map(mapping,i,j,-1,start)
    L[4]=new_map_list(submap,i,j,k-1) + new_map(mapping,-1,-1,k,start)
    L[5]=new_map_list(submap,i,j-1,k) + new_map(mapping,-1,j,-1,start)
    L[6]=new_map_list(submap,i-1,j,k) + new_map(mapping,i,-1,-1,start)

    #print(L)

    if max(L)==L[0]:
        return max(L),(i-1,j-1,k-1)
    elif max(L)==L[1]:
        return max(L),(i,j-1,k-1)
    elif max(L)==L[2]:
        return max(L),(i-1,j,k-1)
    elif max(L)==L[3]:
        return max(L),(i-1,j-1,k)
    elif max(L)==L[4]:
        return max(L),(i,j,k-1)
    elif max(L)==L[5]:
        return max(L),(i,j-1,k)
    elif max(L)==L[6]:
        return max(L),(i-1,j,k)


def dynamic_program(start, end, mapping):
    #submap={}
    #submap[(start[0],start[1],start[2])]=0
    #trace_back={}

    #submap_dim3 = [0 for i in range(end[2]-start[2]+1)]
    #submap_dim2 = [submap_dim3 for j in range(end[1]-start[1]+1)]
    submap = [ [ [ 0 for k in range(end[2]-start[2]+1) ] for j in range(end[1]-start[1]+1) ] for i in range(end[0]-start[0]+1) ]

    #traceback_dim3 = [0 for i in range(end[2]-start[2]+1)]
    #traceback_dim2 = [traceback_dim3 for j in range(end[1]-start[1]+1)]
    #trace_back = [traceback_dim2 for k in range(end[0]-start[0]+1)]
    trace_back = [ [ [ 0 for k in range(end[2]-start[2]+1) ] for j in range(end[1]-start[1]+1) ] for i in range(end[0]-start[0]+1) ]

    #print(end[0]-start[0]+1)
    #print(end[1]-start[1]+1)
    #print(end[2]-start[2]+1)

    #print(submap)
    print(end[0]-start[0]+1)
    print(end[1]-start[1]+1)
    print(end[2]-start[2]+1)

    for i in range(0,end[0]-start[0]+1):
    	print("i = " + str(i))
        for j in range(0,end[1]-start[1]+1):
    	    #print("   j = " + str(j))
            for k in range(0,end[2]-start[2]+1):
    	    #    print("      k = " + str(k))
                if i+j+k>=0:
                    #maxim,come_from = max_of_cube(mapping,submap,[start[0]+i,start[1]+j,start[2]+k])
                    maxim,come_from = max_of_cube(mapping,submap,start,i,j,k)
             #       print(maxim)
                    #submap[(start[0]+i,start[1]+j,start[2]+k)]=maxim
                    submap[i][j][k] = maxim
                    #trace_back[(start[0]+i,start[1]+j,start[2]+k)]=come_from
                    trace_back[i][j][k] = come_from
    #print(submap)
    subpath=[]
    #subpath.extend([end])
    subpath.append((end[0]-start[0],end[1]-start[1],end[2]-start[2]))
    #while list(subpath[0])!=list(start) and list(subpath[0])!=[-1,-1,-1]:
    while list(subpath[0])!=[-1,-1,-1] and list(subpath[0])!=[0,0,0]:
       #print("subpath: ")
       #print(subpath[0])
       #print("trace back:")
       #print(trace_back[subpath[0][0]][subpath[0][1]][subpath[0][2]])
       subpath.insert(0,trace_back[subpath[0][0]][subpath[0][1]][subpath[0][2]])      
    print(subpath)
    #print([ (L[0]+start[0],L[1]+start[1],L[2]+start[2]) for L in subpath ] )
 
    return submap[end[0]-start[0]][end[1]-start[1]][end[2]-start[2]] ,[ (L[0]+start[0],L[1]+start[1],L[2]+start[2]) for L in subpath ]




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
    max1=max([c[0] for c in mapping.keys()])
    max2=max([c[1] for c in mapping.keys()])
    max3=max([c[2] for c in mapping.keys()])
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
    #matching.append((max1,max2,max3))
    #print(matching)
    #print(len(matching))
    matching.append((max1,max2,max3))



    print(max1,max2,max3)
    sum=0
    for i in range(len(matching)-1):
        print(matching[i])
        if (matching[i+1][0]-matching[i][0])*(matching[i+1][1]-matching[i][1])*(matching[i+1][2]-matching[i][2])==1:
            path.append(matching[i])
            path.append(matching[i+1])
            sum=sum+3
        else:
            #if (matching[i+1][0]-matching[i][0] > 20):
            #    continue
            subscore,subpath=dynamic_program(matching[i],matching[i+1],mapping)
            path.extend(subpath)
            sum=sum+subscore
    #print(sum)
    new_path=[]
    [new_path.append(i) for i in path if not i in new_path]
    #print(new_path)
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

sq1,sq2,sq3=dynamic_merge(M_ABC)
#print(sq1)

name1,align1=read_fasta(fileA)
name2,align2=read_fasta(fileB)
name3,align3=read_fasta(fileC)

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
#subprocess.check_call(['rm','temp1','temp2','temp3','temp4','temp5','temp6'])
subprocess.check_call(['rm',fileAB,fileBC,fileCA,fileAB_fix,fileBC_fix,fileCA_fix])
print('finish dynamic merging')


