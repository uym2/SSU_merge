#! /usr/bin/env python
from sequence_lib import gap_propagate, read_fasta, write_fasta
from sys import argv
import numpy as np
import subprocess

aln1 = argv[1]
aln2 = argv[2]
aln3 = argv[3]
newAln_file = argv[4]
'''
aln4 = argv[4]
aln5 = argv[5]
aln6 = argv[6]
'''

def get_seq(aln):
    seq_name=[]
    seq_align=[]
    #i=0    
    with open(aln,'r') as f:
        for x in f:
            if x[0]=='>':
                #seq_name[i]=x.rstrip()
                seq_name.append(x.rstrip())
            else:
                #seq_align[i]=x.rstrip()
                seq_align.append(x.rstrip())                
          
        
    return seq_name,seq_align
      #align aln1 and aln3 through aln2
'''   
    aln1_2 = aln1 + "_" + aln2
    aln2_3 = aln2 + "_" + aln3
    aln_merge_through2 = aln1_2 + "_" + aln2_3
    subprocess.check_call(["java","-jar","$Opal","--in",aln1,"--in2",aln2,"--out",aln1_2])
    subprocess.check_call(["java","-jar","$Opal","--in",aln2,"--in2",aln3,"--out",aln2_3])
    subprocess.check_call(["merge_in.py",aln1_2,aln2_3,aln_merge_through2])
'''  

def merge_matrix(aln1,aln2,aln3,aln1_2,aln2_3):
    '''   
    aln1_2 = 'temp1'
    aln2_3 = 'temp2'
    aln_merge_through2 = 'temp3'	

    subprocess.check_call(["java","-jar","/calab_data/mirarab/home/ziyang96/Tools/opal_2.1.3/Opal.jar","--in",aln1,"--in2",aln2,"--out",aln1_2])
    subprocess.check_call(["java","-jar","/calab_data/mirarab/home/ziyang96/Tools/opal_2.1.3/Opal.jar","--in",aln2,"--in2",aln3,"--out",aln2_3])
    '''
    aln_merge_through2 = 'temp3'
    subprocess.check_call(["merge_in.py",aln1_2,aln2_3,aln_merge_through2])
    
    aln1_name,aln1_align=get_seq(aln1)
    aln2_name,aln2_align=get_seq(aln2)   
    aln3_name,aln3_align=get_seq(aln3)
    aln1_2_name,aln1_2_align=get_seq(aln1_2)
    aln2_3_name,aln2_3_align=get_seq(aln2_3)
    aln_merge_name,aln_merge_align=get_seq(aln_merge_through2)
    I=len(aln1_name)
    J=len(aln2_name)
    K=len(aln3_name)    
    I_J=len(aln1_2_name)
    J_K=len(aln2_3_name)
    I_J_K=len(aln_merge_name)
    mapping_to1=list(range(I))
    mapping_to2=list(range(J))
    mapping_to3=list(range(K))
    for row in range(I_J_K):
        for row1 in range(I):
            if aln1_name[row1]==aln_merge_name[row]:
                mapping_to1[row1]=row
        for row2 in range(J):
            if aln2_name[row2]==aln_merge_name[row]:
                mapping_to2[row2]=row 
        for row3 in range(K):
            if aln3_name[row3]==aln_merge_name[row]:
                mapping_to3[row3]=row
                
    mapping_col_1=[-1]*len(aln_merge_align[0])
    mapping_col_2=[-1]*len(aln_merge_align[0])
    mapping_col_3=[-1]*len(aln_merge_align[0])
    i=0
    j=0
    k=0
    ii=0
    jj=0
    kk=0
    for col in range(ii,len(aln_merge_align[0])):
        for col1 in range(i,len(aln1_align[0])):
            flag=True
            for row1 in range(I):
                if aln1_align[row1][col1].upper()!=aln_merge_align[mapping_to1[row1]][col]:
                    flag=False                    
                    break
            if flag==True:
                mapping_col_1[col]=col1
                i=col1+1
		ii=col+1
                break
    for col in range(jj,len(aln_merge_align[0])):
        for col2 in range(j,len(aln2_align[0])):
            flag=True
            for row2 in range(J):
                if aln2_align[row2][col2].upper()!=aln_merge_align[mapping_to2[row2]][col]:
                    flag=False                    
                    break
            if flag==True:
                mapping_col_2[col]=col2
                j=col2+1
		jj=col+1
                break
    for col in range(kk,len(aln_merge_align[0])):
        for col3 in range(k,len(aln3_align[0])):
            flag=True
            for row3 in range(K):
                if aln3_align[row3][col3].upper()!=aln_merge_align[mapping_to3[row3]][col]:
                    flag=False                    
                    break
            if flag==True:
                mapping_col_3[col]=col3
                k=col3+1
		kk=col+1
                break   
            
    dim1=[]
    dim2=[]
    dim3=[]
    for col in range(len(aln_merge_align[0])):
        #if  mapping_col_1[col]*mapping_col_2[col]*mapping_col_3[col]>0:
        #if  mapping_col_1[col]*mapping_col_2[col]+mapping_col_1[col]*mapping_col_3[col]+mapping_col_2[col]*mapping_col_3[col]>0:
	if  mapping_col_1[col]+mapping_col_2[col]+mapping_col_3[col]>-3:
	    dim1.append(mapping_col_1[col])
            dim2.append(mapping_col_2[col])
            dim3.append(mapping_col_3[col])
    #subprocess.check_call(['rm',aln1_2,aln2_3,aln_merge_through2])
    return dim1,dim2,dim3
                               
                                
                                
def count_rep(lis):
    rep_time=[]
    lis_no_repeat=[]
    #print(lis)
    lis_no_repeat.append(lis[0])    
    for i in range(len(lis)-1):
        if lis[i]!=lis[i+1]:
            lis_no_repeat.append(lis[i+1])
    for i in lis_no_repeat:
        rep_time.append(lis.count(i))
    return lis_no_repeat,rep_time

                    
def tri_matrix(aln1,aln2,aln3):
    subprocess.check_call(["java","-jar","/calab_data/mirarab/home/ziyang96/Tools/opal_2.1.3/Opal.jar","--in",aln1,"--in2",aln2,"--out","1_2"])
    subprocess.check_call(["java","-jar","/calab_data/mirarab/home/ziyang96/Tools/opal_2.1.3/Opal.jar","--in",aln2,"--in2",aln3,"--out","2_3"])
    subprocess.check_call(["java","-jar","/calab_data/mirarab/home/ziyang96/Tools/opal_2.1.3/Opal.jar","--in",aln3,"--in2",aln1,"--out","3_1"])
    # create a merge matrix between those opal referring merge_in in the alignment.py
    dim1_2,dim2_2,dim3_2=merge_matrix(aln1,aln2,aln3,"1_2","2_3")
    dim2_3,dim3_3,dim1_3=merge_matrix(aln2,aln3,aln1,"2_3","3_1")
    dim3_1,dim1_1,dim2_1=merge_matrix(aln3,aln1,aln2,"3_1","1_2")


    vect_1=[[0,0,0] for x in range(len(dim1_1))]
    vect_2=[[0,0,0] for x in range(len(dim1_2))]
    vect_3=[[0,0,0] for x in range(len(dim1_3))]
    for i in range(len(dim1_1)):
        vect_1[i]=[dim1_1[i],dim2_1[i],dim3_1[i]]
    for i in range(len(dim1_2)):
        vect_2[i]=[dim1_2[i],dim2_2[i],dim3_2[i]]
    for i in range(len(dim1_3)):
        vect_3[i]=[dim1_3[i],dim2_3[i],dim3_3[i]]
    vect_1.extend(vect_2)
    vect_1.extend(vect_3)
    vect=sorted(vect_1)
    #print(dim1_1)
    vec,rep_time=count_rep(vect)
    struct=[]     
    #matrix=sorted(vect_1+vect_2+vect_3)
    for i in range(len(vec)):
	struct.append([vec[i],rep_time[i]])
    return struct
    '''
	 dim1_matrix=sorted(dim1_1+dim1_2+dim1_3)
    dim2_matrix=sorted(dim2_1+dim2_2+dim2_3)
    dim3_matrix=sorted(dim3_1+dim3_2+dim3_3)
    
    return dim1_matrix,dim2_matrix,dim3_matrix
    '''
struct1=tri_matrix(aln1,aln2,aln3)
print(struct1)
#print('finish constructing tri_matrix')
#for i in range(len(v1)):
#    print(v1[i],rep1[i])


def conservative_merge(struct):
    dat=[]
    for x in struct:
        #s1,s2,s3,s4=x.split(",")
        n1 = int(x[0][0])
        n2 = int(x[0][1])
        n3 = int(x[0][2])
        n4 = int(x[1])
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
            
    match_site=[] 
    i=0      
    for site in range(1,max(c[0] for c in dat)+1):
        if i <len(dat_3_no_0): 	
   	    if dat_3_no_0[i][0] == site:
            	
	    	match_site.append(dat_3_no_0[i])
            	i=i+1	
	    
            else :
                match_site.append([site,"-","-"])
	else :
	    match_site.append([site,"-","-"])    
        
    site=1
    i=0
    while site<=max(c[1] for c in dat)+1:
	if i < len(match_site):
	        if match_site[i][1]=='-':
	            i=i+1
	        elif match_site[i][1]!=site:
	            match_site.insert(i,['-',site,'-'])
	            site=site+1
	            i=i+1
	        elif match_site[i][1]==site:
	            site=site+1
	            i=i+1
	else:
		match_site.insert(i,['-',site,'-'])
                site=site+1
                i=i+1
    #print(site)            


    site=1
    i=0
    while site<=max(c[2] for c in dat)+1:
	if i < len(match_site):
	        if match_site[i][2]=='-':
	            i=i+1
	        elif match_site[i][2]!=site:
	            match_site.insert(i,['-','-',site])
	            site=site+1
	            i=i+1
	        elif match_site[i][2]==site:
	            site=site+1
	            i=i+1
	else:
		match_site.insert(i,['-','-',site])
                site=site+1
                i=i+1
    #print(site)
    seq1=[c[0] for c in match_site]
    seq2=[c[1] for c in match_site]
    seq3=[c[2] for c in match_site]
    return seq1,seq2,seq3
'''
sq1,sq2,sq3=conservative_merge(struct1)


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
'''
print('finish conservative merging')
