#!/bin/bash
import subprocess
from matching_lib import match_AB_BC, match_AB		
from sequence_lib import gap_propagate, read_fasta, write_fasta

from sys import argv

fileA = argv[1]
fileB = argv[2]
fileC = argv[3]
newAln_file = argv[4]

fileAB = 'temp1'
try:
	subprocess.check_call(["java","-jar","/calab_data/mirarab/home/ziyang96/Tools/opal_2.1.3/Opal.jar","--in",fileA,"--in2",fileB,"--out",fileAB])
except:
	subprocess.check_call(['python','sample_merge.py',fileA,fileB,fileAB])
fileBC = 'temp2'
try:
	subprocess.check_call(["java","-jar","/calab_data/mirarab/home/ziyang96/Tools/opal_2.1.3/Opal.jar","--in",fileB,"--in2",fileC,"--out",fileBC])
except:
        subprocess.check_call(['python','sample_merge.py',fileB,fileC,fileBC])

fileCA = 'temp3'
try:
	subprocess.check_call(["java","-jar","/calab_data/mirarab/home/ziyang96/Tools/opal_2.1.3/Opal.jar","--in",fileC,"--in2",fileA,"--out",fileCA])
except:
        subprocess.check_call(['python','sample_merge.py',fileC,fileA,fileCA])

fileAB_fix='temp4'
fileBC_fix='temp5'
fileCA_fix='temp6'
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


def strict_merging(mapping):
	matching=[]
	for key in sorted(mapping.keys()):
		if mapping[key]==3 and key[0]!=-1 and key[1]!=-1 and key[2]!=-1:
			matching.append(key)
	site=0
	max1=max([c[0] for c in mapping.keys()])
        max2=max([c[1] for c in mapping.keys()])
        max3=max([c[2] for c in mapping.keys()])
	#problematic
	while site <= max1:
		try:
			matching[site][0]
		except:
			matching.insert(site,(site,'-','-'))
                        site=site+1
		else:
			if matching[site][0]==site:
				site=site+1
			else:
				matching.insert(site,(site,'-','-'))
				site=site+1
		#print(matching[site],site)
	site=0
	i=0
	while site <= max2:
		try:
			matching[i][1]
		except:
			matching.insert(i,('-',site,'-'))
                        site=site+1
                        i=i+1
		else:	
			if matching[i][1]==site:
				site=site+1
				i=i+1
			elif matching[i][1]=='-':
				i=i+1
			else:
				matching.insert(i,('-',site,'-'))
				site=site+1
				i=i+1
	print(max1,max2,max3)
        site=0
        i=0
        while site <= max3:
                try:
                        matching[i][2]
                except:
                        matching.insert(i,('-','-',site))
                        site=site+1
                        i=i+1
                else:
                        if matching[i][2]==site:
                                site=site+1
                                i=i+1
                        elif matching[i][2]=='-':
                                i=i+1
                        else:
                                matching.insert(i,('-','-',site))
                                site=site+1
                                i=i+1

	seq1=[c[0] for c in matching]
        seq2=[c[1] for c in matching]
        seq3=[c[2] for c in matching]
	return seq1,seq2,seq3


sq1,sq2,sq3=strict_merging(M_ABC)

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
subprocess.check_call(['rm','temp1','temp2','temp3','temp4','temp5','temp6'])
print('finish conservative merging')


