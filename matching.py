#!/bin/bash

GAP=-1

def printMap(M_ABC):
	keys = sorted(M_ABC.keys(),key=firstkey)
	for key in keys:
		print(str(key) + " " + str(M_ABC[key]))

def firstkey(item):
	return item[0]

def secondkey(item):
	return item[1]

def match_AB(A,B,AB):
	matching = []

	i=0
	j=0

	for k in range(len(AB[0])):
		all_gap = False
		for m in range(len(A)):
			if AB[m][k] != '-':
				while A[m][i] != AB[m][k]:
					#print(A[m][i],AB[m][k])
					i += 1
				break
		if m == len(A)-1:
			i1 = GAP
		else:
			i1 = i
			i += 1

		for m in range(len(A),len(AB)):
			if AB[m][k] != '-':
				while B[m-len(A)][j] != AB[m][k]:
					j += 1
				break

		if m == len(AB)-1:
			# B of AB all gap
			j1 = GAP
		else:
			j1 = j
			j += 1
		
#		if not (i1,j1) in matching: 
#			matching[(i1,j1)] = 0
		matching.append((i1,j1))
	return matching 

def match_idx(i,j,k,idx_type):
	if idx_type == "jki":
		tempj = j
		tempk = k
		tempi = i
		j = tempi
		i = tempk
		k = tempj
	if idx_type == "kij":
		tempk = k
		tempi = i
		tempj = j
		k = tempi
		i = tempj
		j = tempk
	return i,j,k

def match_AB_BC(M_AB,M_BC,M_ABC,idx_type="ijk"):
	M_AB.sort(key=secondkey)
	M_BC.sort(key=firstkey)


	j1 = M_BC[0][0]
	m = 0

	while j1 == -1:
		i = -1
		k = M_BC[m][1]
		i,j,k = match_idx(i,j1,k,idx_type)
		M_ABC[(i,j,k)] = 1 if (i,j,k) not in M_ABC else M_ABC[(i,j,k)]+1		
		m = m+1
		j1 = M_BC[m][0]

	for (i,j) in M_AB:
	#	print("j = " + str(j))
		if j == -1:
			k = -1
		else:
			while j1 != j:
				m = m+1
		#		print("j1 = " + str(j1))
				j1 = M_BC[m][0]
			k = M_BC[m][1]
		i,j,k = match_idx(i,j,k,idx_type)
		M_ABC[(i,j,k)] = 1 if (i,j,k) not in M_ABC else M_ABC[(i,j,k)]+1		
		

from sys import argv
from sequence_lib import read_fasta

fileA = argv[1]
fileB = argv[2]
fileC = argv[3]
fileAB = argv[4]
fileBC = argv[5]
fileCA = argv[6]

nameA, alnA = read_fasta(fileA)
nameB, alnB = read_fasta(fileB)
nameC, alnC = read_fasta(fileC)
nameAB, alnAB = read_fasta(fileAB)
nameBC, alnBC = read_fasta(fileBC)
nameCA, alnCA = read_fasta(fileCA)

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
match_AB_BC(M_BC,M_CA,M_ABC,idx_type="jki")
#print M_ABC[(1948,-1,-1)]
#printMap(M_ABC)
match_AB_BC(M_CA,M_AB,M_ABC,idx_type="kij")
#print M_ABC[(1948,-1,-1)]
printMap(M_ABC)

