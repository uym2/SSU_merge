#!/usr/bin/env python

from sequence_lib import gap_propagate, read_fasta, write_fasta
from sys import argv

cons_file = argv[1]
subAln_file = argv[2]
newAln_file = argv[3]

names,aln = read_fasta(subAln_file)

cons = []
with open(cons_file,'r') as f:
	for line in f:
		cons.append(line.rstrip())

new_aln = []
for seq in aln:
	new_seq = gap_propagate(cons,seq)
	new_aln.append(new_seq)

write_fasta(newAln_file,names,new_aln)

