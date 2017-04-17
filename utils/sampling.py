#! /usr/bin/env python

# usage: sampling.py <file_in> <file_out> <sample_size>

import sys
from os.path import isfile
import random
import sequence_lib

try:
	import cPickle as pickle
except:
	import pickle

file_in = sys.argv[1]
file_out = sys.argv[2]
sample_size = int(sys.argv[3])

f = open(file_in,'r')
file_idx = file_in.split('.fasta')[0]+'.idx'
if not isfile(file_idx):
	sequence_lib.index_fasta(file_in,file_idx)
f_idx = open(file_idx)
seqs = pickle.load(f_idx)
f_idx.close()

seq_pointers = seqs.values()
samples = random.sample(seq_pointers,sample_size)

f_out = open(file_out,'w')

for s in samples:
	f.seek(s)
	seq = f.readline()
	f_out.write(seq)
	L = f.readline()
	while L[0] != '>':
		if L[-1] == '\n':
			f_out.write(L[:-1])
		else:
			f_out.write(L)
		L = f.readline()
		if not L:
			break
	f_out.write('\n')

f.close()
f_out.close()
