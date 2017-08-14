# I/O lib for molecular sequences

from os.path import isfile
from os import remove

try:
	import cPickle as pickle
except:
	import pickle

def get_taxon_list(filename):
	taxon_list = []
	for line in open(filename,'r'):
		if line[0] == '>':
			taxon_list = taxon_list + [line[1:].rstrip()]
	return sorted(taxon_list)


def hash_taxon_seq(filename):
	taxon_dict = {}
	f = open(filename,'r')
	for line in f:
		if line[0] == '>':
			taxon_dict[line[1:-1]] = gap_rm(f.next().rstrip())
	return taxon_dict

def gap_rm(str0,gap='-'):
	str1 = ''
	for c in str0:
		if c != gap:
			str1 =  str1 + c
	return str1

def index_fasta(file_in,file_out=None,store_index_file=True):
	# only work for fasta format
	f = open(file_in,'r')
	seq_pointers = {}
	fp = 0
	while 1:
		line = f.readline()
		if not line:
			break
		if line[0] == '>':
			seq_pointers[line[1:-1]] = fp
		fp = f.tell()

	if not store_index_file:
		return seq_pointers

	if not file_out:
		file_extension = file_in.split('.')[-1]
		file_out = file_in[:-(len(file_extension)+1)]+'.idx'
	fout = open(file_out,'w')
	pickle.dump(seq_pointers,fout)
	f.close()
	fout.close()

	return seq_pointers

def load_index(file_in,store_index_file=True,renew_index_file=False):
	file_extension = file_in.split('.')[-1]
	file_idx = file_in[:-(len(file_extension)+1)] + '.idx'

	if renew_index_file or not isfile(file_idx):
		if renew_index_file:
			remove(file_idx)
		seq_pointers = index_fasta(file_in,store_index_file=store_index_file)
	else:
		with open(file_idx) as f:
			seq_pointers = pickle.load(f)

	return seq_pointers

def sample_from_list(file_in,taxa_list,file_out,store_index_file=True,renew_index_file=False):
	seq_pointers = load_index(file_in,store_index_file=store_index_file,renew_index_file=renew_index_file)
	with open(file_in,'r') as fin:
		with open(file_out,'w') as fout:	 
			for taxon in taxa_list:
				try:
					fin.seek(seq_pointers[taxon])
					fout.write(fin.readline())
					fout.write(fin.readline())
				except:
					print ('taxon inconsistent in query and input files')


def count_gaps(seq_aln):
	N = len(seq_aln[0])
	gap_count = [0]*N
	for seq in seq_aln:
		for i in range(N):
			gap_count[i] += (seq[i] == '-')
	return gap_count

def read_fasta(fas_file):
	taxon_names = []
	seq_aln = []
	with open(fas_file,'r') as f:
		for line in f:
			if line[0] == '>':
				taxon_names.append(line[1:].rstrip())
			else:
				seq_aln.append(line.upper().rstrip())
	return taxon_names, seq_aln	

def write_fasta(output_file,taxon_names,seq_aln):
	with open(output_file,'w') as f:
		T = len(taxon_names)
		for i in range(T):
			f.write(">"+taxon_names[i]+"\n")
			f.write(seq_aln[i]+"\n")

def is_aligned(fas_file):
	taxa,seqs = read_fasta(fas_file)
	l = len(seqs[0])
	for seq in seqs:
		if len(seq) != l:
			return False
	return True

def gap_propagate(cons_seq,targ_seq):
# propagate gaps from cons_seq to targ_seq; this is used after the cons_seq was aligned with another sequence and we want to
# propage the alignment to targ_seq. This is a similar idea with transitivity used in PASTA; the ultimate goal is to merge 2 alignments,
# but in this case we have a consensus sequence for each alignments and we also have a good way (properly using seconday structure) to align them.
# NOTE: gap_propagate is NOT symmetric: only propagate from consensus to target, not in reverse; careful consider which sequence is the cons_seq and which is targ_seq!

	out_seq = ''
	i = 0
	for c in cons_seq:
		if c != '-':
			out_seq += targ_seq[i]
			i += 1
		else:
			out_seq += '-'

	return out_seq

def impose_struct(pri_seq,str_seq):
	out_pri = ''
	out_str = ''

	for i,c in enumerate(pri_seq):
		if c != '-':
			out_pri += c.upper()
			c1 = str_seq[i]
			if c1 == '(' or c1 == '<' or c1 == '{':
				out_str += '('
			elif c1 == ')' or c1 == '>' or c1 == '}':
				out_str += ')'
			else:
				out_str += '.'
	
	return out_pri, out_str
