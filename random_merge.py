#! /usr/bin/env python
# Usage: python merge_in.py <my_file> <her_file> <merge_file>
# (assume FASTA files for now, may extend later if need to)

import sys
from alignment import CompactAlignment

my_file = sys.argv[1]
her_file = sys.argv[2]
merge_file = sys.argv[3]

my_aln = CompactAlignment()
her_aln = CompactAlignment()

my_aln.read_filepath(my_file)
her_aln.read_filepath(her_file)

my_aln.merge_in(her_aln)

my_aln.write_filepath(merge_file)
