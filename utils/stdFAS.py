#! /usr/bin/env python

import sys

ssu_afa = sys.argv[1]
if len(sys.argv) > 2:
	fout = open(sys.argv[2],'w')
else:
	fout = sys.stdout

isFirst = True

with open(ssu_afa,'r') as fin:
	for line in fin:
		if line[0] == '>':
			if not isFirst:			
				fout.write('\n')
			else:
				isFirst = False
			fout.write(line)
		else:
			fout.write(line[:-1])

fout.close()
