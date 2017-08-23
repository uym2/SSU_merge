#! /usr/bin/env python

from dendropy import Tree
from sys import argv
import subprocess


treefile = argv[1]
directory = argv[2]

spTree = Tree.get(path=treefile,schema="newick")

def merge_from_tree(spTree,directory):
	for node in spTree.postorder_node_iter():
		if not node.is_leaf():
			for child in node.child_node_iter():
				p_label = node.label
				c_label = child.taxon.label if child.is_leaf() else child.label
				aln_p = directory + "/" + p_label + ".fasta.ungap"
				aln_c = directory + "/" + c_label + ".fasta.ungap"
				node.label = p_label + "_" + c_label
				aln_merged = directory + "/" + node.label + ".fasta.ungap"
				subprocess.check_call(["python","/calab_data/mirarab/home/ziyang96/Tools/SSU_merge/sample_merge.py",aln_p,aln_c,aln_merged])				

merge_from_tree(spTree,directory)
