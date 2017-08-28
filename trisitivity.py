#! /usr/bin/env python

from dendropy import Tree
from sys import argv
import subprocess


print("Hello. I am doing my job!")

treefile = argv[1]
directory = argv[2]

spTree = Tree.get(path=treefile,schema="newick")

def trisitivity(spTree,directory):
	num_of_nodes=0
	for node in spTree.postorder_node_iter():
		node.process=1
		num_of_nodes= num_of_nodes+1
		node.original=node.taxon.label if node is_leaf() else node.label
	for node in spTree.postorder_node_iter():
		if num_of_nodes==2:
			node_lis=[]
                        for nds in spTree.postorder_node_iter():
				if nds.process==1:
					node_lis.append(nds)
			label1=directory+'/'+(node_lis[0].taxon.label if node_lis[0].is_leaf() else node_lis[0].label)+ '.fasta'
			label2=directory+'/'+(node_lis[1].taxon.label if node_lis[1].is_leaf() else node_lis[1].label)+ '.fasta'
			print(label1,label2)
			aln_merge=directory+ '/' +(node_lis[0].taxon.label if node_lis[0].is_leaf() else node_lis[0].label) + (node_lis[1].taxon.label if node_lis[1].is_leaf() else node_lis[1].label)+ '.fasta'
			try:
				subprocess.check_call(['java','-jar','/calab_data/mirarab/home/ziyang96/Tools/opal_2.1.3/Opal.jar','--in',label1,'--in2',label2,'--out',aln_merge])
			except:
				subprocess.check_call(['python','sample_merge.py',label1,label2,aln_merge])
			if node_lis[0].level()==0:
				node_lis[0].label=(node_lis[0].taxon.label if node_lis[0].is_leaf() else node_lis[0].label) + (node_lis[1].taxon.label if node_lis[1].is_leaf() else node_lis[1].label)
			else: 
				node_lis[1].label=(node_lis[0].taxon.label if node_lis[0].is_leaf() else node_lis[0].label) + (node_lis[1].taxon.label if node_lis[1].is_leaf() else node_lis[1].label) 
			node_lis[0].process=0
			node_lis[1].process=0
			num_of_nodes=num_of_nodes-1
		if num_of_nodes>2:
			if not node.is_leaf():	
				sum=0
				for child in node.child_node_iter():
					sum=sum+child.process
				while node.num_child_nodes() >=2 and sum>=2:
					child_list=node.child_nodes()
					children=[]
					num=0
					for x in child_list:
						if x.process==1 and num<2:
							children.append(x)
							num=num+1
						elif num>=2:
							break
					label1=directory + '/'+ (children[0].taxon.label if children[0].is_leaf() else children[0].label)+'.fasta'
					label2=directory + '/'+ (children[1].taxon.label if children[1].is_leaf() else children[1].label)+'.fasta'
					print('type1')
					print(label1)
					print(label2)
					parent=directory + '/' + node.label + '.fasta'
					print(parent)
					aln_merge=directory+'/'+ (children[0].taxon.label if children[0].is_leaf() else children[0].label) +(children[1].taxon.label if children[1].is_leaf() else children[1].label) + node.label + '.fasta'
					node.label=(children[0].taxon.label if children[0].is_leaf() else children[0].label) + (children[1].taxon.label if children[1].is_leaf() else children[1].label) + node.label 
					subprocess.check_call(['python','matching.py',label1,label2,parent, aln_merge])
					children[0].process=0
					children[1].process=0
					num_of_nodes=num_of_nodes-2
					sum=0
                	       		for child in node.child_node_iter():
                        	       		 sum=sum+child.process

				if node.num_child_nodes()==1 or sum ==1:
					child_list=node.child_nodes()
	                                for x in node.child_nodes():
						if x.process ==1:
							child1=x
							break
					label=directory + '/'+ (child1.taxon.label if child1.is_leaf() else child1.label)+'.fasta'
					itself=directory + '/' + node.label + '.fasta'
					parent=directory + '/' + node.parent_node.label + '.fasta'
					print('type2')
					print(label,itself,parent)
					aln_merge=directory+'/'+ (child1.taxon.label if child1.is_leaf() else child1.label)+node.label +node.parent_node.label+'.fasta'
					node.parent_node.label=(child1.taxon.label if child1.is_leaf() else child1.label)+node.label +node.parent_node.label
					subprocess.check_call(['python','matching.py',label,itself,parent, aln_merge])
					child1.process=0
					node.process=0
					num_of_nodes=num_of_nodes-2
	for node in spTree.postorder_node_iter():
		if node.level()==0:
			path=directory+'/' + node.label + '.fasta'
			return path

						
trisitivity(spTree,directory)
				
				
