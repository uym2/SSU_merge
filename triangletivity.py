#! /usr/bin/env python

from dendropy import Tree
from sys import argv
import subprocess
#from mergeList import mergeList

try:
	import cPickle as pickle
except:
	import pickle

print("Hello. I am doing my job!")

treefile = argv[1]
directory = argv[2]
taskfile = argv[3]

fout=open(taskfile,'w')


spTree = Tree.get(path=treefile,schema="newick")

def triangletivity(spTree,directory,fout):
	mylist=[]
	index=1
	num_of_nodes=0
	for node in spTree.postorder_node_iter():
		node.process=1
		num_of_nodes= num_of_nodes+1
		node.original= node.taxon.label if node.is_leaf() else node.label
	for node in spTree.postorder_node_iter():
		print(num_of_nodes)
		if num_of_nodes==2:
			node_lis=[]
                        for nds in spTree.postorder_node_iter():
				if nds.process==1:
					node_lis.append(nds)
			label1=directory+'/'+node_lis[0].original+ '.fasta'
			label2=directory+'/'+node_lis[1].original+ '.fasta'
			print(label1,label2)
			aln_merge=directory+ '/' +str(index)+ '.fasta'
			#subprocess.check_call(['java','-jar','/calab_data/mirarab/home/ziyang96/Tools/opal_2.1.3/Opal.jar','--in',label1,'--in2',label2,'--out',aln_merge])
			fout.write('java -jar /calab_data/mirarab/home/ziyang96/Tools/opal_2.1.3/Opal.jar --in' + ' '+label1+ ' '+ '--in2' + ' '+ label2 +' --out '+aln_merge + "\n")
			mylist.append(({node_lis[0].original,node_lis[1].original},str(index)+'.fasta'))
			node_lis[0].process=0
			node_lis[1].process=0
			num_of_nodes=num_of_nodes-1
			index=index+1
			
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
					label1=directory + '/'+ children[0].original +'.fasta'
					label2=directory + '/'+ children[1].original+'.fasta'
					print('type1')
					print(label1)
					print(label2)
					parent=directory + '/' + node.original + '.fasta'
					print(parent)
					aln_merge=directory+'/'+ str(index) + '.fasta'
					node.label=children[0].original + children[1].original + node.original 
					#subprocess.check_call(['python','dynamic_programming_symplified.py',label1,label2,parent, aln_merge])
					fout.write('python /calab_data/mirarab/home/ziyang96/Tools/SSU_merge/dynamic_programming_symplified.py ' + label1 + ' ' + label2 + ' ' + parent + ' ' + aln_merge + "\n")
					mylist.append(({children[0].original,children[1].original,node.original},str(index)+'.fasta'))
					children[0].process=0
					children[1].process=0
					num_of_nodes=num_of_nodes-2
					index=index+1
					sum=0
                	       		for child in node.child_node_iter():
                        	       		 sum=sum+child.process
					

				if node.num_child_nodes()==1 or sum ==1:
					child_list=node.child_nodes()
	                                for x in node.child_nodes():
						if x.process ==1:
							child1=x
							break
					label=directory + '/'+ child1.original+'.fasta'
					itself=directory + '/' + node.original + '.fasta'
					parent=directory + '/' + node.parent_node.original + '.fasta'
					print('type2')
					print(label,itself,parent)
					aln_merge=directory+'/'+ str(index)+'.fasta'
					node.parent_node.label=child1.original+node.original +node.parent_node.original
					#subprocess.check_call(['python','dynamic_programming_symplified.py',label,itself,parent, aln_merge])
					fout.write('python /calab_data/mirarab/home/ziyang96/Tools/SSU_merge/dynamic_programming_symplified.py' + ' '+label+ ' '+ itself+ ' '+ parent+ ' '+ aln_merge + "\n")
					mylist.append(({child1.original,node.original,node.parent_node.original},str(index)+'.fasta'))
					index=index+1
					child1.process=0
					node.process=0
					num_of_nodes=num_of_nodes-2
					
						
	return	mylist, index
	

							
mylist,index=triangletivity(spTree,directory,fout)
fout.close()
with open('myList.txt','w') as f:
	#for item in mylist:
	#	f.write(item)
	pickle.dump(mylist,f)
with open('index.txt','w') as f:
	pickle.dump(index,f)
	#f.write(str(index))


def mergeList(L,idxName,directory):
    while len(L) > 1:
        (S,S_name) = L.pop()
        for (S1,S1_name) in L:
            if (len(S & S1)):
                L.remove((S1,S1_name))
                name = str(idxName) + ".fasta"
                subprocess.check_call(["python", "merge_in.py", directory+'/'+ S_name, directory+'/'+S1_name, directory+'/'+name])
                S = (S | S1)
                L.append((S,name))
                idxName = idxName + 1
		print(idxName)
                break
    return L[0]



#merge=mergeList(mylist, index, directory)

				
				
