from sequence_lib import read_fasta
from math import log
import subprocess

class smplMerger:
	def __init__(self,alnFile1,alnFile2):
		name1, aln1 = read_fasta(alnFile1)
		name2, aln2 = read_fasta(alnFile2)
		self.name1 = name1
		self.name2 = name2
		self.aln1 = aln1
		self.aln2 = aln2
		self.alnFile1 = alnFile1
		self.alnFile2 = alnFile2

	def sub_score(self,sub_aln1,sub_aln2,sub_merged,scoring):
		i=0
		j=0


		for k in range(len(sub_merged[0])):
			#print(k)
			all_gap = False
			for m in range(len(sub_aln1)):
				if sub_merged[m][k] != '-':
					while sub_aln1[m][i] != sub_merged[m][k]:
						#print(sub_aln1[m][i],sub_merged[m][k])
						i += 1
					break
			if m == len(sub_aln1)-1:
				i1 = -1
			else:
				i1 = i
				i += 1
	
			for m in range(len(sub_aln1),len(sub_merged)):
				if sub_merged[m][k] != '-':
					while sub_aln2[m-len(sub_aln1)][j] != sub_merged[m][k]:
						j += 1
					break

			if m == len(sub_merged)-1:
				# sub_aln2 of sub_merged all gap
				j1 = -1
			else:
				j1 = j
				j += 1
			
			if not (i1,j1) in scoring: 
				scoring[(i1,j1)] = 0
			scoring[(i1,j1)] += 1
			#print(i1,j1)
			#print(scoring[(i1,j1)])
					

	def sub_merge(self,n1,n2):
		# randomly sample sequences from aln1 and aln2
		print(n1)
		print(n2)
		subprocess.check_call(["sampling.py",self.alnFile1,"temp1.fas",str(n1)])	
		subprocess.check_call(["sampling.py",self.alnFile2,"temp2.fas",str(n2)])	
		# and call opal (or perhaps another merger) to merge them
		subprocess.check_call(["java","-jar","/calab_data/mirarab/home/ziyang96/Tools/SSU_merge/opal_2.1.3/Opal.jar","--in","temp1.fas","--in2","temp2.fas","--out","temp.fas"])
		subprocess.check_call(["/calab_data/mirarab/home/ziyang96/Tools/SSU_merge/utils/stdFAS.py","temp.fas","merged.fas"])
		
		name1,sub_aln1 = read_fasta("temp1.fas")
		name2,sub_aln2 = read_fasta("temp2.fas")
		name,sub_merged = read_fasta("merged.fas")

		subprocess.check_call(["rm","temp1.fas"])
		subprocess.check_call(["rm","temp2.fas"])
		subprocess.check_call(["rm","temp.fas"])
		subprocess.check_call(["rm","merged.fas"])
		
		return sub_aln1, sub_aln2, sub_merged

	def get_score(self,nsmpl,n1,n2):
		scoring = {}
		for i in range(nsmpl):
			print(i)
			sub_aln1,sub_aln2,sub_merged = self.sub_merge(n1,n2)
			self.sub_score(sub_aln1,sub_aln2,sub_merged,scoring)
		
		return scoring

	def merge(self,scoring):
		n = len(self.aln1[0])
		m = len(self.aln2[0])

		aln_score = [[0 for i in range(m+1)] for j in range(n+1)]
		backtrack = [['-' for i in range(m+1)] for j in range(n+1)]

		for i in range(1,m+1):
			backtrack[0][i] = 'L'
			#aln_score[0][i] = aln_score[0][i-1] + scoring[(-1,i-1)] if (-1,i-1) in scoring else aln_score[0][i-1]
			aln_score[0][i] = aln_score[0][i-1]
		for j in range(1,n+1):
			backtrack[j][0] = 'U'
			#aln_score[j][0] = aln_score[j-1][0] + scoring[(j-1,-1)] if (j-1,-1) in scoring else aln_score[j-1][0]
			aln_score[j][0] = aln_score[j-1][0] 

		for j in range(1,n+1):
			for i in range(1,m+1):
				ms = aln_score[j-1][i-1] + scoring[(j-1,i-1)] if (j-1,i-1) in scoring else aln_score[j-1][i-1] 
#				g1 = aln_score[j][i-1] + scoring[(-1,i-1)] if (-1,i-1) in scoring else aln_score[j][i-1]
#				g2 = aln_score[j-1][i] + scoring[(j-1,-1)] if (j-1,-1) in scoring else aln_score[j-1][i]
				g1 = aln_score[j][i-1]
				g2 = aln_score[j-1][i]
		
				if ms >= g1 and ms >= g2:
					aln_score[j][i] = ms
					backtrack[j][i] = 'D'
				elif g1 >= ms and g1 >= g2:
					aln_score[j][i] = g1
					backtrack[j][i] = 'L'
				else:
					aln_score[j][i] = g2
					backtrack[j][i] = 'U'

		i = m
		j = n
		M1 = ""
		M2 = ""
		while (i > 0 or j > 0):
			#print(aln_score[j][i])
			if backtrack[j][i] == 'D':
				#print('D')
				#M1 = str(j-1) + M1
				#M2 = str(i-1) + M2
				M1 = "." + M1
				M2 = "." + M2
				i -= 1
				j -= 1	
			elif backtrack[j][i] == 'L':
				#print('L')
				#M2 = str(i-1) + M2
				M2 = "." + M2
				M1 = "-" + M1
				i -= 1
			else:
				#print('U')
				M2 = "-" + M2
				#M1 = str(j-1) + M1
				M1 = "." + M1
				j -= 1
		return aln_score[n][m], M1, M2		

	def smpl_merge(self,nsmpl=100,maxN=2000,r1=0.1,r2=0.1):
		#n1 = min(maxN,int(len(self.aln1)*r1))
		#n2 = min(maxN,int(len(self.aln2)*r2))
		n1= min(maxN, int(len(self.aln1)*r1), 500)
		n2= min(maxN, int(len(self.aln2)*r2), 500)
		scoring = self.get_score(nsmpl,n1,n2)
		return self.merge(scoring)
