from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import CentriVision.bez as bez
import numpy as np
import pandas as pd
import os
import re

class Get_repeat():
	def __init__(self, options):
		self.workpath = os.getcwd()+'/'
		self.idtag = 'ID'
		self.classtag = 'Classification'
		self.gene_r8_Separator1 = ';'
		self.gene_r8_Separator2 = '='
		for k, v in options:
			setattr(self, str(k), v)
			print(k, ' = ', v)

	def read_r8(self,s,maker):
		dic = {}
		lt = s.split(self.gene_r8_Separator1)
		for i in lt:
			if self.gene_r8_Separator2 in i:
				lt0 = i.split(self.gene_r8_Separator2)
				dic[lt0[0]] = lt0[1]
		return dic[maker]

	def makename(self,chro,index):
		return chro+'r'+str(index+1)

	def run(self):
		genome = SeqIO.to_dict(SeqIO.parse(self.genome, "fasta"))# 提取之后直接返回字典
		repeat  = pd.read_csv(self.repeat_gff3,header = None, sep='\t', comment='#')
		repeat[1] = repeat.apply(lambda x: self.read_r8(x[8],self.idtag), axis=1)
		repeat[8] = repeat.apply(lambda x: self.read_r8(x[8],self.classtag), axis=1)
		repeat.reset_index(inplace=True)
		if not os.path.exists(self.workpath+self.out_path):
			os.makedirs(self.workpath+self.out_path)
		repeat.reset_index(drop=True,inplace=True)
		repeat[9] = repeat.apply(lambda x: self.makename(x[0],x['index']), axis=1)
		del repeat['index']
		repeat.to_csv('new-'+self.repeat_gff3, index=False, header=False,sep='\t',mode='a+')
		# print(repeat)
		dic = repeat.groupby(8).groups
		for key in dic.keys():
			print(" ******** ",key," ******** ")
			local = repeat.loc[dic[key]].sort_values(by=[0,3],ascending= [True,True])
			local.reset_index(drop=True,inplace=True)
			# print(local)
			name = re.sub(r'[/*.,?]', '_', key)
			f = open(self.out_path+'/'+name+'.fasta','w')
			for index,row in local.iterrows():
				seq0 = str(genome[row[0]].seq)[row[3]-1:row[4]].upper()
				f.write('>'+row[9]+'\t'+row[1]+'\n'+seq0+'\n')
			f.close()
