import re,os,math
import sys
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
from tqdm import tqdm
import CentriVision.bez as bez

class EdistDot():
	def __init__(self, options):
		self.similarity = 0.8
		self.evalue = 1e-5
		self.markersize = 5
		self.block = 4
		self.block_gep = 0
		self.figsize = 'default'
		self.position = 'order'
		self.blast_reverse = 'False'
		self.genome_name_size = '30'
		self.chr_name_size = '20'
		self.tandem = 'True'
		self.levels = '1:1:0'
		self.score_levels = 'False'
		self.q_s = "1:1"
		for k, v in options:
			setattr(self, str(k), v)
			print(k, ' = ', v)
		self.q_s = [int(i) for i in self.q_s.split(':')]
		self.levels = [int(i)*self.q_s[1] for i in self.levels.split(':')]
		self.similarity = float(self.similarity)

	def get_chr(self,gene,gff):
		if gene in gff.keys():
			return gff[gene]['chr']
		else:
			return 'None'

	def get_index(self,gene,gff):
		if gene in gff.keys():
			return gff[gene]['order']
		else:
			return 'None'

	def cat_gene(self,gene1,gene2):
		return '_'.join([gene1,gene2])

	def get_tandem(self,x,list0):
		return min([j for j in [abs(x - i) for i in list0] if j != 0])

	def readiblast_large(self,genepairs_,blast_reverse,score_levels_,position_,score,tandem_,gff,lens,chr1_start,chr2_start,step1,step2):
		print("Count the number of lines in the BLAST file")
		start_time = time.time()  # 记录开始时间
		total_lines = sum(1 for _ in open(genepairs_,'r'))
		end_time = time.time()  # 记录结束时间
		print(f"Function execution time: {end_time - start_time:.2f} seconds")
		x1,x2,x3,y1,y2,y3 = [],[],[],[],[],[]
		# 初始化进度条

		with tqdm(total=math.ceil(total_lines/100000), desc="Reading BLAST", unit="block") as pbar:
			# 分块读取文件
			chunk_iter = pd.read_csv(genepairs_,header = 0, sep='\t', comment='#', chunksize=100000)
			processed_lines = 0
			tail = pd.DataFrame()
			for chunk in chunk_iter:
				processed_lines += len(chunk)
				chunk.loc[:, 4] = chunk.apply(lambda x : self.cat_gene(x['Seq1'],x['Seq2']), axis=1)
				chunk = chunk.drop_duplicates(subset=4 , keep='first')
				del chunk[4]
				if blast_reverse == 'True':
					chunk[['Seq1', 'Seq2']] = chunk[['Seq2', 'Seq1']]
				chunk = chunk.drop(chunk[(chunk['Seq1'] == chunk['Seq2'])].index)#.sort_values(by=[0,11],ascending= [True,False])
				last_key = list(chunk['Seq1'].to_list())[-1]
				if not tail.empty:
					chunk = pd.concat([tail, chunk], ignore_index=True)
					chunk.reset_index(drop=True,inplace=True)
				dic = chunk.groupby('Seq1').groups# 按照第几列分组
				# 判断是否为最后一块
				if processed_lines == total_lines:
					# print("This is the last chunk.")
					tail = pd.DataFrame()
				else:
					tail = chunk.loc[dic[last_key]]
					tail.reset_index(drop=True,inplace=True)
					del dic[last_key]

				for key in dic.keys():
					# print(" ******** ",key," ******** ")
					local = chunk.loc[dic[key]]
					local.reset_index(drop=True,inplace=True)

					if key not in gff.keys():
						continue
					for index,row in local.iterrows():
						if row['Similarity'] >= score:
							if row['Seq2'] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[row['Seq2']]['chr'] not in lens.keys():
								continue
							if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
								x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
							else:
								x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
							if 'rev' in lens[gff[row['Seq2']]['chr']].keys() and lens[gff[row['Seq2']]['chr']]['rev'] == '-':
								y = -(chr2_start[gff[row['Seq2']]['chr']] + lens[gff[row['Seq2']]['chr']][position_] - gff[row['Seq2']][position_] + 1)*step2
							else:
								y = -(chr2_start[gff[row['Seq2']]['chr']] + gff[row['Seq2']][position_])*step2
							if row['Direction'] == 'forward':
								x1.append(x)
								y1.append(y)
								if self.genome1_name == self.genome2_name:
									x1.append(-y)
									y1.append(-x)
							else:
								x2.append(x)
								y2.append(y)   
								if self.genome1_name == self.genome2_name:
									x2.append(-y)
									y2.append(-x)

				pbar.update()
		return x1,x2,x3,y1,y2,y3

	def filter_blast(self,blast,gff):
		blast.loc[:, 12] = blast.apply(lambda x : self.get_chr(x[0],gff), axis=1)
		blast.loc[:, 13] = blast.apply(lambda x : self.get_chr(x[1],gff), axis=1)
		# print(blast)
		blast = blast.drop(blast[(blast[12] == 'None') | (blast[13] == 'None')].index)# 删除自身匹配chr_list
		# print(blast)

		# 创建一个空的 DataFrame
		empty_df = pd.DataFrame(columns=[i for i in range(15)])

		dic = blast.groupby(0).groups# 按照第几列分组
		for key in dic.keys():
			local = blast.loc[dic[key]].sort_values(by=[13,11],ascending= [True,False])
			local.reset_index(drop=True,inplace=True)
			# print(local)
			dic0 = local.groupby(13).groups# 按照第几列分组
			for chro in dic0.keys():
				local0 = local.loc[dic0[chro]].sort_values(by=[13,11],ascending= [True,False])
				local0.reset_index(drop=True,inplace=True)
				# print(local0)
				if len(local0) == 1:
					# new.append([local0.iloc[0][0],local0.iloc[0][1]])
					empty_df = pd.concat([empty_df, local0], ignore_index=True)
				else:
					local0[14] = local0.apply(lambda x : self.get_index(x[1],gff), axis=1)
					local0[15] = local0.apply(lambda x : self.get_tandem(x[14],list(local0[14].to_list())), axis=1)
					
					local0 = local0.drop_duplicates(subset=15 , keep='first')
					
					empty_df = pd.concat([empty_df, local0], ignore_index=True)
					# del local0[12]

		print(empty_df)

		exit(0)

	def run(self):
		axis = [0, 1, 1, 0]
		left, right, top, bottom = 0.07, 0.97, 0.93, 0.03
		lens1,chr_list1 = bez.read_lens(self.lens1)# 染色体字典
		lens2,chr_list2 = bez.read_lens(self.lens2)# 染色体字典
		gff1 = bez.read_gff(self.gff1,chr_list1)# 基因字典
		gff2 = bez.read_gff(self.gff2,chr_list2)# 基因字典

		# lens = dict( lens1, **lens2)#只有key是字符串的时候有用
		lens = lens1.copy()
		lens.update(lens2)
		gff = dict( gff1, **gff2)
		chr_list = list(set(chr_list1 + chr_list2))
		lens1sum = sum(lens1[key][self.position] for key in chr_list1)
		lens2sum = sum(lens2[key][self.position] for key in chr_list2)

		step1 = 1 / float(lens1sum)
		step2 = 1 / float(lens2sum)
		chr1_start = {}
		start = 0
		for chro in chr_list1:
			chr1_start[chro] = start
			start = start + lens1[chro][self.position]
		chr2_start = {}
		start = 0
		for chro in chr_list2:
			chr2_start[chro] = start
			start = start + lens2[chro][self.position]

		# print(chr1_start)

		if re.search(r'\d', self.figsize):
			self.figsize = [float(k)*2 for k in self.figsize.split(',')]
		else:
			self.figsize = np.array(
				[1, float(lens1sum)/float(lens2sum)])*10

		plt.figure(figsize=self.figsize, dpi=600)
		path = bez.get_path()
		font_path = os.path.join(path, 'example/arial.ttf')
		from matplotlib.font_manager import FontProperties
		font_prop = FontProperties(fname=font_path)
		if self.genepairsfile_type != 'BLAST' and self.genepairsfile_type != 'EdistAln':
			colineartly = readcolineartly_dotplot(self.genepairs,int(self.block),int(self.block_gep),self.genepairsfile_type,gff,chr_list)
			x,y = [],[]
			for block in colineartly:
				for pair in block:
					# print(pair)
					if gff[pair[0]]['chr'] not in lens1.keys() and gff[pair[0]]['chr'] in lens2.keys():
						pair = pair[::-1]


					if 'rev' in lens[gff[pair[0]]['chr']].keys() and lens[gff[pair[0]]['chr']]['rev'] == '-':
						x0 = (chr1_start[gff[pair[0]]['chr']] + lens[gff[pair[0]]['chr']][self.position] - gff[pair[0]][self.position] + 1)*step1
					else:
						x0 = (chr1_start[gff[pair[0]]['chr']] + gff[pair[0]][self.position])*step1

					# x0 = (chr1_start[gff[pair[0]]['chr']] + gff[pair[0]][self.position])*step1

					if 'rev' in lens[gff[pair[1]]['chr']].keys() and lens[gff[pair[1]]['chr']]['rev'] == '-':
						# print(gff[pair[1]]['chr'],lens[gff[pair[1]]['chr']])
						# exit()
						y0 = -(chr2_start[gff[pair[1]]['chr']] + lens[gff[pair[1]]['chr']][self.position] - gff[pair[1]][self.position] + 1)*step2
					else:
						y0 = -(chr2_start[gff[pair[1]]['chr']] + gff[pair[1]][self.position])*step2

					# y0 = -(chr2_start[gff[pair[1]]['chr']] + gff[pair[1]][self.position])*step2
					x.append(x0)
					y.append(y0)
					if self.genepairsfile_type == 'MCScanX' and self.genome1_name == self.genome2_name:# MCScanX的结果只有a->b没有b->a

						if 'rev' in lens[gff[pair[1]]['chr']].keys() and lens[gff[pair[1]]['chr']]['rev'] == '-':
							x0 = (chr1_start[gff[pair[1]]['chr']] + lens[gff[pair[1]]['chr']][self.position] - gff[pair[1]][self.position] + 1)*step1
						else:
							x0 = (chr1_start[gff[pair[0]]['chr']] + gff[pair[0]][self.position])*step1


						if 'rev' in lens[gff[pair[0]]['chr']].keys() and lens[gff[pair[0]]['chr']]['rev'] == '-':
							y0 = -(chr2_start[gff[pair[0]]['chr']] + lens[gff[pair[0]]['chr']][self.position] - gff[pair[0]][self.position] + 1)*step2
						else:
							y0 = -(chr2_start[gff[pair[0]]['chr']] + gff[pair[0]][self.position])*step2

						x.append(x0)
						y.append(y0)

			plt.scatter(x,y,s = float(self.markersize), edgecolors='none', facecolors='#e2041b', alpha=0.3)
		elif self.genepairsfile_type == 'BLAST':
			x1,x2,x3,y1,y2,y3 = readblast_large(self.genepairs,self.blast_reverse,self.score_levels,self.position,self.levels,self.tandem,gff,lens,chr1_start,chr2_start,step1,step2)
			plt.scatter(x3,y3,s = float(self.markersize), edgecolors='none', facecolors='#494a41', alpha=0.5)
			plt.scatter(x2,y2,s = float(self.markersize), edgecolors='none', facecolors='#1e50a2', alpha=0.4)
			plt.scatter(x1,y1,s = float(self.markersize), edgecolors='none', facecolors='#e2041b', alpha=0.3)
		elif self.genepairsfile_type == 'EdistAln':
			x1,x2,x3,y1,y2,y3 = self.readiblast_large(self.genepairs,self.blast_reverse,self.score_levels,self.position,self.similarity,self.tandem,gff,lens,chr1_start,chr2_start,step1,step2)
			# plt.scatter(x3,y3,s = float(self.markersize), c='#494a41', alpha=0.3)
			plt.scatter(x2,y2,s = float(self.markersize), edgecolors='none', facecolors='#1e50a2', alpha=0.4)
			plt.scatter(x1,y1,s = float(self.markersize), edgecolors='none', facecolors='#e2041b', alpha=0.3)


		plt.axis ('off')
		plt.text((lens1sum/2)*step1, 0.03, self.genome1_name,ha='center',va='center', fontproperties=font_prop,fontsize=int(self.genome_name_size))
		plt.text(-0.03, (-lens2sum/2)*step2, self.genome2_name,rotation=90,ha='center',va='center', fontproperties=font_prop,fontsize=int(self.genome_name_size))
		for key in chr1_start.keys():
			x = chr1_start[key]
			plt.axvline(x=float(x)*step1, c="black", ls="-", lw=0.5)
			plt.text((x+(lens1[key][self.position]/2))*step1, 0.01, key,ha='center',va='center', fontproperties=font_prop,fontsize=int(self.chr_name_size))
		plt.axvline(x=0*step1, c="black", ls="-", lw=2.5)
		plt.axvline(x=float(lens1sum)*step1, c="black", ls="-", lw=2.5)
		for key in chr2_start.keys():
			y = chr2_start[key]
			plt.axhline(y=-float(y)*step2, c="black", ls="-", lw=0.5)
			plt.text(-0.01, (-y-(lens2[key][self.position]/2))*step2, key,ha='center',va='center',rotation=90, fontproperties=font_prop,fontsize=int(self.chr_name_size))#,rotation=30
		plt.axhline(y=-float(lens2sum)*step2, c="black", ls="-", lw=2.5)
		plt.axhline(y=-0*step2, c="black", ls="-", lw=2.5)
		plt.xlim(xmin = 0, xmax=1)
		plt.ylim(ymin = -1, ymax=0)

		plt.savefig(self.savefig, dpi=600, bbox_inches='tight')
		# plt.show()
		sys.exit(0)
