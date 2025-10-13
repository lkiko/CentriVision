# python
# -*- encoding: utf-8 -*-
'''
@File		:bio.py
@Time		:2021/04/27 19:54:08
@Author		:charles kiko
@Version		:1.0
@Contact		:charles_kiko@163.com
@Desc		:python bio.py gff pep cds ncbi/embl xxx
@annotation		:
'''

import sys
import ast
import os
import gc# 内存管理模块
from tqdm import trange
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import CentriVision.bez as bez

class Centrifinder():
	def __init__(self, options):
		self.workpath = os.getcwd()+'/'
		self.workpaths = './TRF/'
		self.chip_seq = 'None'
		self.score = 80
		self.colors = '#38b48b,#1e50a2,#d7003a'
		bez_conf = bez.config()
		for k, v in bez_conf:# 
			setattr(self, str(k), v)
		for k, v in options:
			setattr(self, str(k), v)
			print(str(k), ' = ', v)
		self.windows = int(self.windows)
		self.step = int(self.step)
		self.gap = int(self.gap)
		self.score = float(self.score)
		

	def plot(self,ax,x,y,w,h,c,a,p,n):
		rect = plt.Rectangle((x,y),w,h,facecolor=c,alpha=a, edgecolor='none')# x,y,宽，高
		ax.add_patch(rect)

	def covarage(self,suf_dic,length,chro):
		array_d = np.zeros(length, dtype = int)
		for start,end in suf_dic:
			start = start - 1
			array_d[start:end] = array_d[start:end]+1
		array_d1 = np.int64(array_d>0)
		return array_d1

	def get_overlap(self,gff,dic,lens,key,genome,chip_seq):
		print(" ******** ",key," ******** ")
		local = gff.loc[dic[key]].sort_values(by=[0,3],ascending= [True,True])
		local.reset_index(drop=True,inplace=True)
		print(local)

		dic0 = local.groupby(1).groups# 按照第几列分组
		# lt = ['TRF']
		lt = list(set(list(local[1].to_list())))
		TRFlocal = local.loc[dic0[lt[0]]].sort_values(by=[0,3],ascending= [True,True])
		TRFlocal.reset_index(drop=True,inplace=True)
		intervals = TRFlocal[9]
		intervals = list(map(eval, intervals))
		# print(len(intervals))
		array_d1 = self.covarage(intervals,lens[key],key)
		# print(array_d1)

		index = []
		centri = []
		starts,ends = [],[]
		x = []
		for i in range(0,len(array_d1),self.step):
			if i+self.windows == len(array_d1):
				end = -1
				s = array_d1[i:]
			elif i+self.windows > len(array_d1):
				continue
			else:
				s = array_d1[i:i+self.windows]
			# print(s)
			p = 100*sum(s)/len(s)

			if p >= self.score:
				# print(p)
				index.append(p)
				starts.append(i)
				ends.append(i+self.windows)
				end0 = i+self.windows
			else:
				index.append(p)
				if len(starts) >= 1 and (i+self.windows - end0 >= self.gap*self.windows or i+self.windows>=len(array_d1)):
					s,e = 0,0
					if end0 + self.windows <= len(array_d1):
						e = end0 + self.windows
					else:
						e = end0

					if starts[0]-self.windows < 0:
						s = starts[0]
					else:
						s = starts[0]-self.windows
					centri.append([s,e])
					starts,ends = [],[]
				elif len(starts) >= 1 and i+self.windows - ends[-1] < self.gap*self.windows:
					starts.append(i)
					ends.append(i+self.windows) 
				else:
					pass
			x.append(i+int(self.windows/2))
		if len(starts) > 0:
			s,e = 0,0
			if end0 + self.windows <= len(array_d1):
				e = end0 + self.windows
			else:
				e = end0

			if starts[0]-self.windows < 0:
				s = starts[0]
			else:
				s = starts[0]-self.windows
			centri.append([s,e])
			starts,ends = [],[]
		chro_s = str(genome[key].seq).upper()
		maxindex = max(index)
		# 创建图形和轴对象
		fig, ax = plt.subplots(figsize=(5, 5))
		# 绘制曲线图
		chr_d = {}
		for cen in centri:
			if key not in chr_d.keys():
				chr_d[key] = 1
				name = key+'_1'
			else:
				chr_d[key] += 1
				name = key+'_'+str(chr_d[key])
			f = open(self.workpath+self.centrigff,'a+')
			f.write('\t'.join([name,str(cen[0]),str(cen[1])])+'\n')
			f.close()
			f = open(self.workpath+self.centrifasta,'a+')
			f.write('>'+name+'\n'+chro_s[cen[0]:cen[1]]+'\n')
			f.close()
			self.plot(ax,cen[0],0,cen[1]-cen[0],maxindex*1.05,self.colors.split(',')[1],0.8,False,'chro')

		ax.plot(x, index,alpha = .5,color = self.colors.split(',')[0])
		if self.chip_seq != 'None':
			if key in chip_seq.keys():
				for zsl in chip_seq[key]:
					# print(key,chip_seq[key],zsl[0],-maxindex*0.05,zsl[2],maxindex*0.2)
					# self.plot(ax,zsl[0]-(zsl[2]/2),-maxindex*0.05,(zsl[2]*2),maxindex*0.2,self.colors.split(',')[2],0.8,False,'chro')
					self.plot(ax,zsl[0],-maxindex*0.05,zsl[2],maxindex*0.2,self.colors.split(',')[2],0.8,False,'chro')

		# 设置图形标题和坐标轴标签
		ax.set_title(key +' - predicting centromeres')
		ax.set_xlabel('Chromosome')
		ax.set_ylabel('TRF Coverage')
		# 显示图形
		fig.savefig(key +'-predicting.png',dpi = 1000)
		# fig.close()

	def get_centri(self,genome,chip_seq):

		# 检查文件是否存在
		if os.path.exists(self.workpath+self.centrigff):
			# 如果文件存在，则删除
			os.remove(self.workpath+self.centrigff)
		else:
			pass

		# 检查文件是否存在
		if os.path.exists(self.workpath+self.centrifasta):
			# 如果文件存在，则删除
			os.remove(self.workpath+self.centrifasta)
		else:
			pass
		
		lens = pd.read_csv(self.workpath+self.lens,header = None, sep='\t', comment='#').sort_values(by=[0,1],ascending= [True,True])
		lens = dict(zip(lens[0],lens[1]))
		gff = pd.read_csv(self.workpath+self.trfgff,header = None, sep='\t', comment='#').sort_values(by=[0,3],ascending= [True,True])
		gff[8] = gff[1]
		print("ALL GFF")
		gff[7] = abs(gff[4] - gff[3])
		gff[[3, 4]] = gff[[3, 4]].astype(str)
		gff[9] = '('+gff[3] +','+ gff[4]+')'
		gff[[3, 4]] = gff[[3, 4]].astype(int)
		print(gff)
		dic = gff.groupby(0).groups# 按照第几列分组
		lt = list(lens.keys())
		for key in lt:
			if key == 'Unassemble':
				continue
			self.get_overlap(gff,dic,lens,key,genome,chip_seq)


	def run(self):
		# 检查文件夹是否存在
		if not os.path.exists(self.workpaths):
			# 如果文件夹不存在，则创建文件夹
			os.makedirs(self.workpaths)
			# print("文件夹不存在，已创建:", self.workpaths)
		else:
			pass
			# print("文件夹已存在:", self.workpaths)
		# 改变工作路径
		genome = SeqIO.to_dict(SeqIO.parse(self.genome_file, "fasta"))
		if self.chip_seq != 'None':
			chip_seq = bez.read_zsl(self.chip_seq)
		else:
			chip_seq = ''
		os.chdir(self.workpath+self.workpaths)

		self.get_centri(genome,chip_seq)
