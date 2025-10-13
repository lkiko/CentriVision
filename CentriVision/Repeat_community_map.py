import re
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool
from multiprocessing import cpu_count#读取CPU核心数用于匹配线程数
import CentriVision.bez as bez
import os
from matplotlib import gridspec
import seaborn as sns
from matplotlib.colors import Normalize, ListedColormap,LinearSegmentedColormap

class Repeat_community_map():
	def __init__(self, options):
		self.workpath = os.getcwd()+'/'
		self.top = '10'
		self.dpi = '1000'
		self.windows = '50000'
		self.step = '5000'
		self.col = '2'
		self.gradient = 'seismic'
		self.colors = '#824880,#00a381,#e6b422,#96514d,#bf242a,#cc7eb1,#00a497,#bc763c,#f09199,#007bbb,#ffea00,#007b43,#9079ad,#17184b,#f39800,#ea5506,#474a4d,#00552e'
		for k, v in options:
			setattr(self, str(k), v)
			print(str(k), ' = ', v)
		self.top =int(self.top)
		self.colors = self.colors.split(',')
		self.dpi = int(self.dpi)
		self.col = int(self.col)
		self.row = 0
		self.windows = int(self.windows)
		self.step = int(self.step)

	def value_to_color(self,value, cmap, vmin, vmax):
		norm = Normalize(vmin=vmin, vmax=vmax)
		return cmap(norm(value))

	def get_index(self,a,b):
		return [min(a,b)-1,max(a,b)]

	def covarage(self,suf_dic,length):
		# print(suf_dic)
		array_d = np.zeros(length, dtype = int)
		for start,end in suf_dic:
			start = start
			array_d[start:end] = array_d[start:end]+1
		array_d1 = np.int64(array_d>0)
		x1,y1,y_max = self.line_plot(array_d1,length)
		return x1,y1,y_max

	def line_plot(self,data,length):
		x,y = [],[]
		start,end = -self.step,self.windows-self.step
		while end < length:
			start,end = start+self.step,end+self.step
			index = int((start+end)/2)
			data0 = 100*sum(data[start:end])/self.windows
			x.append(index)
			y.append(data0)
		# print(y)
		y_max = max(y)
		# y = list(map(lambda a: (a/y_max)*h+h0, y))
		return x,y,y_max

	def read_focus(self,file):
		focus = {}
		for line in open(file,'r'):
			lt = line.strip('\n').split('\t')
			lt[0] = lt[0]
			if lt[0] == 'Chr_ID':
				continue
			if lt[0] not in focus.keys():
				focus[lt[0]] = []
				focus[lt[0]].append([int(lt[1]),int(lt[2]),abs(int(lt[1])-int(lt[2]))+1])
			else:
				focus[lt[0]].append([int(lt[1]),int(lt[2]),abs(int(lt[1])-int(lt[2]))+1])
		return focus

	def plot(self,ax,x,y,w,h,c,a,p,n):
		rect = plt.Rectangle((x,y),w,h, facecolor=c,alpha=a, edgecolor='none')# x,y,宽，高
		ax.add_patch(rect)

	def run(self):
		if self.focus_areas == 'None':
			focus = {}
		else:
			lt = self.focus_name.split(',')
			ltfile = self.focus_areas.split(',')
			focus = {}
			focusc = {}
			for names,file in zip(lt,ltfile):
				name = names.split(':')[0]
				focus[name] = self.read_focus(file)
				focusc[name] = names.split(':')[1]

		lens,chr_list = bez.read_lens(self.lens)# 染色体字典
		genegff = pd.read_csv(self.repeat_gff,header = None, sep='\t', comment='#').sort_values(by=[0,3],ascending= [True,True])
		genegff[10] = genegff.apply(lambda x: self.get_index(x[3],x[4]), axis=1)
		gff = dict(zip(genegff[9],genegff[0]))
		gene = dict(zip(genegff[9],genegff[10]))
		genechr = dict(zip(genegff[9],genegff[0]))
		community = {}
		communitys = []
		color = {}

		if self.model == 'global':
			if self.top > len(self.colors):
				# 创建渐变颜色映射
				colors = sns.color_palette(self.gradient, as_cmap=True)
				cmap = LinearSegmentedColormap.from_list("Custom", colors(np.linspace(0, 1, 256)))
				self.colors = []
				for i in range(0,100,int(100/self.top)):
					self.colors.append(self.value_to_color(i,cmap,0,100))
			communityf = open(self.community_file,'r').readlines()[:self.top]
			for i,line in zip(range(self.top), communityf):
				lt = line.strip('\n').split('\t')
				community[lt[0]] = lt[1:]
				color[lt[0]] = self.colors[i]
				communitys.append(lt[0])
		else:
			top_d = {}
			communityf = open(self.community_file,'r').readlines()
			chr_com = {}
			communityd = {}
			for line in communityf:
				lt = line.strip('\n').split('\t')
				communityd[lt[0]] = lt[1:]
				for i in lt[1:]:
					if gff[i] not in chr_com.keys():
						dic = {lt[0]:1}
						chr_com[gff[i]] = dic
					else:
						if lt[0] not in chr_com[gff[i]].keys():
							chr_com[gff[i]][lt[0]] = 1
						else:
							chr_com[gff[i]][lt[0]] += 1
			for chro in (chro for chro in chr_list if chro in chr_com.keys()):
				dic = chr_com[chro]
				top = dict(sorted(dic.items(), key=lambda item: item[1], reverse=True)[:self.top])
				top_d[chro] = top
				# print(top)
				# if 'CENT10' in top.keys():
				# 	continue
				# else:
				# 	top['CENT10'] = dic['CENT10']
				# print(top)
				for i in top.keys():
					communitys.append(i)
					community[i] = communityd[i]
			communitys = list(set(communitys))
			# print(len(communitys))
			if len(communitys) > len(self.colors):
				# 创建渐变颜色映射
				# colors = sns.color_palette("rainbow", as_cmap=True)
				colors = sns.color_palette(self.gradient, as_cmap=True)
				cmap = LinearSegmentedColormap.from_list("Custom", colors(np.linspace(0, 1, 256)))
				self.colors = []
				for i in range(0,100,int(100/len(communitys))):
					self.colors.append(self.value_to_color(i,cmap,0,100))
			if ':' in self.colors[0]:
				for i in self.colors:
					color[i.split(':')[0]] = i.split(':')[1]
			else:
				for i in range(len(communitys)):
					color[communitys[i]] = self.colors[i]

		if len(chr_list)%self.col == 0:
			rows = int(len(chr_list)/self.col)
		else:
			rows = int(len(chr_list)/self.col) + 1
		fig = plt.figure(figsize=(20, 2*rows))
		gs = gridspec.GridSpec(5*rows, 30*self.col)#, width_ratios=[10, 1]
		data = {}
		for i in communitys:
			dic = {}
			for repeat in community[i]:
				if genechr[repeat] not in dic.keys():
					dic[genechr[repeat]] = [gene[repeat]]
				else:
					dic[genechr[repeat]].append(gene[repeat])
			for chro in chr_list:
				if chro not in dic.keys():
					continue
				if self.model != 'global' and i not in top_d[chro]:
					continue
				x,y,ymax = self.covarage(dic[chro],lens[chro]['end'])
				if chro not in data.keys():
					data[chro] = [[x,y,ymax,i]]
				else:
					data[chro].append([x,y,ymax,i])
		for j in range(len(chr_list)):
			row = int(j/self.col)
			col = j%self.col
			print(row,col,(row*5)+2,(row+1)*5,30*col,30*(col+1))
			chro = chr_list[j]
			if chro not in data.keys():
				continue
			if col == 0:
				ax = fig.add_subplot(gs[(row*5)+2:(row+1)*5,30*col:30*(col+1)-2])
			else:
				ax = fig.add_subplot(gs[(row*5)+2:(row+1)*5,30*col+2:30*(col+1)])

			for x,y,ymax,i in data[chro]:
				# print(chro,x,y)
				ax.plot(x,y, color = color[i],linestyle = '-', linewidth=.5,label=i)
			for name in focus.keys():
				if chro in focus[name].keys():
					areas = focus[name][chro]
					# print(areas)
					for area in areas:
						start,end = min(area[0],area[1]),max(area[0],area[1])
						# start0 = start/length
						w = abs(end-start)
						self.plot(ax,start,-5,w,5,focusc[name],0.8,False,chro)
						# plot_(ax,start0,h,w,0.3,'red',0.8,False,chro)

			# 给子图添加坐标轴名字
			# ax.set_xlabel("X轴名称")
			# ax.set_ylabel(chro.split('hifi_')[1][:-2],rotation=90, fontsize=16)
			ax.set_ylabel(chro,rotation=90, fontsize=16)
			# 设置图例及其字体大小
			ax.legend(fontsize=10)  # 设置图例的字体大小
			# 设置刻度文字大小
			# 设置纵轴刻度
			ax.set_yticks([0, 50, 100])

			ax.tick_params(axis='both', which='major', labelsize=14)

		plt.savefig(self.savefile, dpi = self.dpi, bbox_inches = 'tight')
		exit()


