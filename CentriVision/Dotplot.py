import sys
import ast
import os
import seaborn as sns
import gc# 内存管理模块
from tqdm import trange
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import multiprocessing # Step I : 导入模块
from multiprocessing import cpu_count#读取CPU核心数用于匹配线程数
import CentriVision.bez as bez
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib
from scipy.signal import find_peaks
from matplotlib import gridspec
from Bio.Align import PairwiseAligner
from matplotlib.colors import Normalize, ListedColormap,LinearSegmentedColormap
matplotlib.use('Agg')# 强制使用非交互式后端（交互式后端在GUI资源枯竭或满线程情况下报错）



class Dotplot():
	def __init__(self, options):
		self.workpath = os.getcwd()+'/'
		self.workpaths = './Dotplot/'
		self.fastapath = './fasta/'
		self.plotpath = './plot/'
		self.plot = 'False'
		self.temp = 'False'
		self.temppath = './temp/'
		self.pilpath = './PIL/'
		self.histpath = './hist/'
		self.hist = 'False'
		self.cpu = 16
		self.gap = -10
		self.windows_size = 15
		self.windows = 4000
		self.minlength = 10
		self.konbai = 20
		self.poly = 'False'
		self.seed = 0
		self.window = 0
		self.m = 0.8
		# self.score_matrix = pd.read_excel(bez.get_scorefile())

		path = bez.get_path()
		font_path = os.path.join(path, 'example/arial.ttf')
		from matplotlib.font_manager import FontProperties
		self.font_prop = FontProperties(fname=font_path)

		bez_conf = bez.config()
		for k, v in bez_conf:# 
			setattr(self, str(k), v)
		for k, v in options:
			setattr(self, str(k), v)
			print(str(k), ' = ', v)
		self.cpu = int(self.cpu)
		self.windows = int(self.windows)
		self.windows_size = int(self.windows_size)
		self.gap = int(self.gap)
		self.minlength = int(self.minlength)
		self.konbai = int(self.konbai)
		self.m = float(self.m)
		if self.cpu > cpu_count():
			print('Warning:Thread count reduced to '+str(cpu_count())+' due to system limitations.')
			self.cpu = cpu_count()

	def update_best_matrix(self,s1,s2,best_matrix,gap):
		job = len(s2)+1
		for i in range(job):
			p = int(((i+2)/job)*30)
			print("\r["+"*"*p+" "*(30-p)+"] "+str(i+1)+"/"+str(job)+' -> Loading ',end="")
			for j in range(job):
				if i == 0 or j == 0 or i == j:
					best_matrix[i][j] = 0
				else:
					if s2[i-1] != s1[j-1] or s2[i-1] == 'N' or s1[j-1] == 'N':
						best_matrix[i][j] = 0
					else:
						match_score = best_matrix[i-1][j-1]+1
						best_matrix[i][j] = match_score

		return best_matrix

	def clear_matrix(self,s1,s2,data,gap,match,unmatch,name):
		repeat_dic = {}
		print('')
		job = len(data)
		length = job-1
		np_p =  np.zeros(shape= (job,job),dtype = int)
		for i in range(job):
			p = int(((i+1)/job)*30)
			print("\r["+"*"*p+" "*(30-p)+"] "+str(i+1)+"/"+str(job)+ ' -> Cleaning ',end="")
			for j in range(job):
				if data[-i][-j] >= self.minlength and np_p[-i][-j] == 0: #从后往前遍历，如果当前分数大于最小长度
					l = data[-i][-j]
					np_p[-i][-j] = data[-i][-j]
					x,y = i,j
					d = 0
					s1_,s2_ = s1[-i],s2[-j]
					while 1:
						note = data[-x][-y]
						note1 = data[-x-1][-y-1]
						note2 = data[-x-1][-y]
						note3 = data[-x][-y-1]
						if -x-1 > -job and -y-1 > -job and note1 != 0 and (note == note1 + match or note == note1 + unmatch):
							np_p[-x-1][-y-1] = data[-x-1][-y-1]
							s1_,s2_ = s1[-x-1]+s1_,s2[-y-1]+s2_
							x,y = -(-x-1),-(-y-1)
							d += 1
							continue
						else:
							break
					if s1_ not in repeat_dic.keys():
						repeat_dic[s1_] = [length-i-l+1,length-j-l+1]# 输出的gff文件中的位置从0开始
					else:
						repeat_dic[s1_] += [length-i-l+1,length-j-l+1]# 输出的gff文件中的位置从0开始
						repeat_dic[s1_] = list(set(repeat_dic[s1_]))
		out_list = []
		all_repeat = []
		for key in repeat_dic.keys():
			all_repeat.append(len(repeat_dic[key]))
			s = '\t'.join([name,str(len(s1)),key,str(len(key)),str(len(repeat_dic[key])),'_'.join([str(i) for i in repeat_dic[key]])])
			out_list.append(s)
		f0 = open(self.workpaths+self.fastapath+self.outfile+'.repeat.index.csv','a')
		f0.write('\n'.join(out_list)+'\n')
		f0.close()

		return np_p,len(data)-1

	def np2plot(self,s1,s2,best_matrix,gap,name):
		print('')
		np_p =  np.zeros(shape= (len(s2)+1,len(s1)+1),dtype = int)
		np_p[np_p==0] = 255
		x,y = [],[]
		job = len(s2)+1
		for i in range(job):
			p = int(((i+1)/job)*30)
			print("\r["+"*"*p+" "*(30-p)+"] "+str(i+1)+"/"+str(job)+ ' -> Drawing ',end="")
			for j in range(job):
				if best_matrix[i][j] > 0 or i == j:
					x.append(j)
					y.append(i)
					np_p[i][j] = 0
		if self.plot == 'True':
			fig = plt.subplots(figsize = (12,12))
			plt.scatter(x, y,s=0.4,c='black', edgecolor='none')
			plt.xlim(0,job)
			plt.ylim(0,job)

			plt.gca().invert_yaxis()
			plt.title(name, fontsize=12, fontproperties=self.font_prop)# , fontsize=12, fontproperties=self.font_prop
			plt.savefig(self.workpaths+'./plot/'+name+'.png',dpi = 1000, bbox_inches='tight')	
			# 清除画布
			plt.close()

		# plt.imshow(np.array(np_p), cmap='gray')
		# plt.savefig(self.workpaths+'./plot/'+name+'.hist.png', dpi = 1000, bbox_inches = 'tight')# 存储图片
		# plt.close()

	def condense_array(self,original_array, target_length=100):
		old_length = len(original_array)
		new_original_array = []
		for i in original_array:
			new_original_array = new_original_array + [i]*target_length
		condensed_array = []
		for i in range(target_length):
			condensed_array.append(sum(new_original_array[i*old_length:(i+1)*old_length])/float(old_length))
		return np.array(condensed_array)

	def read_best_matrix(self,matrix):
		matrix = np.int64(matrix>0)
		# 对每行求和
		row_sums = np.sum(matrix, axis=1)[1:]#行
		col_sums = np.sum(matrix, axis=0)[1:]#列
		# 将数组浓缩为长度为100的数组
		row_sums = list(self.condense_array(row_sums))
		col_sums = list(self.condense_array(col_sums))
		matrix_s = row_sums+col_sums
		# del matrix,row_sums,col_sums
		return [float(i)/(len(matrix)-1) for i in matrix_s]

	def find_best_split(self,arro):
		# 对数组进行排序
		arr = sorted(arro)
		# 计算数组长度
		n = len(arr)
		# 初始化变量
		min_difference = float('inf')
		best_split_index = 0
		# 初始化前缀和数组
		prefix_sum = np.cumsum(arr)
		# 遍历每一个可能的分割点
		for i in range(1, n):
			# 计算前半部分和后半部分的平均值
			left_sum = prefix_sum[i - 1]
			right_sum = prefix_sum[-1] - prefix_sum[i - 1]
			left_avg = left_sum / i
			right_avg = right_sum / (n - i)
			# 计算两个部分的平均值的差异
			difference = abs(left_avg - right_avg)
			# 更新最小差异和最佳分割点
			if difference < min_difference:
				min_difference = difference
				best_split_index = i
		# 根据最佳分割点分割数组
		left_part = arr[:best_split_index]
		right_part = arr[best_split_index:]
		return left_part, right_part

	def levenshtein(self,s1, s2):
		if len(s1) < len(s2):
			return levenshtein(s2, s1)
		if len(s2) == 0:
			return len(s1)
		previous_row = list(range(len(s2) + 1))
		for i, c1 in enumerate(s1):
			current_row = [i + 1]
			for j, c2 in enumerate(s2):
				insertions = previous_row[j + 1] + 1
				deletions = current_row[j] + 1
				substitutions = previous_row[j] + (c1 != c2)
				current_row.append(min(insertions, deletions, substitutions))
			previous_row = current_row
		return previous_row[-1]

	def AMPD(self, data, threshold):
		"""
		实现AMPD算法，并通过阈值过滤小波动
		:param data: 1-D numpy.ndarray 
		:param threshold: 过滤波峰的振幅阈值，默认为0（即不过滤）
		:return: 波峰所在索引值的列表
		"""
		p_data = np.zeros_like(data, dtype=np.int32)
		count = data.shape[0]
		arr_rowsum = []
		# 计算行和
		for k in range(1, count // 2 + 1):
			row_sum = 0
			for i in range(k, count - k):
				if data[i] > data[i - k] and data[i] > data[i + k]:
					row_sum -= 1
			arr_rowsum.append(row_sum)
		min_index = np.argmin(arr_rowsum)
		max_window_length = min_index
		# 统计局部最大值
		for k in range(1, max_window_length + 1):
			for i in range(k, count - k):
				if data[i] > data[i - k] and data[i] > data[i + k]:
					p_data[i] += 1

		# 获取波峰索引
		peak_indices = np.where(p_data == max_window_length)[0]
		if len(peak_indices) == 0:
			peak_indices = np.where(p_data == max(p_data))[0]
		# 通过阈值过滤小波峰
		if threshold > 0:
			peak_indices = [i for i in peak_indices if data[i] > threshold]
		return np.array(peak_indices)

	def start(self,s1,s2,name):
		length = len(s1)
		best_matrix = np.zeros(shape= (length+1,length+1),dtype = int)	#初始化得分矩阵
		best_matrix = self.update_best_matrix(s1,s2,best_matrix,self.gap)
		gc.collect()
		if self.temp == 'True':
			np.savetxt(self.workpaths+'./temp/'+name+'.csv',best_matrix,fmt='%d',delimiter='\t')
		best_matrix1,length = self.clear_matrix(s1,s2,best_matrix,self.gap,1,0,name)
		del best_matrix
		gc.collect()
		matrix_s = self.read_best_matrix(best_matrix1)
		self.np2plot(s1,s2,best_matrix1,self.gap,name)
		for i in range(len(best_matrix1)):
			best_matrix1[i][i] = 1
		print('\t'+name)
		matrix0 = np.where(np.array(best_matrix1) >= 1, 1, 0)
		if self.temp == 'True':
			np.savetxt(self.workpaths+'./temp/'+name+'.one.csv',matrix0,fmt='%d',delimiter='\t')
		# 神来之笔
		matrix_new = []
		for y in best_matrix1:
			arr = np.array(y)
			if arr.sum() == 0:
				continue
			# 找到第一个非零数字的位置
			nonzero_indices = np.where(arr != 0)  # 获取非零索引
			first_nonzero_index = nonzero_indices[0][0]  # 获取第一个非零索引
			if first_nonzero_index != 0:
				y0_ = list(y)[first_nonzero_index:] + list(y)[:first_nonzero_index]
			else:
				y0_ = list(y)
			matrix_new.append(y0_)
		matrix = np.where(np.array(matrix_new) >= 1, 1, 0)# 冗余的一次数据清洗
		col_sums = np.sum(matrix, axis=0)[1:]#列
		x = [i for i in range(len(col_sums))]
		x1 = np.array(x)  # 示例索引数组
		y1 = np.array(col_sums)  # 示例数据，替换为你的实际数据
		y1 = np.maximum(y1 - (max(y1)*0.05), 0) # 去除部分背景

		fig = plt.figure(figsize=(6, 7))
		gs = gridspec.GridSpec(6, 7)
		ax1 = fig.add_subplot(gs[0:1,0:6])# , aspect=1
		px = self.AMPD(y1,max(y1)*0.25)
		peakout = ''
		segments_split1,segments_split2 = 0,0
		if len(px) == 0:
			print(name,"The peak number is 0.")
			px_c,mode_value,number,block,maxpeak = 0,0,0,0,0
			peakx1,peaky1 = [],[]
			peakx0,peaky0 = [],[]
		else:
			peakx,peaky = px, y1[px]
			peakx1,peaky1 = [],[]
			peakx0,peaky0 = [],[]
			maxpeak = max(peaky)# 最相似的序列
			if maxpeak < length*0.15:
				print(name,"Peak value too low.")
				px_c,mode_value,number,block,maxpeak = 0,0,0,0,0
			else:
				for i in range(len(peakx)):
					peakx0.append(peakx[i])
					peaky0.append(peaky[i])
					if peaky[i] < maxpeak*self.m:
						continue
					else:
						peakx1.append(peakx[i])
						peaky1.append(peaky[i])
				if len(peakx1) > 5:
					y_array = np.array(peaky1)
					top_indices = np.argsort(y_array)[-5:]
					top_x = [peakx1[i] for i in top_indices]
					top_y = [peaky1[i] for i in top_indices]
					px_c = sum(top_y)/len(top_y)
				else:
					top_x = peakx1
					top_y = peaky1
					px_c = sum(top_y)/len(top_y)

				# 使用嵌套循环计算差值
				differences = []
				for i in range(len(peakx0)):
					if i+2 >= len(peakx0):
						continue
					for j in range(i + 1, i + 2):
						seedlength = abs(peakx0[i] - peakx0[j])
						if seedlength >= self.minlength:
							differences.append(seedlength)  # 计算绝对差值

				if len(differences) >= 3:
					from statistics import mode,mean
					# 计算数组的众数
					mode_value = mode(differences)
					# 计算数组的均值
					mean_value = mean(differences)
					# number = int(len(s1)/mode_value)
				elif len(differences) == 0:# 若黑点也只有1或者2时
					print(name,"The peak count is less than 3.")
					mode_value = 0
					# number = 0
				else:
					mode_value = min(differences)
					# number = int(len(s1)/mode_value)
				block = len(peakx1) + 1

			if len(peakx0)>0:
				# 处理分割点（包含首尾）
				split1 = [0] + sorted(peakx0) + [length]
				# 计算第一个数组分割的片段数
				segments_split1 = len(split1) - 1
				if len(peakx1) > 0:
					split2 = [0] + sorted(peakx1) + [length]
					# 计算第二个数组分割的片段数
					segments_split2 = len(split2) - 1
					# 统计每个 split2 区间包含的 split1 片段数
					contained_segments = []
					for i in range(len(split2) - 1):
						start, end = split2[i], split2[i+1]
						# 找到 split1 中属于 [start, end] 的分割点
						sub_points = [p for p in split1 if start <= p <= end]
						# 区间数 = 分割点数 - 1
						contained_segments.append(len(sub_points) - 1)
				else:
					segments_split2 = 0
					contained_segments = []
				peakout = '|'.join([str(segments_split1),str(segments_split2),'/'.join([str(i) for i in contained_segments])])

		ax1.plot(x1, y1, alpha=0.6, label='Phase synchronization', zorder=1)
		ax1.scatter(peakx0,peaky0, color="black",s=1, label='AMPD', alpha=0.8, zorder=3)
		ax1.scatter(peakx1, peaky1, color="red",marker='x',s=10, label='AMPD-F', alpha=0.8, zorder=2)
		ax1.axhline(y=(max(y1)*0.25), color='green', linestyle='-', linewidth=1, zorder=4)
		ax1.axhline(y=(maxpeak*self.m), color='red', linestyle='-', linewidth=1, zorder=4)
		ax1.set_xlim(0, len(s1)+1)
		ax1.grid(True)

		ax2 = fig.add_subplot(gs[1:7,0:6])# , aspect=1

		matrix = np.where(best_matrix1 >= 1, 1, 0)# 冗余的一次数据清洗
		matrix[matrix==0] = 255
		ax2.imshow(matrix, cmap='gray')

		plt.title(name, fontsize=12, fontproperties=self.font_prop)
		plt.savefig(self.workpaths+'./hist/'+name+'.png', dpi = 600, bbox_inches = 'tight')# 存储图片

		plt.close()
		gc.collect()
		GC = (s1.count('G')+s1.count('C'))/len(s1)
		seq_mode = mode_value

		def find_and_extract_sequence(data,length):
			# 将序列转换为字符串
			sequence_length = len(data)
			repeat_start = -1
			found = False
			# 查找第一个至少包含3个连续的"A"或"T"的单碱基重复
			for i in range(sequence_length - 2):
				if (data[i:i+3] == 'AAA' or data[i:i+3] == 'TTT'):
					repeat_start = i  # 找到第一个重复的位置
					while i + 3 < sequence_length and data[i] == data[i+3]:
						i += 1
					found = True
					break
			if not found:
				i = -1
			# 提取从该位置后的100个碱基
			start_pos = i + 1  # 从重复后的下一个碱基开始
			if start_pos + length <= sequence_length:
				# 如果后面有足够的碱基，直接提取
				extracted_sequence = data[start_pos:start_pos + length]
			else:
				# 如果后面的碱基不足100，则从序列前端补充
				missing_length = length - (sequence_length - start_pos)
				extracted_sequence = data[start_pos:]+data[:missing_length]
			return extracted_sequence
		seq0 = find_and_extract_sequence(s1.upper(),seq_mode)
		# return [length,GC,number,px_c,seq_mode,block]+matrix_s+[seq0,peakout]
		return [length,GC,segments_split1,px_c,seq_mode,segments_split2]+matrix_s+[seq0,peakout]

	def merge_single_base_repeats(self,sequence):
		if not sequence:
			return ""
		merged_sequence = sequence[0]
		
		for base in sequence[1:]:
			if base != merged_sequence[-1]:
				merged_sequence += base
		return merged_sequence

	def run(self):
		pool = multiprocessing.Pool(processes = self.cpu)
		genome = SeqIO.to_dict(SeqIO.parse(self.genome_file, "fasta"))
		if not os.path.exists(self.workpaths):
			# 如果文件夹不存在，则创建文件夹
			os.makedirs(self.workpaths)
		if not os.path.exists(self.workpaths+self.fastapath):
			os.makedirs(self.workpaths+self.fastapath)
		if not os.path.exists(self.workpaths+self.plotpath) and self.plot == 'True':
			os.makedirs(self.workpaths+self.plotpath)

		if not os.path.exists(self.workpaths+self.temppath) and self.temp == 'True':
			os.makedirs(self.workpaths+self.temppath)

		# if not os.path.exists(self.workpaths+self.pilpath):
		# 	os.makedirs(self.workpaths+self.pilpath)
		if not os.path.exists(self.workpaths+self.histpath):
			os.makedirs(self.workpaths+self.histpath)

		# if not os.path.exists(self.workpaths+'./plotcsv/'):
		# 	os.makedirs(self.workpaths+'./plotcsv/')
		# if not os.path.exists(self.workpaths+'./clearcsv/'):
		# 	os.makedirs(self.workpaths+'./clearcsv/')
		f0 = open(self.workpaths+self.fastapath+self.outfile+'.repeat.index.csv','w')
		f0.write('\t'.join(['split_name','split_length','seq','length','number','index(Index starts from 0)'])+'\n')
		f0.close()
		f = open(self.workpaths+self.fastapath+self.outfile+'.split.fasta','w')
		f1 = open(self.workpaths+self.fastapath+self.outfile+'.split.gff','w')
		f2 = open(self.workpaths+self.fastapath+self.outfile+'.split.lens','w')
		results = {}
		for key in genome.keys():
			s = str(genome[key].seq).upper()
			job = int(len(s)/self.windows)
			if job == 0:
				job = 1
			f2.write('\t'.join([key,str(len(s)),str(job)])+'\n')
			for i in range(job):
				if i == job - 1:
					s0 = s[self.windows*i:]
				else:
					s0 = s[self.windows*i:self.windows*(i+1)]
				if self.poly.upper() == 'TRUE': 
					s0 = self.merge_single_base_repeats(s0)
				s1,s2 = s0,s0
				if 'N' in s1:
					continue
				name = key+'_s'+str(i)
				f.write('>'+name+'\n'+s1+'\n')		
				f1.write('\t'.join([key,'CentriVision',name,str(self.windows*i),str(self.windows*(i+1)),str(len(s1)),'+',str(i+1),'',name])+'\n')	
				dic = {}
				dic['gff'] = [key,'dotplot',name,str(self.windows*i+1),'.']
				dic['value'] = pool.apply_async(self.start, (s1,s2,name),  error_callback=bez.print_error)
				results[name] = dic
		pool.close() # Step IV : 准备结束
		pool.join() # Step IV : 完全结束
		f.close()
		f1.close()
		f2.close()
		print("********************输出********************")

		names = list(results.keys())
		maker = ['length','GC','NM(Number of monomers)','top5','LM(Length of monomers)','N_HOR(Number of HORs)']+['M'+str(i) for i in range(1,201)]+['monomer','peakx']
		f = open(self.workpath+self.outfile,'w')
		f.write('\t'.join(['centromere','dotplot','name','start','end']+maker)+'\n')
		value = [[] for _ in range(len(maker))]
		for j in range(len(names)):
			name = names[j]
			lt = results[name]['value'].get()
			if lt == None:
				continue
			for i in range(len(lt)):
				value[i].append(lt[i])
			lt = results[name]['gff'] + lt
			f.write('\t'.join([str(s) for s in lt])+'\n')
		f.close()
