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
		self.windows = 4000
		self.minlength = 8
		self.konbai = 20
		self.poly = 'False'
		self.seed = 0
		self.window = 0
		self.m = '0.3'
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
		self.gap = int(self.gap)
		self.minlength = int(self.minlength)
		self.konbai = int(self.konbai)
		self.m = float(self.m)
		if self.cpu > cpu_count():
			self.cpu = cpu_count()

	# def get_match_score(self,s1,s2,score_matrix):	#获取替换打分
	# 	# global  score_matrix
	# 	s = ['A','G','C','T'].index(s2)
	# 	score = score_matrix[s1][s]
	# 	return score

	def update_best_matrix(self,s1,s2,best_matrix,gap):
		job = len(s2)+1
		for i in range(len(s2)+1):
			p = int(((i+2)/job)*30)
			print("\r["+"*"*p+" "*(30-p)+"]  \033[0;31;42m"+str(i+1)+"\033[0m / "+str(job)+' -> update_matrix ',end="")
			for j in range(len(s1)+1):
				if i == 0 or j == 0 or i == j:
					best_matrix[i][j] = 0
				else:
					if s2[i-1] != s1[j-1] or s2[i-1] == 'N' or s1[j-1] == 'N':
						best_matrix[i][j] = 0
					else:
						match = 1
						if j >= 2 and i >= 2 and self.poly == 'True':# 舍弃单碱基重复
							if s2[i-2] == s1[j-2] and s2[i-1] == s2[i-2]:
								match = 0
						match_score = best_matrix[i-1][j-1]+match
						best_matrix[i][j] = match_score
		# print(best_matrix)
		return best_matrix

	def clear_matrix(self,s1,s2,data,gap,match,unmatch,name):
		cleardata = []
		repeat_dic = {}
		print('')
		job = len(data)
		length = len(s1)
		np_p =  np.zeros(shape= (job,job),dtype = int)
		for i in range(job):
			p = int(((i+1)/job)*30)
			print("\r["+"*"*p+" "*(30-p)+"]  \033[0;31;42m"+str(i+1)+"\033[0m / "+str(job)+ ' -> clear ',end="")
			for j in range(job):
				if data[-i][-j] >= self.minlength and np_p[-i][-j] == 0: #从后往前遍历，如果当前分数大于最小长度
					# print('score',data[-i][-j])
					l = data[-i][-j]
					cleardata.append(data[-i][-j])
					np_p[-i][-j] = data[-i][-j]
					x,y = i,j
					d = 0
					s1_,s2_ = s1[-i],s2[-j]
					while 1:
						note = data[-x][-y]
						note1 = data[-x-1][-y-1]
						note2 = data[-x-1][-y]
						note3 = data[-x][-y-1]
						# print('\n',note,note1,note2,note3,x,y)
						if -x-1 > -job and -y-1 > -job and note1 != 0 and (note == note1 + match or note == note1 + unmatch):
							np_p[-x-1][-y-1] = data[-x-1][-y-1]
							# print(s1[-x-1],s2[-y-1],data[-x-1][-y-1])
							s1_,s2_ = s1[-x-1]+s1_,s2[-y-1]+s2_
							x,y = -(-x-1),-(-y-1)
							d += 1
							continue
						else:
							break
							# exit()
					# print(s1_,s2_)
					if s1_ not in repeat_dic.keys():
						repeat_dic[s1_] = [length-i-l+1,length-j-l+1]# 输出的gff文件中的位置从0开始
					else:
						repeat_dic[s1_] += [length-i-l+1,length-j-l+1]# 输出的gff文件中的位置从0开始
						repeat_dic[s1_] = list(set(repeat_dic[s1_]))
		# print(repeat_dic)
		out_list = []
		all_repeat = []

		for key in repeat_dic.keys():
			all_repeat.append(len(repeat_dic[key]))
			s = '\t'.join([name,str(len(s1)),key,str(len(key)),str(len(repeat_dic[key])),'_'.join([str(i) for i in repeat_dic[key]])])
			out_list.append(s)

		# seed_number = len(repeat_dic)
		# means = sum(all_repeat)/seed_number

		f0 = open(self.workpaths+self.fastapath+self.outfile+'.repeat.index.csv','a')
		f0.write('\n'.join(out_list)+'\n')
		f0.close()


		if len(all_repeat) > 10:# 防止数据集为空
			from statistics import mode,median,mean,variance
			# 计算数组的众数
			mode_value = mode(all_repeat)
			# 计算数组的均值
			mean_value = mean(all_repeat)
			# 计算数组的中位数
			median_value = median(all_repeat)
			# 计算数组的方差
			variance_value = variance(all_repeat)
			max_count = max(all_repeat)
		else:
			mode_value,mean_value,median_value,variance_value,max_count = 0,0,0,0,0

		return np_p,len(data)-1,mode_value,mean_value,median_value,variance_value,max_count,cleardata

	def np2plot(self,s1,s2,best_matrix,gap,name):
		print('')
		np_p =  np.zeros(shape= (len(s2)+1,len(s1)+1),dtype = int)
		np_p[np_p==0] = 255
		x,y = [],[]
		job = len(s2)+1
		for i in range(len(s2)+1):
			p = int(((i+1)/job)*30)
			print("\r["+"*"*p+" "*(30-p)+"]  \033[0;31;42m"+str(i+1)+"\033[0m / "+str(job)+ ' -> drawing ',end="")
			for j in range(len(s1)+1):
				if best_matrix[i][j] > 0 or i == j:
					x.append(j)
					y.append(-i)
					np_p[i][j] = 0
		if self.plot == 'True':
			fig = plt.subplots(figsize = (12,12))
			plt.scatter(x, y,s=0.1,c='black', edgecolor='none')
			plt.xlim(0,job)
			plt.ylim(-job,0)
			# plt.axis('off')
			plt.savefig(self.workpaths+'./plot/'+name+'.png',dpi = 1000, bbox_inches='tight')	
			# 清除画布
			plt.close()

		plt.imshow(np.array(np_p), cmap='gray')
		plt.savefig(self.workpaths+'./PIL/'+name+'.png', dpi = 1000, bbox_inches = 'tight')# 存储图片
		plt.close()

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

	def start(self,s1,s2,name):
		windows = [s1[i:i+30] for i in range(len(s1) - 29)]
		y = []
		for i in range(len(windows) - 1):
			y.append(self.levenshtein(windows[0], windows[i+1]))
		x = list(range(len(y)))
		# # x 是索引数组，y 是编辑距离数组
		x0 = np.array(x)  # 示例索引数组
		y0 = np.array(y)  # 示例数据，替换为你的实际数据
		distances = y0

		# 计算离群值（IQR 方法）
		q1, q3 = np.percentile(distances, [20, 80])
		iqr = q3 - q1
		threshold = q1 - 1.5 * iqr  # 低于此值为离群的小距离
		outlier_mask = distances < threshold
		outlier_indices = np.where(outlier_mask)[0]  # 获取索引

		# plt.figure(figsize=(10, 5))
		# plt.plot(outlier_indices, 30-distances[outlier_mask], color='blue', label="Normal distance", alpha=0.7)
		# plt.axhline(30-threshold, color='gray', linestyle='--', label=f"Outlier threshold: {threshold:.2f}")

		# 计算离群值（IQR 方法）
		q1s,q3s = 0.1, 99.9
		q1, q3 = np.percentile(distances, [q1s,q3s])
		outlier_mask_ = distances < q1
		outlier_indices = np.where(outlier_mask_)[0]  # 获取索引
		q = 0.1
		for i in range(20):
			q1s,q3s = q1s+q,q3s-q
			if len(outlier_indices) <= 2 and q1s < 10:
				q1, q3 = np.percentile(distances, [q1s,q3s])
				outlier_mask_ = distances < q1
				outlier_indices = np.where(outlier_mask_)[0]
			else:
				outlier_indices = np.where(outlier_mask_)[0]  # 获取索引
				break

		# plt.scatter(outlier_indices, 30-distances[outlier_mask_], color='green', label="Normal distance", alpha=0.7)
		# plt.axhline(30-q1, color='red', linestyle='--', label=f"Outlier threshold: {q1:.2f}")

		# plt.xlabel("Index")
		# plt.ylabel("distances")
		# plt.legend()
		# plt.savefig(self.workpaths+'./hist/'+name+'score.png', dpi = 600, bbox_inches = 'tight')# 存储图片

		block = len(outlier_indices)

		# 设置参数：height=None 表示不限制高度，distance 表示峰之间的最小距离
		min_peaks, _ = find_peaks(-y0, distance=self.minlength)
		min_peaks = min_peaks[y0[min_peaks] <= threshold]
		peakx1 = x0[min_peaks]
		# 使用嵌套循环计算差值
		differences = []
		for i in range(len(peakx1)):
			if i+2 >= len(peakx1):
				continue
			for j in range(i + 1, i + 2):
				seedlength = abs(peakx1[i] - peakx1[j])
				if seedlength >= self.minlength:
					differences.append(seedlength)  # 计算绝对差值
		if len(differences) >= 2:
			from statistics import mode,mean
			# 计算数组的众数
			mode_value = mode(differences)
			# 计算数组的均值
			mean_value = mean(differences)
			number = differences.count(mode_value)
			for i in range(1,int(mode_value*0.1)+1):
				number += differences.count(mode_value+i)
				number += differences.count(mode_value-i)
			number += 2
		else:
			mode_value = 0
			number = 0

		def dynamic_variable(p):
			if p <= 10:
				return 3
			elif p <= 50:
				return 3
			elif p >= 100:
				return 10
			else:
				# 线性插值计算
				return 3 + (p - 50) * (10 - 3) / (100 - 50)

		# self.minlength = dynamic_variable(mode_value)

		best_matrix = np.zeros(shape= (len(s2)+1,len(s1)+1),dtype = int)	#初始化得分矩阵
		best_matrix = self.update_best_matrix(s1,s2,best_matrix,self.gap)
		gc.collect()
		if self.temp == 'True':
			np.savetxt(self.workpaths+'./temp/'+name+'.csv',best_matrix,fmt='%d',delimiter='\t')
		best_matrix1,length,mode_value1,mean_value1,median_value1,variance_value1,max_count,cleardata = self.clear_matrix(s1,s2,best_matrix,self.gap,1,0,name)
		gc.collect()
		matrix_s = self.read_best_matrix(best_matrix1)
		self.np2plot(s1,s2,best_matrix1,self.gap,name)
		for i in range(len(best_matrix1)):
			best_matrix1[i][i] = 1
		print('\t'+name)
		matrix0 = np.where(np.array(best_matrix1) >= 1, 1, 0)
		if self.temp == 'True':
			np.savetxt(self.workpaths+'./temp/'+name+'.one.csv',matrix0,fmt='%d',delimiter='\t')
		matrix0[matrix0==0] = 255

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

		if self.temp == 'True':
			np.savetxt(self.workpaths+'./temp/'+name+'.phase.csv',matrix,fmt='%d',delimiter='\t')

		matrix_p = matrix.copy()
		matrix_p[matrix_p==0] = 255

		fig = plt.figure(figsize=(6, 7))
		gs = gridspec.GridSpec(6, 7)
		ax1 = fig.add_subplot(gs[0:1,0:6])# , aspect=1
		ax1.plot(x0, y0, label='Edit Distance')
		# 标记朝下的峰
		ax1.scatter(x0[min_peaks], y0[min_peaks], color='red', zorder=5, label='Downward Peaks',s=5)

		outlier_indicesx = x0[outlier_mask_]
		outlier_indicesy = y0[outlier_mask_]
		ax1.scatter(outlier_indicesx, outlier_indicesy, color='green', zorder=5, label='Downward Peaks',s=6)
		ax1.set_xlim(0, len(s1)+1)
		# 设置第一个子图的 x 轴标题
		ax1.grid(True)

		col_sums = np.sum(matrix, axis=0)[1:]#列
		x = [i for i in range(len(col_sums))]
		ax2 = fig.add_subplot(gs[1:7,0:6])# , aspect=1
		ax2.imshow(matrix_p, cmap='gray')
		ax2.plot(x, col_sums, alpha=0.6, label='Phase synchronization')

		# # x 是索引数组，y 是编辑距离数组
		x1 = np.array(x)  # 示例索引数组
		y1 = np.array(col_sums)  # 示例数据，替换为你的实际数据

		# 设置参数：height=None 表示不限制高度，distance 表示峰之间的最小距离
		min_peaks, _ = find_peaks(y1, distance= 1 if int(mode_value*0.9) <= 1 else int(mode_value*0.9))
		# 标记朝下的峰
		ax2.scatter(x1[min_peaks], col_sums[min_peaks], color='red', zorder=2, label='Peaks',s=5)

		# outlier_indicesx = x1[outlier_mask_]
		# outlier_indicesy = col_sums[outlier_mask_]
		# ax2.scatter(outlier_indicesx, outlier_indicesy, color='green', zorder=5, label='Downward Peaks',s=6)

		peakx1 = x1[min_peaks]
		peaky1 = col_sums[min_peaks]

		if len(peakx1) > 10:
			y_array = np.array(peaky1)
			top_indices = np.argsort(y_array)[-5:]
			top_x = [peakx1[i] for i in top_indices]
			top_y = [peaky1[i] for i in top_indices]
			px_c = sum(top_y)/len(top_y)
		else:
			top_x = peakx1
			top_y = peaky1
			px_c = sum(top_y)/len(top_y)
		# for x,y in zip(top_x,top_y):
		# 	ax2.text(x, y, str(x)+'_'+str(y), fontsize=12, fontproperties=self.font_prop, verticalalignment='bottom',horizontalalignment='center', rotation=90)
		ax2.grid(True)
		# 使用嵌套循环计算差值
		differences = []
		for i in range(len(peakx1)):
			if i+2 >= len(peakx1):
				continue
			for j in range(i + 1, i + 2):
				seedlength = abs(peakx1[i] - peakx1[j])
				if seedlength >= self.minlength:
					differences.append(seedlength)  # 计算绝对差值
		if len(differences) >= 2:
			from statistics import mode,mean
			# 计算数组的众数
			mode_value = mode(differences)
			# 计算数组的均值
			mean_value = mean(differences)
			number = differences.count(mode_value)
			for i in range(1,int(mode_value*0.1)+1):
				number += differences.count(mode_value+i)
				number += differences.count(mode_value-i)
			number += 2
			# print(name,'众数:',mode_value,differences.count(mode_value)+2,number)
		else:
			mode_value = 0
			number = 0
		plt.xlabel(name, fontsize=12)  # 设置 x 轴标题为 'Index'，字体大小为 12
		plt.savefig(self.workpaths+'./hist/'+name+'.png', dpi = 600, bbox_inches = 'tight')# 存储图片
		plt.close()
		gc.collect()
		GC = (s1.count('G')+s1.count('C'))/len(s1)
		seq_mode = mode_value
		seq0 = s1[:seq_mode]
		return [length,GC,number,mode_value1,mean_value1,median_value1,variance_value1,max_count,px_c,seq_mode,block]+matrix_s+[seq0]

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

		if not os.path.exists(self.workpaths+self.pilpath):
			os.makedirs(self.workpaths+self.pilpath)
		if not os.path.exists(self.workpaths+self.histpath):
			os.makedirs(self.workpaths+self.histpath)

		# if not os.path.exists(self.workpaths+'./plotcsv/'):
		# 	os.makedirs(self.workpaths+'./plotcsv/')
		# if not os.path.exists(self.workpaths+'./clearcsv/'):
		# 	os.makedirs(self.workpaths+'./clearcsv/')
		f0 = open(self.workpaths+self.fastapath+self.outfile+'.repeat.index.txt','w')
		f0.write('\t'.join(['split_name','split_length','seq','length','number','index(Index starts from 0)'])+'\n')
		f0.close()
		f = open(self.workpaths+self.fastapath+self.outfile+'.split.fasta','w')
		f1 = open(self.workpaths+self.fastapath+self.outfile+'.split.gff','w')
		results = {}
		for key in genome.keys():
			s = str(genome[key].seq).upper()
			job = int(len(s)/self.windows)
			if job == 0:
				job = 1
			for i in range(job):
				if i == job - 1:
					s0 = s[self.windows*i:]
				else:
					s0 = s[self.windows*i:self.windows*(i+1)]
				s0 = self.merge_single_base_repeats(s0)
				s1,s2 = s0,s0
				if 'N' in s1:
					continue
				name = key+'_s'+str(i)
				f.write('>'+name+'\n'+s1+'\n')		
				f1.write('\t'.join([key,name,str(self.windows*i),str(self.windows*(i+1))])+'\n')	
				# if os.path.exists(self.pilpath+name+'.png') and os.path.exists(self.plotpath+name+'.png'):
				# 	continue
				dic = {}
				dic['gff'] = [key,'dotplot',name,str(self.windows*i+1),'.']
				dic['value'] = pool.apply_async(self.start, (s1,s2,name),  error_callback=bez.print_error)
				
				results[name] = dic
		pool.close() # Step IV : 准备结束
		pool.join() # Step IV : 完全结束
		f.close()
		f1.close()
		# print("")
		print("********************输出********************")

		names = list(results.keys())
		maker = ['length','GC','number','mode_value','mean_value','median_value','variance_value','max_count','AMPD_top','seed_l','block']+['M'+str(i) for i in range(1,201)]+['seed']
		f = open(self.workpath+self.outfile,'w')
		f.write('\t'.join(['chr','dotplot','name','start','end']+maker)+'\n')
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


