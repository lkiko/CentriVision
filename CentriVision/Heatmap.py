from tqdm import trange
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import Normalize, ListedColormap,LinearSegmentedColormap
from io import StringIO
import subprocess
import tempfile
import shutil
import os
import multiprocessing
from Bio import AlignIO
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import CentriVision.bez as bez
import sys
from multiprocessing import cpu_count#读取CPU核心数用于匹配线程数


class Heatmap():
	def __init__(self, options):
		self.workpath = os.getcwd()+'/'
		self.worktemp = 'hmap_temp/'
		self.temp_score = 'temp_score.txt'
		self.aligned_cpu = cpu_count()
		self.reverse_complement = "False"
		self.system = sys.platform
		self.gradient = 'rainbow'
		self.colors = "#7a4171,#aa4c8f,#cc7eb1,#674598,#5654a2,#4d5aaf,#8491c3,#007bbb,#0095d9,#59b9c6,#c1e4e9,#ebf6f7,#eaf4fc,#f8f4e6,#f2f2b0,#f8e58c,#ddbb99,#ebd842,#fef263,#ffd900,#f6ad49,#f8b500,#ec6d51,#eb6101,#e2041b"
		self.min_score = '0'
		self.split = '1100'
		self.segment_length = '0'
		self.background = 'white'
		# 初始化
		self.chip = {}
		self.chipc = ''
		self.trfgff = ''
		self.trfdic =''
		self.trfc = ''
		self.tegff = ''
		self.tedic = ''
		self.tec = ''
		self.genegff = ''
		self.genedic = ''
		self.genec = ''
		self.hpc = 'False'
		self.gffcdic = {}
		self.model = 'global'# "local"
		self.windows = 0
		self.step = 0

		bez_conf = bez.config()
		for k, v in bez_conf:# 
			setattr(self, str(k), v)
		for k, v in options:
			setattr(self, str(k), v)
			print(str(k), ' = ', v)
		self.aligned_cpu = int(self.aligned_cpu)
		self.min_score = float(self.min_score)
		self.segment_length = int(self.segment_length)
		self.split = int(self.split)
		self.windows = int(self.windows)
		self.step = int(self.step)

	def plot_similarity_distribution_subplot(self,score, ax,cmap):

		if self.background != 'white':
			# 修改子图背景色
			ax.set_facecolor('black')
			# 修改坐标轴颜色
			ax.spines['bottom'].set_color('white')
			ax.spines['top'].set_color('white') 
			ax.spines['right'].set_color('white')
			ax.spines['left'].set_color('white')
			# 修改刻度文字颜色
			ax.tick_params(axis='x', colors='white')
			ax.tick_params(axis='y', colors='white')
			# 修改坐标轴标签颜色
			ax.xaxis.label.set_color('white')
			ax.yaxis.label.set_color('white')
			axis_color = 'white'
		else:
			# 修改子图背景色
			ax.set_facecolor('white')
			# 修改坐标轴颜色
			ax.spines['bottom'].set_color('black')
			ax.spines['top'].set_color('black') 
			ax.spines['right'].set_color('black')
			ax.spines['left'].set_color('black')
			# 修改刻度文字颜色
			ax.tick_params(axis='x', colors='black')
			ax.tick_params(axis='y', colors='black')
			# 修改坐标轴标签颜色
			ax.xaxis.label.set_color('black')
			ax.yaxis.label.set_color('black')
			axis_color = 'black'

		# 绘制 kdeplot
		ax = sns.kdeplot(score, color="white", linewidth=1, bw_adjust=0.5, alpha=0, ax=ax)
		line = ax.lines[0]  # 获取曲线对象
		# 获取曲线的坐标
		x_curve = line.get_xdata()
		y_curve = line.get_ydata()

		# 创建 Normalize 对象，映射 score 到 [0, 1]
		# norm = Normalize(vmin=min(score), vmax=max(score))
		norm = Normalize(vmin=-100, vmax=100)
		# 计算 y 轴上限
		y_upper_limit = np.max(y_curve) * 1.3  # 将上限设置为曲线 y_curve 的 130%
		# 绘制渐变色条
		gradient = np.linspace(0, 1, 256).reshape(1, -1)
		# ax.imshow(gradient, aspect="auto", cmap=cmap, extent=(min(x_curve), max(x_curve), 0, y_upper_limit))  # 设置 extent 保持 x 轴和 y 轴范围一致
		ax.imshow(gradient, aspect="auto", cmap=cmap, extent=(-100, 100, -(y_upper_limit/10), y_upper_limit))  # 设置 extent 保持 x 轴和 y 轴范围一致

		# 取消上框线和右框线
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		# 保留底部轴线可见
		ax.spines['left'].set_visible(True)
		ax.spines['bottom'].set_visible(True)
		# ax.spines['bottom'].set_linewidth(2)
		# 设置 x 轴刻度
		ax.tick_params(axis='x', which='both', bottom=True, top=False, 
					   labelbottom=True, length=4, width=1, color=axis_color)
		ax.tick_params(axis='y', which='both', bottom=True, top=False, 
					   labelbottom=True, length=4, width=1, color=axis_color)
		ax.axis('on')
		ax.grid(False)

		# 绘制 kdeplot
		if self.background == 'white':
			ax = sns.kdeplot(score, color="white", linewidth=1, bw_adjust=0.5, ax=ax)
		else:
			ax = sns.kdeplot(score, color="black", linewidth=1, bw_adjust=0.5, ax=ax)
		line = ax.lines[0]  # 获取曲线对象
		line.set_alpha(0)  # 将曲线设置为透明

		# 获取曲线的坐标
		x_curve = list(line.get_xdata())
		y_curve = list(line.get_ydata())
		# print(x_curve,y_curve)
		x,y = [],[]
		for i in range(-100,int(x_curve[0])):
			x.append(i)
			y.append(0)
		x_curve = x + x_curve
		y_curve = y + y_curve
		x,y = [],[]
		for i in range(int(x_curve[-1]),101):
			x.append(i)
			y.append(0)
		x_curve = x_curve + x
		y_curve = y_curve + y

		# 设置图形轴标签
		ax.set_xlabel('identity')
		ax.set_ylabel('')
		
		if self.background == 'white':
			# 填充曲线上方的区域
			ax.fill_between(x_curve, y1=y_curve, y2=y_upper_limit, color="white", alpha=1)
			ax.axhline(y=0, ls='-', c='white')
		else:
			# 填充曲线上方的区域
			ax.fill_between(x_curve, y1=y_curve, y2=y_upper_limit, color="black", alpha=1)
			ax.axhline(y=0, ls='-', c='black')

	def fill_area(self,plt,position,r,color,text,parameter):
		# print(position)
		# r = r/2
		x0,y0 = position[0],position[1]
		# 定义四个坐标点
		x = [x0, x0-r, x0, x0+r, x0]
		y = [y0+r, y0, y0-r, y0, y0+r]
		# # 绘制四个坐标点之间的线段
		# plt.plot(x, y, color='black')
		# 填充四边形区域
		# plt.axis('equal')  # 保持 x 和 y 轴的相同比例
		plt.fill(x, y, color=color, edgecolor='none')# , alpha=0.5
		if parameter:
			plt.text(x0, y0, text, color='black', ha='center', va='center', fontsize=12)

	def calculate_similarity(self,seq1, seq2):# 同位置比对上则得分
		s1,s2 = "",""
		for i in range(len(seq1)):
			if seq1[i] == '-' and seq2[i] == '-':
				continue
			else:
				s1 += seq1[i]
				s2 += seq2[i]
		return 100*sum(a == b for a, b in zip(s1, s2))/len(s1)

	def align(self, input_file, output_dir):
		output_file = os.path.join(output_dir, "align_output.fasta")
		if self.align_software == 'mafft':
			command = [self.mafft_path, '--auto', input_file]
			result = subprocess.run(command, capture_output=True, text=True)
			align = AlignIO.read(StringIO(result.stdout), "fasta")
			AlignIO.write(align, output_file, "fasta")
			return output_file
		elif self.align_software == 'muscle':
			command = [self.muscle_path, '-in', input_file, '-out', output_file]
			result = subprocess.run(command, capture_output=True, text=True)
			return output_file
		else:
			print('多序列比对软件支持muscle/mafft 请重新选择！')
			exit()

	def align1(self, input_file, output_dir):
		output_file = os.path.join(output_dir, "align_output.fasta")
		if self.align_software == 'mafft':
			try:
				command = [self.mpirun_path,'-np',str(self.aligned_cpu),self.mafft_path, '--auto', input_file]
				result = subprocess.run(command, capture_output=True, text=True)
			except:
				command = [self.mafft_path, '--auto', input_file]
				result = subprocess.run(command, capture_output=True, text=True)
			align = AlignIO.read(StringIO(result.stdout), "fasta")
			AlignIO.write(align, output_file, "fasta")
			return output_file
		elif self.align_software == 'muscle':
			try:
				command = [self.mpirun_path,'-np',str(self.aligned_cpu),self.muscle_path, '-in', input_file, '-out', output_file]
				result = subprocess.run(command, capture_output=True, text=True)
			except:
				command = [self.muscle_path, '-in', input_file, '-out', output_file]
				result = subprocess.run(command, capture_output=True, text=True)
			return output_file
		else:
			print('多序列比对软件支持muscle/mafft 请重新选择！')
			exit()

	def get_complementary_dna(self,dna):
		complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
		return ''.join([complement[base] for base in dna])

	def pairwise_alignment_hamming_distance(self,num,job,seq01, seq02, file,index1,index2,path):
		score1 = self.calculate_similarity(seq01, seq02)
		if self.reverse_complement == 'True':# 反向互补
			seq03 = self.get_complementary_dna(str(seq02)[::-1])
		else:# 反向
			seq03 = str(seq02)[::-1]
		score2 = -self.calculate_similarity(seq01, seq03)	
		if score1 >= abs(score2):
			score = score1
		else:
			score = score2
		out = ('1', '2', score,(int(index1),int(index2)))
		p = int(((num+1)/job)*30)
		print("\r["+"*"*p+" "*(30-p)+"] "+str(num)+"/"+str(job)+' add alignment',end="")
		return out

	def pairwise_alignment(self,num,job,seq01, seq02, file,index1,index2,path,dicscore):
		if str(int(index1))+'_'+str(int(index2)) in dicscore.keys():
			lt = dicscore[str(int(index1))+'_'+str(int(index2))]
			out = (lt[0], lt[1], float(lt[2]),(int(lt[3]),int(lt[4])))
			p = int(((num+1)/job)*30)
			print("\r["+"*"*p+" "*(30-p)+"] "+str(num)+"/"+str(job)+' add alignment',end="")
			return out

		os.chdir(path)
		os.makedirs(file, exist_ok=True)
		os.chdir(path+file)

		# # 代码加速
		# seq1 = SeqIO.SeqRecord(Seq(seq01), id="Seq1",name="Seq1",description='')
		# seq2 = SeqIO.SeqRecord(Seq(seq02), id="Seq2",name="Seq2",description='')
		# # 将序列写入临时文件
		# 
		# SeqIO.write([seq1, seq2], input_file, "fasta")

		input_file = os.path.join(path+file, "input.fasta")
		f = open(input_file,'w')
		f.write('\n'.join(['>Seq1',seq01,'>Seq2',seq02]))
		f.close()

		align_output = self.align(input_file, path+file)
		# 读取 MAFFT 输出的比对结果
		alignment = AlignIO.read(align_output, "fasta")
		# 获取比对后的序列
		aligned_seq1, aligned_seq2 = str(alignment[0].seq), str(alignment[1].seq)
		# 计算得分
		score1 = self.calculate_similarity(aligned_seq1, aligned_seq2)

		# if self.reverse_complement == 'True':# 反向互补
		# 	seq3 = SeqIO.SeqRecord(Seq(self.get_complementary_dna(str(seq02)[::-1])), id="Seq3",name="Seq3",description='')
		# else:# 反向
		# 	seq3 = SeqIO.SeqRecord(Seq(str(seq02)[::-1]), id="Seq3",name="Seq3",description='')
		# # 将序列写入临时文件
		# 
		# SeqIO.write([seq1, seq3], input_file, "fasta")

		if self.reverse_complement == 'True':# 反向互补
			seq03 = self.get_complementary_dna(str(seq02)[::-1])
		else:# 反向
			seq03 = str(seq02)[::-1]
		# 将序列写入临时文件
		input_file = os.path.join(path+file, "input2.fasta")
		f = open(input_file,'w')
		f.write('\n'.join(['>Seq1',seq01,'>Seq3',seq03]))
		f.close()

		align_output = self.align(input_file, path+file)
		# 读取 MAFFT 输出的比对结果
		alignment = AlignIO.read(align_output, "fasta")
		# 获取比对后的序列
		aligned_seq1, aligned_seq2 = str(alignment[0].seq), str(alignment[1].seq)
		# 计算得分
		score2 = -self.calculate_similarity(aligned_seq1, aligned_seq2)	
		if score1 >= abs(score2):
			score = score1
		else:
			score = score2
		os.chdir(path)
		shutil.rmtree(file)
		out = ('1', '2', score,(int(index1),int(index2)))

		f0 = open(path+self.temp_score,'a')
		# print(str(int(index1))+'_'+str(int(index2)))
		f0.write('\t'.join([str(int(index1))+'_'+str(int(index2)),'1','2',str(score),str(int(index1)),str(int(index2))])+'\n')
		f0.close()

		p = int(((num+1)/job)*30)
		print("\r["+"*"*p+" "*(30-p)+"] "+str(num)+"/"+str(job)+' add alignment',end="")
		return out

	def pairwise_alignment_c(self,num,job,seq01, seq02, file,index1,index2,path,dicscore):

		import ctypes
		class AlignResult(ctypes.Structure):
			_fields_ = [("aligned_seq1", ctypes.c_char_p),
						("aligned_seq2", ctypes.c_char_p),
						("score", ctypes.c_int)]

		def load_library():
			lib_name = "myalign"
			base_dir = os.path.join(os.path.dirname(__file__), 'Lib')

			if sys.platform.startswith("win"):  # Windows
				lib_path = os.path.join(base_dir, lib_name + ".dll")
			else:  # Linux / Mac
				lib_path = os.path.join(base_dir, lib_name + ".so")

			if not os.path.exists(lib_path):
				raise FileNotFoundError(f"Library not found: {lib_path}")

			lib = ctypes.CDLL(lib_path)

			# 设置函数接口
			lib.nw_align_wrapper.restype = ctypes.POINTER(AlignResult)
			lib.nw_align_wrapper.argtypes = [ctypes.c_char_p, ctypes.c_char_p,
											 ctypes.c_int, ctypes.c_int, ctypes.c_int]
			lib.free_align_result.argtypes = [ctypes.POINTER(AlignResult)]
			return lib
		lib = load_library()
		def align(seq1: str, seq2: str, match: int = 1, mismatch: int = -1, gap: int = -2):
			"""
			对齐两个序列，返回 aligned_seq1, aligned_seq2, score(0~100)
			"""
			seq1_bytes = seq1.encode('utf-8')
			seq2_bytes = seq2.encode('utf-8')
			# 调用 C 函数
			res_ptr = lib.nw_align_wrapper(seq1_bytes, seq2_bytes, match, mismatch, gap)
			res = res_ptr.contents
			aligned_seq1 = res.aligned_seq1.decode()
			aligned_seq2 = res.aligned_seq2.decode()
			score = res.score
			# 映射到 0~100
			max_score = min(len(seq1), len(seq2)) * match
			normalized_score = (score / max_score) * 100
			# 释放内存
			lib.free_align_result(res_ptr)
			return aligned_seq1, aligned_seq2, normalized_score

		a1, a2, score1 = align(seq01, seq02)
		# print(score1,'输出', flush=True)
		# 反向互补比对
		if self.reverse_complement == 'True':
			seq03 = self.get_complementary_dna(seq02[::-1])
		else:
			seq03 = seq02[::-1]
		a1, a2, score2 = align(seq01, seq03)
		score2 = -score2

		if score1 >= abs(score2):
			score = score1
		else:
			score = score2

		out = ('1', '2', score,(int(index1),int(index2)))
		f0 = open(path+self.temp_score,'a')
		f0.write('\t'.join([str(int(index1))+'_'+str(int(index2)),'1','2',str(score),str(int(index1)),str(int(index2))])+'\n')
		f0.close()

		p = int(((num+1)/job)*30)
		print("\r["+"*"*p+" "*(30-p)+"] "+str(num)+"/"+str(job)+' add alignment',end="")
		return out

	def compare_sequence_segments_step_mpi(self, sequence, segment_length):

		dicscore = {}
		# 判断文件是否存在
		if os.path.exists(self.temp_score):
			# out = ('1', '2', score,(int(index1),int(index2)))
			# score_ed = 
			# 读取文件并处理
			# 创建空字典
			with open(self.temp_score, 'r') as file:
				for line in file:
					# 移除行末的换行符并按空格分割
					parts = line.strip().split()
					# 第一列作为键，其余列合并为列表作为值
					key = parts[0]
					values = parts[1:]
					# 将键值对添加到字典中
					dicscore[key] = values
		else:
			f = open(self.temp_score,'w')
			f.close()


		from mpi4py import MPI
		from itertools import combinations
		comm = MPI.COMM_WORLD
		rank = comm.Get_rank()
		size = comm.Get_size()
		num_segments = int(len(sequence) / segment_length)
		all_pairs = list(combinations(range(num_segments), 2))
		# Distribute work among processes
		my_pairs = np.array_split(all_pairs, size)[rank]
		similarities = []
		for pair in my_pairs:
			i, j = pair
			segment1 = sequence[i * segment_length:(i + 1) * segment_length]
			segment2 = sequence[j * segment_length:(j + 1) * segment_length]
			name = f"{i * segment_length - (segment_length // 2)}_{j * segment_length - (segment_length // 2)}"
			result = self.pairwise_alignment(i, j, segment1, segment2, name, i * segment_length - (segment_length // 2), j * segment_length - (segment_length // 2),self.workpath+self.worktemp,dicscore)
			similarities.append(result)
		# Gather results from all processes
		all_similarities = comm.gather(similarities, root=0)
		if rank == 0:
			flattened_similarities = [item for sublist in all_similarities for item in sublist]
			similarity_scores = [similarity[2] for similarity in flattened_similarities]
			min_score = min(similarity_scores)
			max_score = max(similarity_scores)
		else:
			similarity_scores = None
			min_score = None
			max_score = None
		similarity_scores = comm.bcast(similarity_scores, root=0)
		min_score = comm.bcast(min_score, root=0)
		max_score = comm.bcast(max_score, root=0)

		# 检查文件是否存在并删除
		if os.path.exists(self.temp_score):
			os.remove(self.temp_score)

		return all_similarities, min_score, max_score, similarity_scores

	def compare_sequence_segments_step(self, sequence, segment_length):
		dicscore = {}
		# 判断文件是否存在
		if os.path.exists(self.temp_score):
			# 读取文件并处理
			# 创建空字典
			with open(self.temp_score, 'r') as file:
				for line in file:
					# 移除行末的换行符并按空格分割
					parts = line.strip().split()
					# 第一列作为键，其余列合并为列表作为值
					key = parts[0]
					values = parts[1:]
					# 将键值对添加到字典中
					dicscore[key] = values
		else:
			f = open(self.temp_score,'w')
			f.close()

		num_segments = int(len(sequence) / segment_length)
		similarities = []
		job = sum(range(num_segments))
		num = 1
		with multiprocessing.Pool(processes=self.aligned_cpu) as pool:
			for i in range(num_segments):
				segment1 = sequence[i * segment_length:(i + 1) * segment_length]
				for j in range(i + 1, num_segments):
					segment2 = sequence[j * segment_length:(j + 1) * segment_length]
					name = f"{i * segment_length - (segment_length // 2)}_{j * segment_length - (segment_length // 2)}"
					similarities.append(pool.apply_async(self.pairwise_alignment, 
														 (num, job, segment1, segment2, name, 
														  i * segment_length - (segment_length // 2), 
														  j * segment_length - (segment_length // 2),self.workpath+self.worktemp,dicscore)))
					num += 1
			pool.close()
			pool.join()

		print('\n')
		similarities = [similarities[i].get() for i in trange(len(similarities))]
		similarity_scores = [similarities[i][2] for i in trange(len(similarities))]

		# 检查文件是否存在并删除
		if os.path.exists(self.temp_score):
			os.remove(self.temp_score)

		return similarities, min(similarity_scores), max(similarity_scores), similarity_scores

	def compare_sequence_segments_step_my(self, sequence, segment_length):

		dicscore = {}
		# 判断文件是否存在
		if os.path.exists(self.temp_score):
			# out = ('1', '2', score,(int(index1),int(index2)))
			# score_ed = 
			# 读取文件并处理
			# 创建空字典
			with open(self.temp_score, 'r') as file:
				for line in file:
					# 移除行末的换行符并按空格分割
					parts = line.strip().split()
					# 第一列作为键，其余列合并为列表作为值
					key = parts[0]
					values = parts[1:]
					# 将键值对添加到字典中
					dicscore[key] = values
		else:
			f = open(self.temp_score,'w')
			f.close()

		num_segments = int(len(sequence) / segment_length)
		similarities = []
		job = sum(range(num_segments))
		num = 1
		with multiprocessing.Pool(processes=self.aligned_cpu) as pool:
			for i in range(num_segments):
				segment1 = sequence[i * segment_length:(i + 1) * segment_length]
				for j in range(i + 1, num_segments):
					segment2 = sequence[j * segment_length:(j + 1) * segment_length]
					name = f"{i * segment_length - (segment_length // 2)}_{j * segment_length - (segment_length // 2)}"
					similarities.append(pool.apply_async(self.pairwise_alignment_c, 
														 (num, job, segment1, segment2, name, 
														  i * segment_length - (segment_length // 2), 
														  j * segment_length - (segment_length // 2),self.workpath+self.worktemp,dicscore)))
					num += 1
			pool.close()
			pool.join()

		print('\n')
		similarities = [similarities[i].get() for i in trange(len(similarities))]
		similarity_scores = [similarities[i][2] for i in trange(len(similarities))]

		# 检查文件是否存在并删除
		if os.path.exists(self.temp_score):
			os.remove(self.temp_score)

		return similarities, min(similarity_scores), max(similarity_scores), similarity_scores

	def compare_sequence_segments_step_hamming_distance(self, sequence, segment_length):
		num_segments = int(len(sequence) / segment_length)
		similarities = []
		job = sum(range(num_segments))
		num = 1
		with multiprocessing.Pool(processes=self.aligned_cpu) as pool:
			for i in range(num_segments):
				segment1 = sequence[i * segment_length:(i + 1) * segment_length]
				for j in range(i + 1, num_segments):
					segment2 = sequence[j * segment_length:(j + 1) * segment_length]
					name = f"{i * segment_length - (segment_length // 2)}_{j * segment_length - (segment_length // 2)}"
					similarities.append(pool.apply_async(self.pairwise_alignment_hamming_distance, 
														 (num, job, segment1, segment2, name, 
														  i * segment_length - (segment_length // 2), 
														  j * segment_length - (segment_length // 2),self.workpath+self.worktemp)))
					num += 1

			pool.close()
			pool.join()

		print('\n')
		similarities = [similarities[i].get() for i in trange(len(similarities))]
		similarity_scores = [similarities[i][2] for i in trange(len(similarities))]
		return similarities, min(similarity_scores), max(similarity_scores), similarity_scores

	def compare_sequence_segments_step_linux(self, sequence, segment_length):
		# print(segment_length)
		num_segments = int(len(sequence) / segment_length)
		similarities = []
		job = sum(range(num_segments))
		num = 1
		f = open('align.fasta','w')
		for i in range(num_segments):
			segment1 = sequence[i * segment_length:(i + 1) * segment_length]
			name = str(i)+'_0'
			f.write(">"+name+"\n"+segment1+"\n")
			name = str(i)+'_1'
			if self.reverse_complement == 'True':# 反向互补
				f.write(">"+name+"\n"+self.get_complementary_dna(str(segment1)[::-1])+"\n")
			else:# 反向
				f.write(">"+name+"\n"+str(segment1)[::-1]+"\n")
		f.close()
		align_output = self.align1('align.fasta', './')
		genome = {}
		for seq_record in SeqIO.parse('align_output.fasta', "fasta"):# biopython循环读取,去除冗余
			if seq_record.id not in genome.keys():
				genome[seq_record.id] = str(seq_record.seq)
			else:
				pass
			# if seq_record.id == '1':
			# 	print(str(seq_record.seq))
		for i in range(num_segments):
			segment1 = str(genome[str(i)+'_0']).upper()
			for j in range(i + 1, num_segments):
				name = f"{i * segment_length - (segment_length // 2)}_{j * segment_length - (segment_length // 2)}"
				segment2 = str(genome[str(j)+'_0']).upper()# 正向比对
				# 计算得分
				score1 = self.calculate_similarity(segment1, segment2)

				segment2 = str(genome[str(j)+'_1']).upper()# 反向比对
				# 计算得分
				score2 = -self.calculate_similarity(segment1, segment2)
				if score1 >= abs(score2):
					score = score1
				else:
					score = score2

				out = ('1', '2', score,(int(i * segment_length - (segment_length // 2)),int(j * segment_length - (segment_length // 2))))
				similarities.append(out)
		similarities = [similarities[i] for i in trange(len(similarities))]
		similarity_scores = [similarities[i][2] for i in trange(len(similarities))]
		return similarities, min(similarity_scores), max(similarity_scores), similarity_scores

	def compare_sequence_segments(self,sequence, segment_length):
		similaritys = []
		num_segments = len(sequence) - segment_length + 1
		similarities = []
		for i in trange(num_segments):
			segment1 = sequence[i:i+segment_length]
			for j in range(i+1, num_segments):
				segment2 = sequence[j:j+segment_length]
				similarity = self.calculate_similarity(segment1, segment2)
				# relative_position = (i, j)
				relative_position = (i*segment_length-(segment_length/2), j*segment_length-(segment_length/2))
				similarities.append((segment1, segment2, similarity, relative_position))
				similaritys.append(similarity)
		return similarities,min(similaritys),max(similaritys),similaritys

	def value_to_color(self,value, cmap, vmin, vmax):
		norm = Normalize(vmin=vmin, vmax=vmax)
		return cmap(norm(value))


	def plot(self,ax,x,y,w,h,c,a,p,n,key):
		rect = plt.Rectangle((x,y),w,h,facecolor=c,alpha=a, label=key, edgecolor='none')# x,y,宽，高
		ax.add_patch(rect)


	def line_plot(self,data,length,yindex,h,length0,segment_length):
		x,y = [],[]
		if self.windows == 0:
			windows = int(length/100)
			step = int(windows/2)
		else:
			windows = self.windows
			step = self.step
		start,end = -step,windows-step
		while end < length:
			start,end = start+step,end+step
			index = (start+end)/2
			data0 = 100*sum(data[start:end])/(end-start)
			x.append(index)
			y.append(data0)
		# y_max = max(y)
		y_max = 100

		x = [i-(segment_length/2) for i in x]
		y = list(map(lambda a: (a/y_max)*h+yindex, y))
		x = list(map(lambda a: (a/length)*length0, x))
		return x,y

	def covarageGC(self,chrseq,length0,chro,yindex,h,segment_length):
		length = len(chrseq)
		array_d = np.zeros(length, dtype = int)
		for i in range(length):
			if chrseq[i] == 'C' or chrseq[i] == 'G':
				array_d[i] = 1
		x1,y1 = self.line_plot(array_d,length,yindex,h,length0,segment_length)
		return x1,y1,array_d

	def covarage(self,suf_dic,length,length0,chro,yindex,h,segment_length):
		array_d = np.zeros(length, dtype = int)
		for start,end in suf_dic:
			start = int(start - 1)
			end = int(end)
			# print(start,end,length)
			array_d[start:end] = array_d[start:end]+1
		array_d1 = np.int64(array_d>0)
		x1,y1 = self.line_plot(array_d1,length,yindex,h,length0,segment_length)
		return x1,y1,array_d1


	def read_overlap(self,ax,suf_dic,length,chro,chrseq,segment_length):
		# h0 = -int(length/((2**0.5)*10))
		# h0 = -int(len(chrseq)/((2**0.5)*10))
		h0 = int(len(chrseq)/(2*10.5))
		h1 = int(len(chrseq)/(2*10))
		print('centromere')
		# ax.text(-0.1, h, 'stq'+chro[3:], fontdict=font0,color = 'black')
		self.plot(ax,-segment_length/2,-h1,length,h0,'#d3cfd9',0.5,False,chro,'')
		# ax.plot([-segment_length/2, length-segment_length/2], [0,0], color='r', linestyle='--')
		# ax.plot([-segment_length/2, length-segment_length/2], [h0,h0], color='r', linestyle='--')
		
		if chro in self.chip.keys():
			for centromere in self.chip[chro]:
				start,end = min(centromere[0],centromere[1]),max(centromere[0],centromere[1])
				start0 = (start/len(chrseq))*length-segment_length/2
				w = ((abs(end-start)/len(chrseq))*length)
				self.plot(ax,start0,-h1,w,h0,self.chipc,0.5,False,chro,'chip')

		x_d1,y_d1,array_d1 = self.covarageGC(chrseq,length,chro,-h1,h0,segment_length)
		cells = ax.plot(x_d1,y_d1, color = self.gcc, linestyle = '-', linewidth=.5, label='GC', alpha=0.5)# , linewidth=.5

		for key in ['TRF','TE','GENE']:
			print(key)
			if key in suf_dic.keys():
				x_d1,y_d1,array_d1 = self.covarage(suf_dic[key],len(chrseq),length,chro,-h1,h0,segment_length)
				cells = ax.plot(x_d1,y_d1,  color = self.gffcdic[key], linestyle = '-', linewidth=.8, label=key, alpha=0.7)# , linewidth=2

	def run_c(self,sequence,segment_length,name,length):
		fig, ax = plt.subplots(figsize=(10, 7.5))

		if self.background != 'white':
			# 修改图表背景色
			fig.patch.set_facecolor('black')
			# 修改子图背景色
			ax.set_facecolor('black')
			# 修改坐标轴颜色
			ax.spines['bottom'].set_color('white')
			ax.spines['top'].set_color('white') 
			ax.spines['right'].set_color('white')
			ax.spines['left'].set_color('white')

			# 修改刻度文字颜色
			ax.tick_params(axis='x', colors='white')
			ax.tick_params(axis='y', colors='white')

			# 修改坐标轴标签颜色
			ax.xaxis.label.set_color('white')
			ax.yaxis.label.set_color('white')
			axis_color = 'white'
		else:
			# 修改图表背景色
			fig.patch.set_facecolor('white')
			# 修改子图背景色
			ax.set_facecolor('white')
			# 修改坐标轴颜色
			ax.spines['bottom'].set_color('black')
			ax.spines['top'].set_color('black') 
			ax.spines['right'].set_color('black')
			ax.spines['left'].set_color('black')

			# 修改刻度文字颜色
			ax.tick_params(axis='x', colors='black')
			ax.tick_params(axis='y', colors='black')

			# 修改坐标轴标签颜色
			ax.xaxis.label.set_color('black')
			ax.yaxis.label.set_color('black')
			axis_color = 'black'

		import time
		start = time.perf_counter()
		if self.model == 'global':
			similarities,vmin,vmax,score = self.compare_sequence_segments_step_linux(sequence, segment_length)
		else:
			if self.align_software == 'hamming':
				similarities,vmin,vmax,score = self.compare_sequence_segments_step_hamming_distance(sequence, segment_length)
			else:
				if self.hpc != 'True':
					if self.align_software == 'ialign':
						similarities,vmin,vmax,score = self.compare_sequence_segments_step_my(sequence, segment_length)
					else:
						similarities,vmin,vmax,score = self.compare_sequence_segments_step(sequence, segment_length)
				else:
					similarities,vmin,vmax,score = self.compare_sequence_segments_step_mpi(sequence, segment_length)
		# 获取结束时间
		end = time.perf_counter()
		# 计算运行时间
		runTime = end - start
		runTime_ms = runTime * 1000
		# 输出运行时间
		print("多序列比对运行时间：", runTime, "秒")
		print("多序列比对运行时间：", runTime_ms, "毫秒")

		if self.color_mode == "Discrete":
			# 定义颜色列表，你可以根据需要调整颜色顺序和数量
			# 创建 ListedColormap 对象
			cmap = ListedColormap(self.colors.split(','))
		elif self.color_mode == "Gradient":
			# 创建渐变颜色映射
			# colors = sns.color_palette("rainbow", as_cmap=True)
			colors = sns.color_palette(self.gradient, as_cmap=True)
			cmap = LinearSegmentedColormap.from_list("Custom", colors(np.linspace(0, 1, 256)))
		else:
			print('当前软件仅支持Gradient和Discrete两种颜色模式！请重新选择！')
			exit()

		# r = 0.5
		r = segment_length/2
		for i in trange(len(similarities)):
			segment1, segment2, similarity, relative_position = similarities[i][0], similarities[i][1], similarities[i][2], similarities[i][3]
			if abs(similarity) <= self.min_score:
				continue

			x0,y0 = (relative_position[1] + relative_position[0])/2,(relative_position[1] - relative_position[0])/2
			# self.fill_area(plt,(x0,y0),r*0.98,self.value_to_color(similarity,cmap,vmin,vmax),str(similarity),False)
			self.fill_area(plt,(x0,y0),r*0.98,self.value_to_color(similarity,cmap,-100,100),str(similarity),False)

		# 取消上框线和右框线
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		# 关闭 y 轴坐标刻度
		ax.set_yticks([])
		ax.spines['left'].set_visible(False)
		ax.grid(False)  # 确保网格线关闭

		# 保留底部轴线可见
		ax.spines['bottom'].set_visible(True)
		ax.spines['bottom'].set_linewidth(2)

		# 设置 x 轴刻度
		ax.tick_params(axis='x', which='both', bottom=True, top=False, 
					   labelbottom=True, length=4, width=1, color=axis_color)

		# 设置坐标轴标签
		plt.xlabel(name+' (bp)')
		# plt.ylabel('Y')
		if self.annotation == 'True':
			plt.ylim(ymin=-int(len(sequence)/(2*10)),ymax=len(sequence)*0.75)
		else:
			plt.ylim(ymin=-2*r,ymax=len(sequence)*0.75)
		plt.xlim(xmin=-2*r,xmax=len(sequence))

		# 提前读取文件会导致降速
		if self.annotation == 'True':
			if self.chip_file != 'None':
				self.chip = self.read_chip("../"+self.chip_file.split(':')[0])
				self.chipc = self.chip_file.split(':')[1]
			if self.trf_gff != 'None':
				self.trfgff = pd.read_csv("../"+self.trf_gff.split(':')[0],header = None, sep='\t', comment='#').sort_values(by=[0,3],ascending= [True,True])
				self.trfdic = self.trfgff.groupby(0).groups# 按照第几列分组
				self.trfc = self.trf_gff.split(':')[1]
				self.trfgff[7] = abs(self.trfgff[4] - self.trfgff[3])
				self.trfgff[9] = self.trfgff.apply(lambda row: '('+str(row[3]) + ',' + str(row[4])+')', axis=1)
			if self.te_gff != 'None':
				self.tegff = pd.read_csv("../"+self.te_gff.split(':')[0],header = None, sep='\t', comment='#').sort_values(by=[0,3],ascending= [True,True])
				self.tedic = self.tegff.groupby(0).groups# 按照第几列分组
				self.tec = self.te_gff.split(':')[1]
				self.tegff[7] = abs(self.tegff[4] - self.tegff[3])
				self.tegff[9] = self.tegff.apply(lambda row: '('+str(row[3]) + ',' + str(row[4])+')', axis=1)
			if self.gene_gff != 'None':
				self.genegff = pd.read_csv("../"+self.gene_gff.split(':')[0],header = None, sep='\t', comment='#').sort_values(by=[0,3],ascending= [True,True])
				self.genedic = self.genegff.groupby(0).groups# 按照第几列分组
				self.genec = self.gene_gff.split(':')[1]
				self.genegff[7] = abs(self.genegff[3] - self.genegff[2])
				self.genegff[9] = self.genegff.apply(lambda row: '('+str(row[2]) + ',' + str(row[3])+')', axis=1)

		if self.annotation == 'True':
			suf_dic = {}
			if self.trf_gff != 'None':
				self.gffcdic['TRF'] = self.trfc
				trflocal = self.trfgff.loc[self.trfdic[name]].sort_values(by=[0,3],ascending= [True,True])
				trflocal.reset_index(drop=True,inplace=True)
				# print(trflocal)
				intervals = trflocal[9]
				intervals = list(map(eval, intervals))
				suf_dic['TRF'] = intervals
				# print(intervals)
				# exit()

			if self.te_gff != 'None':
				self.gffcdic['TE'] = self.tec
				telocal = self.tegff.loc[self.tedic[name]].sort_values(by=[0,3],ascending= [True,True])
				telocal.reset_index(drop=True,inplace=True)
				# print(telocal)
				intervals = telocal[9]
				intervals = list(map(eval, intervals))
				suf_dic['TE'] = intervals

			if self.gene_gff != 'None':
				self.gffcdic['GENE'] = self.genec
				genelocal = self.genegff.loc[self.genedic[name]].sort_values(by=[0,2],ascending= [True,True])
				genelocal.reset_index(drop=True,inplace=True)
				# print(telocal)
				intervals = genelocal[9]
				intervals = list(map(eval, intervals))
				suf_dic['GENE'] = intervals

			num_segments = int(len(sequence)/segment_length)-1
			length = num_segments*segment_length
			self.read_overlap(ax,suf_dic,length,name,sequence,segment_length)

			if self.background != 'white':
				# 显示图例并设置图例背景色
				legend = ax.legend(facecolor='black')
				# 设置图例文字颜色
				plt.setp(legend.get_texts(), color='white')
			else:
				# 显示图例
				plt.legend()

		# 初始化
		self.chip = {}
		self.chipc = ''
		self.trfgff = ''
		self.trfdic =''
		self.trfc = ''
		self.tegff = ''
		self.tedic = ''
		self.tec = ''
		self.gffcdic = {}

		#内嵌图
		a = plt.axes([0.2, 0.6, .2, .2], facecolor="white")
		self.plot_similarity_distribution_subplot(score, a,cmap)
		# plt.grid(False)
		# 显示图形
		plt.savefig(self.workpath+self.out_path+'/'+name+'-hmap.png', dpi=1000, bbox_inches='tight')
		print(name,'draw over!')
		# ashell = os.system('python /mycode/lan_email.py lanmeifang@ibcas.ac.cn p '+ self.workpath+self.out_path+'/'+name+'-hmap.png')

	def read_chip(self,file):
		zsl = {}
		for line in open(file,'r'):
			lt = line.strip('\n').split('\t')
			# lt[0] = 'chr'+lt[0][3:]
			if lt[0] == 'Chr_ID':
				continue
			if lt[0] not in zsl.keys():
				zsl[lt[0]] = []
				zsl[lt[0]].append([int(lt[1]),int(lt[2]),abs(int(lt[1])-int(lt[2]))+1])
			else:
				zsl[lt[0]].append([int(lt[1]),int(lt[2]),abs(int(lt[1])-int(lt[2]))+1])
		return zsl

	def run(self):
		if not os.path.exists(self.workpath+self.out_path):
			os.makedirs(self.workpath+self.out_path)

		if not os.path.exists(self.workpath+self.worktemp):
			os.makedirs(self.workpath+self.worktemp)

		genome = SeqIO.to_dict(SeqIO.parse(self.centromere_file, "fasta"))# 提取之后直接返回字典
		os.chdir(self.worktemp)
		
		if self.model == 'global':
			lens = []
			for key in genome.keys():
				l = len(genome[key].seq)
				lens.append(l)

			if self.segment_length == 0:
				segment_length = max(100, min(5000, int(len(sequence) / self.split)))
			else:
				segment_length = max(100, min(5000, self.segment_length))

			for key in genome.keys():
				name = key
				if os.path.exists(self.workpath+self.out_path+'/'+name+'-hmap.png'):
					continue
				sequence = str(genome[key].seq).upper()
				print(min(lens),max(lens),'segment_length:',segment_length)
				print(name,segment_length)
				self.run_c(sequence,segment_length,name,len(sequence))
		else:
			lens = []
			for key in genome.keys():
				l = len(genome[key].seq)
				lens.append(l)
			for key in genome.keys():
				name = key
				if os.path.exists(self.workpath+self.out_path+'/'+name+'-hmap.png'):
					continue
				sequence = str(genome[key].seq).upper()
				# segment_length = int(len(sequence)/self.split)

				if self.segment_length == 0:
					segment_length = max(100, min(5000, int(len(sequence) / self.split)))
				else:
					segment_length = max(100, min(5000, self.segment_length))

				print(min(lens),max(lens),'segment_length:',segment_length)
				print(name,segment_length)
				self.run_c(sequence,segment_length,name,len(sequence))
		os.chdir(self.workpath)
		shutil.rmtree(self.worktemp)

