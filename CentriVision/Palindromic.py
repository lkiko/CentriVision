# python
# -*- encoding: utf-8 -*-
'''
@File		:Palindromic.py
@Time		:2021/04/27 19:54:08
@Author		:charles kiko
@Version		:1.0
@Contact		:charles_kiko@163.com
@Desc		:python Palindromic.py oldgff oldpep oldcds oldgene genome xxx
@annotation		:
'''

import sys
import os
import ctypes
import gc# 内存管理模块
from tqdm import trange
from tqdm import tqdm
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np
import pandas as pd
import CentriVision.bez as bez
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
import multiprocessing # Step I : 导入模块
from multiprocessing import cpu_count#读取CPU核心数用于匹配线程数

def worker_chromosome(record, length, reach, code, queue, libpath=None):
	try:
		chrom_name = record.id
		seq = str(record.seq).upper()
		p = Palindromic([])
		if code == 'c':
			palindromes = p.run_palindrome_c_with_queue_chunk(seq, length, reach, queue, libpath)
		else:
			palindromes = p.find_palindromesL(chrom_name, seq, queue)
		return chrom_name, palindromes
	except Exception as e:
		import traceback
		tb = traceback.format_exc()
		if queue is not None:
			queue.put(("error", tb))
		raise




class Palindromic():
	def __init__(self, options):
		self.length = 10
		self.reach = 2000
		self.step = 5000
		self.windows = 10000
		self.coln = 3
		self.width = 15
		self.height = 10
		self.code = 'c'
		self.dpi = 1000
		self.cpu = cpu_count() - 1
		for k, v in options:
			setattr(self, str(k), v)
			print(k, ' = ', v)

		# self.maxp = int(self.maxp)
		# self.minp = int(self.minp)
		self.coln = int(self.coln)
		self.width = int(self.width)
		self.height = int(self.height)
		self.length = int(self.length)
		self.reach = int(self.reach)
		self.step = int(self.step)
		self.dpi = int(self.dpi)
		self.windows = int(self.windows)
		self.cpu = int(self.cpu)

	# def write_gff_file(self, chromosome_name, palindromes):
	# 	with open(self.gff_file, 'a') as f:
	# 		for start, end, seq0, length,types in palindromes:
	# 			f.write(f"{chromosome_name}\tCentriVision\t{types}\t{start+1}\t{end}\t.\t.\t.\tLength={length};seq={seq0}\n")

	def write_gff_file(self, chromosome_name, palindromes):
		fh = open(self.gff_file, 'a+')
		lines = []
		total = len(palindromes)  # 用于显示进度
		for start, end, seq0, length, types in tqdm(
			palindromes, total=total, desc="Writing GFF", unit="record"
		):
			fields = [
				chromosome_name,
				"CentriVision",
				str(types),
				str(start + 1),
				str(end),
				".",
				".",
				".",
				f"Length={length};seq={seq0}"
			]
			lines.append("\t".join(fields))
		fh.write("\n".join(lines) + "\n")
		fh.close()

	def readseq(self,seq1,seq2):
		d = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
		seq2_ = ''.join([d[i] for i in seq2])
		if seq1 == seq2:
			v = 'EQ'
			p = True
		elif seq1 == seq2[::-1]:
			v = 'PS'
			p =True
		elif seq1 == seq2_:
			v = 'RC'
			p =True
		else:
			v = ''
			p = False
		return p,v

	def find_palindromesL(self,chromosome_name,sequence, queue=None):
		seqlength = len(sequence)
		start,palindromes = [],[]
		last_percent = -1
		for i in range(0,seqlength-(2*self.length),self.length):
		# for i in trange(seqlength-(2*self.length)):
			seq1 = sequence[i:i+self.length]
			if 'N' in seq1 or len(set(list(seq1))) == 1:
				continue
			for j in range(i+self.length,min(i+self.reach,seqlength-self.length)):
				seq2 = sequence[j:j+self.length]
				if 'N' in seq2 or len(set(list(seq2))) == 1:
					continue
				p,v = self.readseq(seq1,seq2)
				if p:
					if i not in start:
						start.append(i)
						palindromes.append((i, i+self.length,seq1, self.length,v))
					if j not in start:
						start.append(j)
						palindromes.append((j, j+self.length,seq2, self.length,v))

			if queue is not None:
				percent = int(100 * i / seqlength)
				if percent != last_percent:
					queue.put(("chrom_progress", i / seqlength))
					last_percent = percent
		if queue is not None:
			queue.put(("chrom_done", 1))
		return palindromes

	def run_palindrome_c_with_queue_chunk(self, sequence: str, k: int, reach: int, queue=None, libpath: str = None, chunk_size: int = 500000):
		import ctypes
		import numpy as np
		import os

		base_dir = os.path.dirname(__file__)
		if libpath is None:
			libpath = os.path.join(base_dir, "Lib", "palindrome.so") if os.name != "nt" else os.path.join(base_dir, "Lib", "palindrome.dll")

		lib = ctypes.CDLL(libpath)
		lib.find_palindromes_chunk.argtypes = [
			ctypes.c_char_p,
			ctypes.c_int,
			ctypes.c_int,
			ctypes.c_int,
			ctypes.c_int,
			ctypes.c_int,
			ctypes.POINTER(ctypes.c_int)
		]
		lib.find_palindromes_chunk.restype = ctypes.c_int

		seqlength = len(sequence)
		results = (ctypes.c_int * (chunk_size * 5))()
		palindromes = []

		last_percent = -1  # 上一次发送的百分比
		for start_idx in range(0, seqlength, chunk_size):
			end_idx = min(start_idx + chunk_size, seqlength)
			count = lib.find_palindromes_chunk(sequence.encode("utf-8"),
											   seqlength,
											   k,
											   reach,
											   start_idx,
											   end_idx,
											   results)
			if count > 0:
				arr = np.ctypeslib.as_array(results, shape=(count * 5,))
				arr = arr[:count*5].reshape(count, 5)
				for r in arr:
					s, e, t, length, _ = r
					if t == 1: v = "EQ"
					elif t == 2: v = "PS"
					elif t == 3: v = "RC"
					else: v = ""
					palindromes.append((s, e, sequence[s:e], length, v))

			if queue is not None:
				percent = int(100 * end_idx / seqlength)
				if percent != last_percent and percent % 1 == 0:  # 每 1% 发送一次
					queue.put(("chrom_progress", min(1.0, end_idx / seqlength)))
					last_percent = percent

		start_set = set()
		dedup_palindromes = []
		for s, e, seq0, length, v in palindromes:
			if s not in start_set:
				start_set.add(s)
				dedup_palindromes.append((s, e, seq0, length, v))

		if queue is not None:
			queue.put(("chrom_done", 1))

		return dedup_palindromes



	def line_plot(self,data,length):
		x,y = [],[]
		start,end = -self.step,self.windows-self.step
		while end < length:
			start,end = start+self.step,end+self.step
			index = (start+end)/2
			data0 = 100*sum(data[start:end])/self.windows
			x.append(index)
			y.append(data0)
		y_max = max(y)
		return x,y,y_max

	def covarage(self,local,length,chro):
		array_dall = np.zeros(length, dtype = int)
		array_deq = np.zeros(length, dtype = int)
		array_dps = np.zeros(length, dtype = int)
		array_drc = np.zeros(length, dtype = int)
		local[3] = local[3]-1
		for index,row in local.iterrows():
			array_dall[row[3]:row[4]] = array_dall[row[3]:row[4]]+1
			if row[2] == 'EQ':
				array_deq[row[3]:row[4]] = array_deq[row[3]:row[4]]+1
			elif row[2] == 'PS':
				array_dps[row[3]:row[4]] = array_dps[row[3]:row[4]]+1
			elif row[2] == 'RC':
				array_drc[row[3]:row[4]] = array_drc[row[3]:row[4]]+1
		array_dall1 = np.int64(array_dall>0)
		xall1,yall1,yall_max = self.line_plot(array_dall1,length)
		array_deq1 = np.int64(array_deq>0)
		xeq1,yeq1,yeq_max = self.line_plot(array_deq1,length)
		array_dps1 = np.int64(array_dps>0)
		xps1,yps1,yps_max = self.line_plot(array_dps1,length)
		array_drc1 = np.int64(array_drc>0)
		xrc1,yrc1,yrc_max = self.line_plot(array_drc1,length)
		return (xall1,yall1,yall_max,xeq1,yeq1,yeq_max,xps1,yps1,yps_max,xrc1,yrc1,yrc_max)

	def vision(self,lens):

		# 读取CSV文件
		df = pd.read_csv(self.gff_file,header = None, sep='\t', comment='#').sort_values(by=[0,3],ascending= [True,True])
		print(df)
		# exit()
		# 设置绘图样式
		sns.set(style="whitegrid")

		# 创建一个新的图形和子图
		fig = plt.figure(figsize=(self.width, self.height))

		if (df[0].nunique() / self.coln) - (df[0].nunique() // self.coln) > 0:
			row = df[0].nunique() // self.coln + 1
		else:
			row = df[0].nunique() // self.coln
		gs = GridSpec(row, self.coln, figure=fig)

		# 分别绘制不同染色体的重复序列长度分布
		for i, (chro, data) in enumerate(df.groupby(0)):
			row = (i // self.coln)
			col = i % self.coln
			# print(row,col)
			ax = fig.add_subplot(gs[row, col])
			local = self.covarage(data,lens[chro],chro)
			sns.lineplot(x='index', y='covarage', data={'index':local[0],'covarage':local[1]}, color='red', marker='d', alpha=1, label='ALL')
			sns.lineplot(x='index', y='covarage', data={'index':local[3],'covarage':local[4]}, color='black', marker='o',linestyle='dashdot', alpha=0.6, label='EQ')
			sns.lineplot(x='index', y='covarage', data={'index':local[6],'covarage':local[7]}, color='blue', marker='+',linestyle='dashdot', alpha=0.6, label='PS')
			sns.lineplot(x='index', y='covarage', data={'index':local[9],'covarage':local[10]}, color='green', marker='^',linestyle='dashdot', alpha=0.6, label='RC')
			# 设置轴标签
			ax.set_xlabel(chro + 'index (bp)')
			ax.set_ylabel('covarage (%)')

			ax.set_title(chro+' Length: '+str(self.length))
			# if i == len(lens.keys())-1:
			#	 plt.legend()

		# 调整子图之间的间距
		plt.tight_layout()

		ax.legend()
		plt.savefig(self.savefile, dpi = self.dpi, bbox_inches = 'tight')

	# 多进程错误打印
	def print_error(self,value):
		print("Process pool error, the cause of the error is :", value)


	# def run(self):
	# 	with open(self.gff_file, 'w') as f:  # 清空或创建输出文件
	# 		f.write("") 
	# 	lens = {}
	# 	pool = multiprocessing.Pool(processes = self.cpu) # Step II : 进程池
	# 	result = []
	# 	for record in SeqIO.parse(self.genome_file, "fasta"):
	# 		chromosome_name = record.id
	# 		chromosome_sequence = str(record.seq)
	# 		lens[chromosome_name] = len(chromosome_sequence)
	# 		# palindromes = self.find_palindromesL(chromosome_name, chromosome_sequence.upper())
	# 		if self.code == 'c':
	# 			# result.append([chromosome_name,pool.apply_async(self.run_palindrome_c, args=(chromosome_sequence.upper(),self.length,self.reach))])
	# 			result.append([chromosome_name,pool.apply_async(self.run_palindrome_c_with_progress, args=(chromosome_sequence.upper(),self.length,self.reach))])
	# 		else:
	# 			result.append([chromosome_name,pool.apply_async(self.find_palindromesL, args=(chromosome_name, chromosome_sequence.upper()))])

	# 	pool.close()
	# 	pool.join()
	# 	for item in result:
	# 		chromosome_name, palindromes = item[0],item[1].get()
	# 		self.write_gff_file(chromosome_name, palindromes)
	# 	self.vision(lens)


	def run(self):
		from multiprocessing import Pool, Manager
		from Bio import SeqIO
		from tqdm import tqdm
		import time

		# 清空输出文件
		with open(self.gff_file, 'w') as f:
			f.write("")

		chromosomes = list(SeqIO.parse(self.genome_file, "fasta"))
		total = len(chromosomes)
		lens = {rec.id: len(rec.seq) for rec in chromosomes}

		manager = Manager()
		queue = manager.Queue()
		pool = Pool(processes=self.cpu)
		async_results = []

		for record in chromosomes:
			async_results.append(pool.apply_async(
				worker_chromosome,
				args=(record, self.length, self.reach, self.code, queue, None),
				error_callback=self.print_error
			))

		pbar_genome = tqdm(total=total, desc="Genome", position=0)
		pbar_chr = None
		finished = 0

		chrom_name_current = ""
		while finished < total:
			try:
				msg = queue.get(timeout=1)
			except:
				continue
			if msg[0] == "chrom_progress":
				if pbar_chr is None:
					pbar_chr = tqdm(total=1.0, desc=f"Chromosome", position=1, leave=False)
				pbar_chr.n = msg[1]
				pbar_chr.refresh()
			elif msg[0] == "chrom_done":
				if pbar_chr is not None:
					pbar_chr.n = 1.0
					pbar_chr.refresh()
					pbar_chr.close()
					pbar_chr = None
				pbar_genome.update(1)
				finished += 1
			elif msg[0] == "error":
				print("Subprocess error:\n", msg[1])
			time.sleep(0.01)

		pool.close()
		pool.join()
		print('write GFF!')
		# 写入 GFF
		for res in async_results:
			chrom_name, palindromes = res.get()
			print(chrom_name)
			self.write_gff_file(chrom_name, palindromes)

		# 可视化
		self.vision(lens)
