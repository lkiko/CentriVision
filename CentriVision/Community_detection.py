import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import gc
import multiprocessing # Step I : 导入模块
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import random
import os
import shutil
import CentriVision.bez as bez
from multiprocessing import cpu_count#读取CPU核心数用于匹配线程数

class Community_detection():
	def __init__(self, options):
		self.workpath = os.getcwd()+'/'
		self.cpu = cpu_count()
		self.gap = '10'
		self.identity = '75'
		self.alignment = '75'
		self.convergence_p = 'convergence.png'
		self.blast = 'True'
		bez_conf = bez.config()
		for k, v in bez_conf:# 
			setattr(self, str(k), v)
		for k, v in options:
			setattr(self, str(k), v)
			print(str(k), ' = ', v)
		self.cpu = int(self.cpu)
		self.gap = int(self.gap)
		self.identity = int(self.identity)
		self.alignment = int(self.alignment)
		

	def communitys_found(self,communitys0,Round_,code_convergence,threadi,bodys):
		# global thread,bodys
		communitys = communitys0
		print(len(communitys))
		Rounds,number,convergence,convergences = [],[],[],[]
		Round = 0
		while 1:
			Round += 1
			if Round_:
				job = int(len(communitys)/2)
			else:
				job = len(communitys)
			for i in range(job):
				if (Round_ and threadi == self.cpu - 1) or not Round_:
					p = int(((i+1)/job)*50)
					print("\r["+"*"*p+" "*(50-p)+"] "+str(i+1)+"/"+str(job)+' Round '+str(Round)+' ',end="")
				lt = np.linspace(0,len(communitys)-1,len(communitys),dtype='int')
				lt_ = np.random.choice(lt,2,replace=False)
				if len(list(set(communitys[lt_[0]]) & set(communitys[lt_[1]]))) > 0:
					lt1,lt2 = communitys[lt_[0]], communitys[lt_[1]]
					new_community = lt1 + lt2
					new_community = list(set(new_community))
					communitys.remove(lt1)
					communitys.remove(lt2)
					communitys.append(new_community)
			Rounds.append(Round)
			number.append(len(communitys))
			if len(convergence)<10:
				convergence.append(len(communitys))
				if (Round_ and threadi == self.cpu - 1) or not Round_:
					print(str(len(communitys))+' -> std -> 000')
			else:
				del convergence[0]
				convergence.append(len(communitys))
				if np.std(convergence) <= code_convergence:
					convergences.append(np.std(convergence))
					if Round_:
						print('\033[0;31;42m'+str(len(communitys))+'\033[0m -> std -> '+str(int(np.std(convergence))))
						break
					else:
						bodys2 = []
						for i in communitys:
							bodys2 = bodys2 + i
						print('communitys > \033[0;31;42m'+str(len(communitys))+'\033[0m -> std > '+str(int(np.std(convergence)))+' -> bodys > '+str(len(bodys))+' -> bodys2 > \033[0;31;42m'+str(len(bodys2))+'\033[0m')
						if len(bodys2) == len(bodys):
							print('\ncommunitys -> ',len(communitys))
							break
				else:
					convergences.append(np.std(convergence))
					if (Round_ and threadi == self.cpu - 1) or not Round_:
						print('\033[0;31;42m'+str(len(communitys))+'\033[0m -> std -> '+str(int(np.std(convergence))))
					pass
		if Round_:
			return communitys
		else:
			return communitys,Rounds,number,Rounds[10:],convergences

	def get_length(self,seq1,seq_):
		return len(str(seq_[seq1].seq))

	def read_blast(self,blast,seq_,number):
		blast[12] = blast.apply(lambda x: self.get_length(x[0],seq_), axis=1)
		blast[13] = blast.apply(lambda x: self.get_length(x[1],seq_), axis=1)
		# blast[12] = blast[0].map(self.get_length)# 序列长度
		# blast[13] = blast[1].map(self.get_length)
		blast[14] = 100*(blast[7]-blast[6])/blast[12] #1比对上占比
		blast[15] = 100*(blast[9]-blast[8])/blast[13]
		blast[16] = 100*(blast[5])/blast[3] # gap占比
		blast = blast.drop(blast[(blast[14] < self.alignment) & (blast[15] < self.alignment) | (blast[0] == blast[1]) | (blast[2] < self.identity) | (blast[16] > self.gap)].index)# The operators are: | for or, & for and, and ~ for not.
		del blast[2]
		del blast[3]
		del blast[4]
		del blast[5]
		del blast[6]
		del blast[7]
		del blast[8]
		del blast[9]
		del blast[10]
		del blast[11]
		del blast[12]
		del blast[13]
		del blast[14]
		del blast[15]
		del blast[16]
		# gc.collect()
		print('***************************************** .%s end. *****************************************' % number)
		return blast


	def get_communitys(self,blast):
		
		communitys = [] #社区
		bodys = list(set(list(blast[0].to_list())+list(blast[1].to_list())))
		dic0 = blast.groupby(0).groups
		for key in dic0.keys():
			# print(" ******** ",key," ******** ")
			local = blast.loc[dic0[key]].sort_values(by=[1],ascending= [True])
			local.reset_index(drop=True,inplace=True)
			# print(local)
			communitys_local = list(set(list(local[0].to_list())+list(local[1].to_list())))
			communitys.append(communitys_local)
		return communitys,bodys

	def parse_sam(self,file):
		print('''There are sequences with a length of less than or equal to 20bp in the file.
Bowtie2 is used for alignment by default. 
Please ensure that the 'pysam' module is correctly installed or complete the 
installation using the command : 
			pip install pysam
Additionally, please configure Bowtie2.''')
		import pysam
		samfile = pysam.AlignmentFile(file, "r")

		# 存储比对结果
		results = []

		# 遍历每个比对条目
		for read in samfile.fetch():
			if not read.is_unmapped:
				# 获取查询序列和目标序列
				query = read.query_name
				target = samfile.getrname(read.reference_id)
				
				# 获取比对相关的信息
				query_start = read.query_alignment_start + 1
				query_end = read.query_alignment_end
				target_start = read.reference_start + 1
				target_end = read.reference_end

				# 获取比对分数和匹配信息
				alignment_length = read.query_alignment_length
				mismatches = read.get_tag("NM")
				identity = (alignment_length - mismatches) / alignment_length if alignment_length > 0 else 0
				gap_opens = 0  # pysam 不直接提供 gap 信息，假设为 0，或根据需要自定义

				# 假设 e-value 和 bit score （或从其他信息中推导）
				e_value = 0.0  # 可根据需求调整
				bit_score = read.mapping_quality

				# 将结果存储为类 BLAST 格式
				result = {
					"query": query,
					"target": target,
					"identity": round(identity * 100, 2),  # 百分比表示
					"alignment_length": alignment_length,
					"mismatches": mismatches,
					"gap_opens": gap_opens,
					"query_start": query_start,
					"query_end": query_end,
					"target_start": target_start,
					"target_end": target_end,
					"e_value": e_value,
					"bit_score": bit_score
				}
				results.append(result)
		
		# 关闭 SAM 文件
		samfile.close()
		
		return results

	def print_blast_format(self,results):
		# 打印符合 BLAST 标准的11列格式
		f = open(self.fasta_file + ".blast",'w')
		for res in results:
			s = str(res['query'])+'\t'+str(res['target'])+'\t'+str(res['identity'])+'\t'+\
			str(res['alignment_length'])+'\t'+str(res['mismatches'])+'\t'+str(res['gap_opens'])+'\t'+\
			str(res['query_start'])+'\t'+str(res['query_end'])+'\t'+str(res['target_start'])+'\t'+\
			str(res['target_end'])+'\t'+str(res['e_value'])+'\t'+str(res['bit_score'])
			f.write(s+'\n')
		f.close()

	def run(self):
		def old_rev(input_fasta,output_fasta):
			# 读取原始FASTA文件，并生成新的序列
			new_records = []
			for record in SeqIO.parse(input_fasta, "fasta"):
				name = str(record.id)
				# 原序列
				original_record = record[:]
				original_record.id = name + "_old"
				original_record.description = ""
				new_records.append(original_record)
				# 反向互补序列
				rev_comp_record = record[:]
				rev_comp_record.id = name + "_rev"
				rev_comp_record.seq = record.seq.reverse_complement()
				rev_comp_record.description = ""
				new_records.append(rev_comp_record)
			# 写入新的FASTA文件
			SeqIO.write(new_records, output_fasta, "fasta")

		code_convergence1 = 0 #第一轮的标准差
		code_convergence2 = 0
		seq_0 = SeqIO.to_dict(SeqIO.parse(self.fasta_file, "fasta"))
		old_rev(self.fasta_file,self.fasta_file+'.new.fasta')
		self.fasta_file = self.fasta_file+'.new.fasta'
		seq_ = SeqIO.to_dict(SeqIO.parse(self.fasta_file, "fasta"))# 提取之后直接返回字典
		minseq = min([len(seq_[seq0].seq) for seq0 in seq_.keys()])

		if self.blast == 'True':
			if minseq > 20:
				d = os.popen(self.blast_path+"/makeblastdb -in %s -dbtype nucl -out ./%s/%s" % (self.fasta_file,'db',self.fasta_file+'.db')).read().strip()
				d = os.popen(self.blast_path+"/blastn -query %s -db ./%s/%s -out %s -outfmt 6 -evalue 1e-5 -num_threads %s" % (self.fasta_file,'db',self.fasta_file+'.db',self.fasta_file + ".blast",str(self.cpu))).read().strip()
				shutil.rmtree('db')
			else:
				os.makedirs('db')
				d = os.popen(self.bowtie2_path+"/bowtie2-build %s ./db/%s.index" % (self.fasta_file,self.fasta_file)).read().strip()
				d = os.popen(self.bowtie2_path+"/bowtie2 -x ./db/%s.index -f %s -S ./db/output.sam --threads %s -N %s --very-sensitive --all" % (self.fasta_file,self.fasta_file,str(self.cpu),'1')).read().strip()
				results = self.parse_sam("./db/output.sam")
				self.print_blast_format(results)
				shutil.rmtree('db')

		blast = pd.read_csv(self.fasta_file + ".blast",header = None, sep='\t',chunksize=200000)

		result = []
		local_n = 0
		pool = multiprocessing.Pool(processes = self.cpu) # Step II : 进程池
		for local in blast:
			local[0] = local[0].str[:-4]
			local[1] = local[1].str[:-4]
			local_n += 1
			print('***************************************** local number of %s *****************************************' % str(local_n))
			result.append(pool.apply_async(self.read_blast, (local,seq_0,str(local_n)), ))  # Step III : 异步（并行）计算 返回值！！！！
		pool.close() # Step IV : 准备结束
		pool.join() # Step IV : 完全结束
		del seq_,seq_0
		gc.collect()
		blast = pd.DataFrame()
		for item in result:
			blast = pd.concat([blast, item.get()])
		communitys,bodys = self.get_communitys(blast)
		# print(bodys)
		
		# print(len(communitys))
		print("###################### Step 1 ######################")
		one = int(len(communitys)/self.cpu)
		result = []
		pool = multiprocessing.Pool(processes = self.cpu) # Step II : 进程池
		for i in range(self.cpu):#################
			if i == self.cpu - 1:
				communitys_ = communitys[i*one:]
			else:
				communitys_ = communitys[i*one:(i+1)*one]
			result.append(pool.apply_async(self.communitys_found, (communitys_,True,code_convergence1,i,bodys), ))  # Step III : 异步（并行）计算 返回值！！！！
		pool.close() # Step IV : 准备结束
		pool.join() # Step IV : 完全结束
		communitys = []
		for i in result:

			communitys = communitys + list(i.get())

		print(len(communitys))
		print("###################### Step 2 ######################")
		one = int(len(communitys)/self.cpu)
		result = []
		pool = multiprocessing.Pool(processes = self.cpu) # Step II : 进程池
		for i in range(self.cpu):#################
			if i == self.cpu - 1:
				communitys_ = communitys[i*one:]
			else:
				communitys_ = communitys[i*one:(i+1)*one]
			result.append(pool.apply_async(self.communitys_found, (communitys_,True,code_convergence1,i,bodys,), ))  # Step III : 异步（并行）计算 返回值！！！！
		pool.close() # Step IV : 准备结束
		pool.join() # Step IV : 完全结束
		communitys = []
		for i in result:
			# print(len(i.get()))
			communitys = communitys + i.get()

		print('next step jobe -- >',len(communitys))
		print("\n###################### Step 3 ######################")
		communitys,Rounds2,number2,Rounds_2,convergences2 = self.communitys_found(communitys,False,code_convergence2,None,bodys)
		fig,ax1 = plt.subplots()
		ax2 = ax1.twinx()
		ax1.plot(Rounds2,number2,'g-', label = 'Num_C')
		ax2.plot(Rounds_2,convergences2,'b--', label = 'Std')
		ax1.set_xlabel('Step2')	#设置x轴标题
		ax1.set_ylabel('Number of Community',color = 'g')   #设置Y1轴标题
		ax2.set_ylabel('Standard Deviation',color = 'b')   #设置Y2轴标题
		plt.legend()
		plt.savefig(self.convergence_p,dpi = 1000, bbox_inches='tight')
		plt.clf()
		out = open(self.out_file,'w')
		communitys1 = []
		for i in communitys:
			i.append(len(i))
			communitys1.append(i)
		sorted_communitys = sorted(communitys1, key=lambda x:x[-1], reverse=True)
		num = 0
		for i in sorted_communitys:
			out.write('community_'+str(num+1)+'\t'+'\t'.join(i[:-1])+'\n')
			num += 1
		out.close()
		# print('\n',len(bodys),len(bodys2))
		print('\n code end!')