from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
# import Levenshtein
from rapidfuzz.distance import Levenshtein
import itertools
import multiprocessing
from tqdm import tqdm
import CentriVision.bez as bez


class EdistAln():
	def __init__(self, options):
		self.window = 2000
		self.cpu = 8
		bez_conf = bez.config()
		for k, v in bez_conf:# 
			setattr(self, str(k), v)
		for k, v in options:
			setattr(self, str(k), v)
			print(k, ' = ', v)
		self.window = int(self.window)
		self.cpu = int(self.cpu)

		self.gff = 'centri.gff'
		self.lens = 'centri.lens'
		self.fasta = 'centri.fasta'

	def compute_similarity(self,seqA, seqB):
		"""计算 A vs B 和 A'（反向互补） vs B 的相似度，并返回匹配方向"""
		# dist_AB = Levenshtein.distance(seqA, seqB)
		dist_AB = Levenshtein.normalized_similarity(seqA, seqB)
		# rev_comp_A = str(Seq(seqA).reverse_complement())
		rev_comp_A = seqA[::-1].translate(str.maketrans("ATCGatcg", "TAGCtagc"))
		# dist_ApB = Levenshtein.distance(rev_comp_A, seqB)
		dist_ApB = Levenshtein.normalized_similarity(rev_comp_A, seqB)
		
		max_len = max(len(seqA), len(seqB))
		sim_AB = 1 - dist_AB / max_len
		sim_ApB = 1 - dist_ApB / max_len

		if sim_AB >= sim_ApB:
			return sim_AB, "forward"
		else:
			return sim_ApB, "reverse_complement"

	def compare_pair(self,pair):
		"""对给定的两个序列计算相似度"""
		(id1, seq1), (id2, seq2) = pair
		similarity, direction = self.compute_similarity(seq1, seq2)
		return (id1, id2, similarity, direction)

	def process_fasta(self,input_fasta1, input_fasta2=None, output_file="output.tsv", num_processes=4):
		"""读取FASTA文件，计算相似度，并排序后输出结果"""
		
		# 读取第一个基因组
		sequences1 = {record.id: str(record.seq) for record in SeqIO.parse(input_fasta1, "fasta") if len(record.seq) >= 10}

		if input_fasta2:  
			# **两个基因组之间比对**
			sequences2 = {record.id: str(record.seq) for record in SeqIO.parse(input_fasta2, "fasta") if len(record.seq) >= 10}
			sequence_pairs = itertools.product(sequences1.items(), sequences2.items())
			total_pairs = len(sequences1) * len(sequences2)
		else:
			# **单基因组内比对**
			sequence_pairs = itertools.combinations(sequences1.items(), 2)
			from math import comb
			total_pairs = comb(len(sequences1), 2)

		# 边计算边写文件，避免堆积结果
		with open(output_file, "w") as f:
			f.write("Seq1\tSeq2\tSimilarity\tDirection\n")
			with multiprocessing.Pool(num_processes) as pool:
				with tqdm(total=total_pairs, desc="Comparing Pairs", unit="pair") as pbar:
					for id1, id2, similarity, direction in pool.imap_unordered(self.compare_pair, sequence_pairs, chunksize=100):
						f.write(f"{id1}\t{id2}\t{similarity:.4f}\t{direction}\n")
						pbar.update()

	def cat_fasta(self):
		# 切割片段长度
		fragment_size = self.window
		lens_records = []
		gff_records = []
		split_records = []

		# 读取序列并进行切割
		for record in SeqIO.parse(self.centri_sequence, 'fasta'):
			chrom = record.id
			seq = str(record.seq)
			chrom_length = len(seq)
			fragment_count = (chrom_length + fragment_size - 1) // fragment_size  # 计算切片个数

			lens_records.append(f"{chrom}\t{chrom_length}\t{fragment_count}")

			for i in range(fragment_count):
				start = i * fragment_size + 1
				end = min((i + 1) * fragment_size, chrom_length)
				fragment_name = f"{chrom}_part{i + 1}"
				gff_records.append(f"{chrom}\t{fragment_name}\t{start}\t{end}\t+\t{i + 1}\t{fragment_name}")
				fragment_seq = seq[start-1:end]
				split_records.append(SeqRecord(Seq(fragment_seq), id=fragment_name, description=""))

		# 写入lens文件
		with open(self.lens, 'w') as lens_out:
			lens_out.write('\n'.join(lens_records) + '\n')

		# 写入gff文件
		with open(self.gff, 'w') as gff_out:
			gff_out.write('\n'.join(gff_records) + '\n')

		# 写入切片fasta文件
		SeqIO.write(split_records, self.fasta, 'fasta')

		print("文件生成完成！")
	def run(self):
		print('加速版！')
		self.cat_fasta()
		import time
		start_time = time.time()
		self.process_fasta(self.fasta, output_file=self.out_file, num_processes=self.cpu)
		# ... 你的代码运行 ...
		end_time = time.time()

		print(f"运行时间: {end_time - start_time:.2f} 秒")