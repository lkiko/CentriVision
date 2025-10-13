import pandas as pd
from tqdm import tqdm
import CentriVision.bez as bez
import os

class Get_centgff():
	def __init__(self, options):
		# self.centromere_file = 'new-myai-extend.txt'
		# self.gff_file = 'trf.gff'
		self.start = 3
		self.end = 4
		self.locmin = 9
		# self.output_file = 'cent-trf.gff'
		for k, v in options:
			setattr(self, str(k), v)
			print(str(k), ' = ', v)
		self.start = int(self.start)
		self.end = int(self.end)
		self.locmin = int(self.locmin)

	def load_centromere_positions(self,centromere_file):
		"""
		加载着丝粒位置文件到 DataFrame 中。
		"""
		centromeres = pd.read_csv(centromere_file, sep='\t', comment='#', header=None)
		return centromeres

	def filter_gff_in_centromere(self,gff_file, centromeres, start_col, end_col,locmin):
		"""
		过滤 GFF 文件中的行，保留仅位于着丝粒区域内的行，并显示进度条。
		参数：
		- gff_file: GFF 文件路径
		- centromeres: 着丝粒位置的 DataFrame
		- start_col: GFF 文件中起始位置的列索引（从 0 开始）
		- end_col: GFF 文件中终止位置的列索引（从 0 开始）
		"""
		filtered_lines = []

		grouped_centromeres  = centromeres.groupby(0)

	 # 读取 GFF 文件
		gff_df = pd.read_csv(gff_file, sep='\t', comment='#', header=None)
		# 按染色体进行分组
		# 使用 tqdm 显示进度条
		for chrom, group in tqdm(grouped_centromeres , desc="Processing Chromosomes", unit="chromosome"):

			gff_chrom = gff_df[gff_df[0] == chrom]
			# 如果没有对应的着丝粒区域，跳过
			if gff_chrom.empty:
				continue

			# 对当前组进行过滤
			for _, row in group.iterrows():
				centromere_start  = row[1]
				centromere_end = row[2]
				filtered_df = gff_chrom[(gff_chrom[start_col] >= centromere_start) & (gff_chrom[end_col] <= centromere_end)].copy()

				# 如果过滤后的 DataFrame 不为空
				if not filtered_df.empty:
					col = filtered_df.shape[1]
					filtered_df.loc[:, col] = filtered_df[start_col]
					filtered_df.loc[:, col+1] = filtered_df[end_col]
					filtered_df.loc[:, start_col] = filtered_df[start_col] - centromere_start
					filtered_df.loc[:, end_col] = filtered_df[end_col] - centromere_start
					filtered_df.to_csv(self.output_file, mode='a', header=False,index=False,sep='\t')


	def save_filtered_gff(self,output_file, filtered_lines):
		"""
		将过滤后的 GFF 行保存到输出文件中。
		"""
		with open(output_file, 'w') as out:
			for line in filtered_lines:
				out.write(line)

	def run(self):
		try:
			if os.path.exists(self.output_file):
				os.remove(self.output_file)
		except Exception as e:
			print(f"err...: {e}")

		# 加载着丝粒位置文件
		centromeres = self.load_centromere_positions(self.centromere_file)
		# print(centromeres)
		# 过滤 GFF 文件
		filtered_lines = self.filter_gff_in_centromere(self.gff_file, centromeres,self.start,self.end,self.locmin)
