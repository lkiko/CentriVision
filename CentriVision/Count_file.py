import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import CentriVision.bez as bez
import os

class Count_file():
	def __init__(self, options):
		self.workpath = os.getcwd() + '/'
		self.lmmin = 8
		self.lmmax = 1000
		self.bin_size = 15
		self.y_break_min = 0
		self.y_break_max = 0
		self.peak_index = 1
		self.peak_indices = ''
		
		# 添加字体大小参数
		self.title_fontsize = 22
		self.label_fontsize = 18
		self.tick_fontsize = 15
		self.annotation_fontsize = 15

		bez_conf = bez.config()
		for k, v in bez_conf:
			setattr(self, str(k), v)
		for k, v in options:
			setattr(self, str(k), v)
			print(str(k), ' = ', v)
		self.lmmin = int(self.lmmin)
		self.lmmax = int(self.lmmax)
		self.bin_size = int(self.bin_size)
		self.peak_index = max(int(self.peak_index), 1)
		self.y_break_max = max(int(self.y_break_max), 0)
		self.y_break_min = max(int(self.y_break_min), 0)
		self.peak_indices = [int(i) for i in self.peak_indices.split(',') if i] if self.peak_indices != 'None' else []
		self.title_fontsize = int(self.title_fontsize)
		self.label_fontsize = int(self.label_fontsize)
		self.tick_fontsize = int(self.tick_fontsize)
		self.annotation_fontsize = int(self.annotation_fontsize)

	def run(self):
		# 读取和处理数据
		df = pd.read_csv(self.dot_file, sep='\t')
		cols_to_drop = [f'M{i}' for i in range(1, 201)]
		df = df.drop(columns=[col for col in cols_to_drop if col in df.columns])
		df.to_csv(self.out_file, sep='\t', index=False)

		# 过滤数据
		filtered_df = df[(df['LM(Length of monomers)'] >= self.lmmin) & 
						 (df['LM(Length of monomers)'] <= self.lmmax)].copy()
		filtered_df['repeat_number'] = (filtered_df['length'] / 
										filtered_df['LM(Length of monomers)']).round().astype(int)
		result_array = np.repeat(filtered_df['LM(Length of monomers)'].values, 
								 filtered_df['repeat_number'].values)

		# 创建直方图
		bins = np.arange(self.lmmin, self.lmmax, self.bin_size)
		counts, edges = np.histogram(result_array, bins=bins)

		# 标记峰值
		peak_indices = np.argsort(counts)[-self.peak_index:]
		if self.peak_indices:
			peak_indices = np.append(peak_indices, self.peak_indices)
		
		bar_colors = ['#949495'] * len(counts)
		for peak_index in peak_indices:
			bar_colors[peak_index] = 'red'

		if self.y_break_min == 0:
			# 创建普通图表
			# plt.figure(figsize=(12, 6))
			plt.figure(figsize=(5, 5))
			plt.bar(edges[:-1], counts, width=np.diff(edges), 
					color=bar_colors, align='edge', edgecolor=None)
			
			ax = plt.gca()
			xlim = ax.get_xlim()
			ylim = ax.get_ylim()

			# 标记峰值
			for peak_index in peak_indices:
				peak_value = counts[peak_index]
				peak_position = (edges[peak_index] + edges[peak_index + 1]) / 2
				peak_range = f'[{edges[peak_index]}, {edges[peak_index + 1]}]'
				
				if peak_position < (xlim[0] + xlim[1]) / 2:
					a = +max((xlim[1] - xlim[0]) * 0.05, 5)
				else:
					a = -max((xlim[1] - xlim[0]) * 0.05, 5)
					
				if peak_value < (ylim[0] + ylim[1]) / 2:
					b = +max((ylim[1] - ylim[0]) * 0.1, 5)
				else:
					b = -max((ylim[1] - ylim[0]) * 0.1, 5)
					
				x_text = min(max(peak_position + a, xlim[0] + 5), xlim[1] - 5)
				y_text = min(max(peak_value + b, ylim[0] + 5), ylim[1] - 5)
				
				plt.annotate(f'Peak: {peak_value}\nRange: {peak_range}', 
							 xy=(peak_position, peak_value), 
							 xytext=(x_text, y_text),
							 arrowprops=dict(facecolor='black', arrowstyle="->"), 
							 fontsize=self.annotation_fontsize)

			# 使用字体大小变量
			plt.title('Distribution of Repeat Unit Lengths (bp)', fontsize=self.title_fontsize)
			plt.xlabel('Repeat Unit Length (bp)', fontsize=self.label_fontsize)
			plt.ylabel('Frequency', fontsize=self.label_fontsize)
			plt.xticks(fontsize=self.tick_fontsize)
			plt.yticks(fontsize=self.tick_fontsize)

		else:
			# 创建断轴图表
			from matplotlib.gridspec import GridSpec
			fig = plt.figure(figsize=(10, 6))
			gs = GridSpec(2, 1, height_ratios=[1, 3], hspace=0.1)
			ax1 = fig.add_subplot(gs[0])
			ax2 = fig.add_subplot(gs[1], sharex=ax1)
			
			# 绘制直方图
			ax1.bar(edges[:-1], counts, width=np.diff(edges), 
					 color=bar_colors, align='edge', edgecolor=None)
			ax2.bar(edges[:-1], counts, width=np.diff(edges), 
					 color=bar_colors, align='edge', edgecolor=None)
			
			# 设置轴范围
			ax1.set_ylim(self.y_break_max, max(counts) * 1.1)
			ax2.set_ylim(0, self.y_break_min)
			
			# 设置轴样式
			ax1.spines['bottom'].set_visible(False)
			ax2.spines['top'].set_visible(False)
			ax1.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

			# 计算断裂符号的位置
			ax1_pos = ax1.get_position()
			ax2_pos = ax2.get_position()
			y_break = (ax1_pos.y0 + ax2_pos.y1) / 2  # 断裂符号的 Y 位置
			x_break = ax1_pos.x0  # 断裂符号的 X 位置（左对齐）

			# 添加断轴标记
			fig.text(x_break, y_break, '//', ha='center', va='center', fontsize=10, 
					 transform=fig.transFigure, rotation=90)
			
			# 标记峰值
			for peak_index in peak_indices:
				peak_value = counts[peak_index]
				peak_position = (edges[peak_index] + edges[peak_index + 1]) / 2
				peak_range = f'[{edges[peak_index]}, {edges[peak_index + 1]}]'
				
				if peak_value > self.y_break_max:
					ax = ax1
					ylim = ax1.get_ylim()
				else:
					ax = ax2
					ylim = ax2.get_ylim()
					
				xlim = ax.get_xlim()
				
				if peak_position < (xlim[0] + xlim[1]) / 2:
					a = +max((xlim[1] - xlim[0]) * 0.05, 5)
				else:
					a = -max((xlim[1] - xlim[0]) * 0.05, 5)
					
				if peak_value < (ylim[0] + ylim[1]) / 2:
					b = +max((ylim[1] - ylim[0]) * 0.1, 5)
				else:
					b = -max((ylim[1] - ylim[0]) * 0.1, 5)
					
				x_text = min(max(peak_position + a, xlim[0] + 5), xlim[1] - 5)
				y_text = min(max(peak_value + b, ylim[0] + 5), ylim[1] - 5)
				
				ax.annotate(f'Peak: {peak_value}\nRange: {peak_range}', 
							xy=(peak_position, peak_value), 
							xytext=(x_text, y_text),
							arrowprops=dict(facecolor='black', arrowstyle="->"), 
							fontsize=self.annotation_fontsize)
			
			# 使用字体大小变量
			ax1.set_title('Distribution of Repeat Unit Lengths (bp)', 
						  fontsize=self.title_fontsize)
			ax2.set_xlabel('Repeat Unit Length (bp)', fontsize=self.label_fontsize)
			ax1.tick_params(axis='both', labelsize=self.tick_fontsize)
			ax2.tick_params(axis='both', labelsize=self.tick_fontsize)
		plt.ylabel('Frequency', fontsize=self.label_fontsize)
		plt.savefig(self.savefile, dpi=1000, bbox_inches='tight')