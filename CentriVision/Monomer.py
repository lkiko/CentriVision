from Bio import SeqIO
import matplotlib.pyplot as plt
# import Levenshtein
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from Bio.Align import PairwiseAligner
import seaborn as sns
from matplotlib.colors import Normalize, ListedColormap,LinearSegmentedColormap
import os
import CentriVision.bez as bez

# 比对函数，使用Levenshtein距离
# 如果两个字符串完全相同，Levenshtein距离为0。
# 如果两个字符串完全不同，Levenshtein距离为两个字符串长度之和。

class Monomer():
    def __init__(self, options):
        self.workpath = os.getcwd()+'/'
        self.workpaths = './Monomer_output/'
        self.horfasta = 'Monomer.fasta'
        self.horgff = 'Monomer.gff'
        self.decision_fig = 'True'
        self.Decision_file = './Decision_plot/'
        self.seed = 321
        self.window = 20
        path = bez.get_path()
        font_path = os.path.join(path, 'example/arial.ttf')
        from matplotlib.font_manager import FontProperties
        self.font_prop = FontProperties(fname=font_path)
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)
        self.seed = int(self.seed)
        self.window = int(self.window)

    def levenshtein_numpy(self,str1, str2):
        m, n = len(str1), len(str2)
        dp = np.zeros((m + 1, n + 1), dtype=int)
        dp[:, 0] = np.arange(m + 1)
        dp[0, :] = np.arange(n + 1)
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                cost = 0 if str1[i - 1] == str2[j - 1] else 1
                dp[i, j] = min(dp[i - 1, j] + 1,  # 删除
                               dp[i, j - 1] + 1,  # 插入
                               dp[i - 1, j - 1] + cost)  # 替换
        return dp[m, n]

    def compare_sequences(self, seed, seq):
        return self.levenshtein_numpy(seed, seq)

    # 比较并绘图
    def compare_and_plot(self,seqneme, sequence, start,length0, window_start, window_end, step=1):
            # 提取前320bp作为种子序列
        seed_sequence = sequence[start:start+length0]
        distances = []
        indices = []
        # 从320bp后开始逐步提取片段，步长为step
        for i in range(start+length0, start+length0+1, step): # len(sequence) - window_end
            for length in range(window_start, window_end + 1):
                if i+length >= len(sequence):
                    return 0,0,0
                sub_sequence = sequence[i:i + length]
                distance = self.compare_sequences(seed_sequence, sub_sequence)
                # print(distance)
                if distance >= int((len(seed_sequence)+len(sub_sequence))*0.5):
                    continue
                distances.append(distance)
                indices.append(length)
        if len(distances) < 1:
            return 'None',0,0
        # print(distances)
        # 找到最近距离及其对应位置
        min_distance = min(distances)
        min_index = distances.index(min_distance)
        min_position = indices[min_index]
        # # 绘制折线图
        plt.plot(indices, distances)
        plt.xlabel('Position in Sequence', fontproperties=self.font_prop)
        plt.ylabel('Distance to Seed Sequence', fontproperties=self.font_prop)
        plt.title('Distance of Sequence Segments to Seed', fontproperties=self.font_prop)
        plt.savefig(self.workpath+self.workpaths+self.Decision_file+seqneme+'_'+str(i)+'score.png', dpi = 600, bbox_inches = 'tight')# 存储图片
        plt.close()
        return min_distance,min_position,i


    def calculate_similarity(self, seq1, seq2):
        """计算两个序列的相似度，返回相似度分值（0到1之间）"""
        # 使用 Needleman-Wunsch 算法进行全局比对
        aligner = PairwiseAligner()
        # 设置全局比对
        aligner.mode = 'global'
        # aligner.mode = 'local'  # 设置局部比对

        # 执行比对
        alignment = aligner.align(seq1, seq2)
        best_alignment = alignment[0]  # 获取最佳比对结果
        aligned_seq1 = best_alignment[0]
        aligned_seq2 = best_alignment[1]

        # 统计相同碱基的数量
        total_length = 0
        matches = 0

        for base1, base2 in zip(aligned_seq1, aligned_seq2):
            if base1 == '-' and base2 == '-':
                continue
            else:
                total_length += 1  # 计入总长度
                if base1 == base2:
                    matches += 1  # 计入匹配数量

        # 计算相似度分值
        similarity_score = matches / total_length if total_length > 0 else 0
        # print(similarity_score,matches,total_length)
        return similarity_score

    def value_to_color(self, value, cmap, vmin, vmax):
        norm = Normalize(vmin=vmin, vmax=vmax)
        return cmap(norm(value))

    def plot_rectangles(self, seqneme, start_positions, lengths, long_sequence_l, long_sequence):
        fig, ax = plt.subplots()
        # 提取第一个单元的序列
        first_unit_seq = long_sequence[start_positions[0]:start_positions[0] + lengths[0]]
        # 存储分值
        scores = []

        # 第一个循环：计算所有单元的相似度分值
        for start, length in zip(start_positions, lengths):
            # 切割出当前单元的序列
            current_unit_seq = long_sequence[start:start + length]
            
            # 计算相似度分值
            score = self.calculate_similarity(first_unit_seq, current_unit_seq)
            scores.append(score)

        colors = sns.color_palette("viridis", as_cmap=True)
        cmap = LinearSegmentedColormap.from_list("Custom", colors(np.linspace(0, 1, 256)))

        # 存储中点和归一化宽度
        midpoints = []
        normalized_widths = []

        # 第二个循环：根据分值绘制矩形
        for i, (start, length) in enumerate(zip(start_positions, lengths)):
            # 绘制矩形
            rect = patches.Rectangle((start, 0), length, 1, facecolor=self.value_to_color(scores[i],cmap,0,1), edgecolor=None, alpha=0.7)
            ax.add_patch(rect)

            # 计算中点和归一化宽度
            midpoints.append(start + length / 2)
            normalized_widths.append((length - min(lengths)) / (max(lengths) - min(lengths)))  # 归一化宽度
        
        # 绘制折线图
        plt.plot(midpoints, normalized_widths, marker='o', color='red', linestyle='-', linewidth=0.5, markersize=1, label='Normalized Widths')

        # 添加颜色条
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
        sm.set_array([])  # 使得颜色条可以显示
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label('Similarity Score')

        # 设置图像显示范围
        plt.xlim(0, long_sequence_l)  # 根据最大起始位置和长度设置X轴范围
        plt.xticks(rotation=45)
        plt.ylim(0, 1)  # Y轴范围设置为1.5，这里矩形块在Y=0的水平线上
        plt.xlabel('Position', fontproperties=self.font_prop)
        plt.ylabel('Repetitive Units', fontproperties=self.font_prop)
        plt.title('Repetitive Unit Visualization', fontproperties=self.font_prop)
        plt.savefig(self.workpath+self.workpaths+seqneme+'.png', dpi=600, bbox_inches='tight')  # 存储图片
        plt.close()

    def run(self):
        if not os.path.exists(self.workpath+self.workpaths):
            os.makedirs(self.workpath+self.workpaths)
        if not os.path.exists(self.workpath+self.workpaths+self.Decision_file) and self.decision_fig.upper() == 'TRUE': 
            os.makedirs(self.workpath+self.workpaths+self.Decision_file)

        f = open(self.workpaths+self.horfasta,'w')
        f1 = open(self.workpaths+self.horgff,'w')
        f1.write('\t'.join(['Hor_name','monomer_name','start','end','minimum_distance'])+'\n')

        # 读取序列
        cent = SeqIO.to_dict(SeqIO.parse(self.centri_sequence, "fasta"))
        for seqneme in cent.keys():
            long_sequence = str(cent[seqneme].seq).upper()
            long_sequence_l = len(long_sequence)

            # 设置315到325bp范围的窗口
            window_start = self.seed - self.window
            window_end = self.seed + self.window

            start = 0
            end = self.seed
            start_index = [start]
            length_list = [end]
            while 1:
                min_distance,min_position,order = self.compare_and_plot(seqneme,long_sequence, start, end, window_start, window_end)
                if min_distance == 0 and min_position == 0:
                    break
                elif min_distance == 'None':

                    if len(start_index) >= 2 and start - length_list[-2] == start_index[-2]:# 如果是和后面的序列断开
                        start_index = start_index[:-1]
                        start_index.append(start)
                        start += 1
                        continue
                    elif len(start_index) >= 2 and start - length_list[-2] != start_index[-2]:# 如果是和后面的序列断开
                        start += 1
                        continue

                    start += 1
                    start_index = start_index[:-1]
                    start_index.append(start)
                    continue
                start += end
                end = min_position
                start_index.append(start)
                length_list.append(min_position)
                print(f"Minimum distance: {min_distance}, Position: {min_position}, start: {start}, end: {start+min_position}")
                f.write('>'+seqneme+'_'+str(order)+'\n'+long_sequence[start:start+min_position]+'\n')
                f1.write('\t'.join([seqneme,seqneme+'_'+str(order),str(start),str(start+min_position),str(min_distance)])+'\n')
            # 调用函数绘制矩形块和折线图
            self.plot_rectangles(seqneme, start_index, length_list, long_sequence_l, long_sequence)
        f.close()
        f1.close()

