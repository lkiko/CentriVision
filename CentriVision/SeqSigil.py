from Bio import AlignIO, SeqIO
from io import StringIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import subprocess
import pandas as pd
import logomaker as lm
import matplotlib.pyplot as plt
import os
import math
import seaborn as sns
import CentriVision.bez as bez
from multiprocessing import cpu_count

class SeqSigil():
    def __init__(self, options):
        self.missing_threshold = 0.5
        self.split_position = 75
        self.title_fontsize = 14
        self.axis_fontsize = 10
        self.consensus_sequence = 'consensus_sequence.fasta'
        self.savefig_ic = "ic_profile.png"
        self.savefig_ic_file = "ic_values.tsv"
        self.cpu = cpu_count()

        bez_conf = bez.config()
        for k, v in bez_conf:
            setattr(self, str(k), v)
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)

        self.split_position = int(self.split_position)
        self.title_fontsize = int(self.title_fontsize)
        self.axis_fontsize = int(self.axis_fontsize)
        self.missing_threshold = float(self.missing_threshold)
        self.workpath = os.getcwd() + '/'
        path = bez.get_path()
        font_path = os.path.join(path, 'example/arial.ttf')
        from matplotlib.font_manager import FontProperties
        self.font_prop = FontProperties(fname=font_path)

    def align(self, input_file, output_dir):
        output_file = os.path.join(output_dir, "align_output.fasta")
        if self.align_software == 'mafft':
            command = [self.mafft_path, '--auto', '--adjustdirectionaccurately', '--thread', str(self.cpu), '--clustalout', input_file]
            result = subprocess.run(command, capture_output=True, text=True)
            align = AlignIO.read(StringIO(result.stdout), "clustal")
            AlignIO.write(align, output_file, "clustal")
            return output_file
        elif self.align_software == 'muscle':
            command = [self.muscle_path, '-in', input_file, '-out', output_file, '-clw']
            subprocess.run(command, capture_output=True, text=True)
            return output_file
        elif self.align_software == 'clustalw':
            command = [self.clustalw_path, '-INFILE=' + input_file, '-OUTFILE=' + output_file, '-OUTPUT=CLUSTAL', '-ALIGN']
            subprocess.run(command, capture_output=True, text=True)
            return output_file
        elif self.align_software == 'clustalo':
            command = [self.clustalo_path, '-i', input_file, '-o', output_file, '-t', 'dna', '--threads', str(self.cpu), '--outfmt=clustal']
            subprocess.run(command, capture_output=True, text=True)
            return output_file
        else:
            print('多序列比对软件支持muscle/mafft 请重新选择！')
            exit()

    def run(self):
        align_output = self.align(self.monomer_seq, './')
        alignment = AlignIO.read(align_output, "clustal")

        # 过滤缺失比例过高的列
        filtered_columns = []
        for i in range(alignment.get_alignment_length()):
            column = alignment[:, i]
            missing_ratio = column.count('-') / len(column)
            if missing_ratio < self.missing_threshold:
                filtered_columns.append(column)

        filtered_alignment = MultipleSeqAlignment(
            [SeqRecord(Seq("".join(col)), id=rec.id) for rec, col in zip(alignment, zip(*filtered_columns))]
        )

        with open("filtered_alignment.fasta", "w") as f:
            AlignIO.write(filtered_alignment, f, "clustal")
        alignment = filtered_alignment

        # 构建保守性矩阵
        def build_conservation_matrix(alignment):
            consensus_matrix = {base: [0] * alignment.get_alignment_length() for base in "ACGT-"}
            for record in alignment:
                for i, base in enumerate(record.seq):
                    if base in consensus_matrix:
                        consensus_matrix[base][i] += 1
                    else:
                        consensus_matrix.setdefault(base, [0] * alignment.get_alignment_length())[i] += 1
            df = pd.DataFrame(consensus_matrix).fillna(0)
            df = df.div(df.sum(axis=1), axis=0)
            return df

        conservation_matrix = build_conservation_matrix(alignment)
        color_scheme = {'A': 'red', 'T': 'blue', 'C': 'green', 'G': 'orange'}
        conservation_matrix = conservation_matrix.drop('-', axis=1, errors='ignore')

        # 计算信息量(IC)
        def calculate_ic(df):
            ic_list = []
            for _, row in df.iterrows():
                probs = row.values
                entropy = -sum(p * math.log2(p) if p > 0 else 0 for p in probs)
                ic = math.log2(4) - entropy
                ic_list.append(ic)
            return ic_list

        ic_scores = calculate_ic(conservation_matrix)
        overall_ic = sum(ic_scores) / len(ic_scores)
        print(f"\n整体保守性（平均IC）: {overall_ic:.3f} bits")

        # 区域拆分
        split_len = self.split_position
        total_len = conservation_matrix.shape[0]
        split_indices = list(range(0, total_len, split_len))
        if total_len % split_len < split_len * 0.25 and len(split_indices) > 1:
            split_indices[-1] = total_len
        else:
            split_indices.append(total_len)

        num_rows = len(split_indices) - 1
        fig, axes = plt.subplots(num_rows, 1, figsize=(25, 2 * num_rows))
        if num_rows == 1:
            axes = [axes]

        font_name = self.font_prop.get_name()
        consensus_sequence = []

        for i in range(num_rows):
            start = split_indices[i]
            end = split_indices[i + 1]
            matrix_part = conservation_matrix.iloc[start:end, :]

            # 区域平均IC
            ic_part = ic_scores[start:end]
            avg_ic = sum(ic_part) / len(ic_part)
            print(f"区域 {start}-{end} bp 平均IC: {avg_ic:.3f} bits")

            # 一致性序列
            for _, row in matrix_part.iterrows():
                max_base = row.idxmax()
                consensus_sequence.append(max_base)

            lm.Logo(matrix_part, ax=axes[i], shade_below=0.5, fade_below=0.5, color_scheme=color_scheme)
            axes[i].set_title(f"Alignment Logo ({start}–{end} bp)", fontsize=self.title_fontsize)
            axes[i].tick_params(labelsize=self.axis_fontsize)

        plt.tight_layout()
        plt.savefig(self.savefig, dpi=1000, bbox_inches='tight')
        plt.close()

        # 绘制IC曲线

        # 保存 IC 值到文件
        ic_df = pd.DataFrame({
            "Position(bp)": range(1, len(ic_scores) + 1),
            "IC(bits)": ic_scores
        })
        ic_df.to_csv(self.savefig_ic_file, sep="\t", index=False)

        # 美化 IC 曲线图
        plt.figure(figsize=(12, 3))
        sns.set_style("whitegrid")
        sns.lineplot(x=range(1, len(ic_scores) + 1), y=ic_scores, color='black', lw=1)
        plt.ylim(0, 2)
        plt.xlabel("Position (bp)", fontsize=12)
        plt.ylabel("Information Content (bits)", fontsize=12)
        plt.title("Sequence Conservation Profile", fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig(self.savefig_ic, dpi=1000, bbox_inches='tight')
        plt.close()


        # 输出一致性序列
        consensus_string = ''.join(consensus_sequence)
        revcomp_string = str(Seq(consensus_string).reverse_complement())
        with open(self.consensus_sequence, "w") as f:
            f.write(">Consensus_sequence\n")
            f.write(consensus_string + "\n")
            f.write(">Reverse_complement\n")
            f.write(revcomp_string + "\n")
