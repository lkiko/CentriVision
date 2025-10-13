import os, subprocess, shutil, random
from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import MinMaxScaler
from tqdm import tqdm
import multiprocessing as mp
import CentriVision.bez as bez

class SatAge:
    def __init__(self, options):

        self.blast = 'True'
        self.age = 'True'

        self.blast_dir="blast_tmp"

        self.matched_prefix="matched_regions"

        self.blast_evalue="1e-5"
        self.window_size=10000
        self.k=5
        self.mismatches = 5
        self.distance_threshold=0.6
        self.chrom_label= 20
        self.xlabel= 18
        self.xtick= 16
        self.colorbar_label= 18
        self.legend_fontsize= 18
        self.discrete_colormap='False'
        self.n_bins=8

        bez_conf = bez.config()
        for k, v in bez_conf:# 
            setattr(self, str(k), v)
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)

        self.window_size=int(self.window_size)
        self.k=int(self.k)
        self.mismatches = int(self.mismatches)
        self.distance_threshold=float(self.distance_threshold)
        self.chrom_label= int(self.chrom_label)
        self.xlabel= int(self.xlabel)
        self.xtick= int(self.xtick)
        self.colorbar_label= int(self.colorbar_label)
        self.legend_fontsize= int(self.legend_fontsize)
        self.n_bins=int(self.n_bins)

    def make_blast_db(self):
        os.makedirs(self.blast_dir, exist_ok=True)
        db_prefix = os.path.join(self.blast_dir, 'genome_db')
        subprocess.run([
            self.blast_path+'makeblastdb', '-in', self.genome_fa,
            '-dbtype', 'nucl', '-out', db_prefix
        ], check=True)

    def run_blast(self):
        db_prefix = os.path.join(self.blast_dir, 'genome_db')
        subprocess.run([
            self.blast_path+'blastn', '-query', self.monomers_fa,
            '-db', db_prefix,
            '-out', self.blast_output,
            '-outfmt', '6', '-evalue', self.blast_evalue
        ], check=True)

    def get_chrom_lengths(self, fasta_file):
        return {rec.id: len(rec.seq) for rec in SeqIO.parse(fasta_file, "fasta")}

    def generate_color_map(self, seq_ids):
        random.seed(42)
        return {sid: (random.random(), random.random(), random.random()) for sid in seq_ids}

    def filter_blast_results(self, df, monomer_lens):
        return df[df.apply(lambda r: abs(r['length'] - monomer_lens.get(r['qseqid'], 0)) <= self.mismatches, axis=1)]

    def extract_matched_regions(self, df, chrom_lens):
        genome_dict = SeqIO.to_dict(SeqIO.parse(self.genome_fa, "fasta"))
        records, gff_lines, lens_data = [], [], defaultdict(lambda: {'length': 0, 'count': 0})
        hit_counts = defaultdict(int)

        for _, row in df.iterrows():
            chrom = row['sseqid']
            start, end = min(row['sstart'], row['send']) - 1, max(row['sstart'], row['send'])
            strand = '+' if row['sstart'] < row['send'] else '-'
            hit_counts[chrom] += 1
            monomer_id = row['qseqid']
            ID = f"{chrom}_{monomer_id}_{hit_counts[chrom]}"
            gff_lines.append(f"{chrom}\t{ID}\t{start+1}\t{end}\t{strand}\t{hit_counts[chrom]}\tID={ID}")
            lens_data[chrom]['length'] = len(genome_dict[chrom])
            lens_data[chrom]['count'] += 1
            seq = genome_dict[chrom].seq[start:end]
            if strand == '-':
                seq = seq.reverse_complement()
            records.append(SeqIO.SeqRecord(seq, id=ID, description=""))

        prefix = self.matched_prefix
        with open(f"{prefix}.gff", "w") as f: f.write("\n".join(gff_lines))
        with open(f"{prefix}.lens", "w") as f:
            for chrom, v in lens_data.items():
                f.write(f"{chrom}\t{v['length']}\t{v['count']}\n")
        with open(f"{prefix}.fasta", "w") as f:
            SeqIO.write(records, f, "fasta")

    def plot_blast_hits(self, df, color_map, chrom_lens):
        chroms = list(chrom_lens.keys())
        fig, ax = plt.subplots(figsize=(14, 1 * len(chroms)))  # 更紧凑的布局

        bar_height = 0.6
        spacing = 1.0
        yticks = []
        ylabels = []

        for i, chrom in enumerate(chroms):
            y = i * (bar_height + spacing)
            yticks.append(y + bar_height / 2)
            ylabels.append(chrom)

            chrom_df = df[df['sseqid'] == chrom]
            for _, row in tqdm(chrom_df.iterrows(), total=len(chrom_df), desc=f"{chrom}", leave=False):
                s, e = sorted([row['sstart'], row['send']])
                ax.broken_barh(
                    [(s, e - s)],
                    (y, bar_height),
                    facecolors=color_map.get(row['qseqid'], '#999999'),
                    edgecolors='none'
                )

        # 设置 y 轴为染色体名
        ax.set_yticks(yticks)
        ax.set_yticklabels(ylabels, fontsize=self.chrom_label)
        ax.tick_params(axis='y', length=0)

        # x 轴样式
        ax.set_xlabel("Genomic position (bp)", fontsize=self.xlabel)
        ax.tick_params(axis='x', labelsize=self.xtick)

        # 边框和范围
        ax.set_xlim(0, max(chrom_lens.values()))
        ax.set_ylim(-spacing, yticks[-1] + bar_height + spacing)

        # 去除边框
        for side in ['top', 'right', 'left']:
            ax.spines[side].set_visible(False)

        ax.spines['bottom'].set_linewidth(0.5)

        # 设置图例（颜色-标签映射）
        legend_handles = [mpatches.Patch(color=color, label=label) for label, color in color_map.items()]
        ax.legend(handles=legend_handles, fontsize=self.legend_fontsize, loc='upper right')


        plt.tight_layout()
        plt.savefig(self.blast_hits_file, dpi=1000, bbox_inches='tight')
        print(f"✅ 图像已保存为 {self.blast_hits_file}")


    def compute_age(self):
        k = self.k
        win = self.window_size
        dist_thres = self.distance_threshold

        fasta = SeqIO.to_dict(SeqIO.parse(f"{self.matched_prefix}.fasta", "fasta"))
        gff = pd.read_csv(f"{self.matched_prefix}.gff", sep="\t", header=None,
                          names=["chr", "monomer_id", "start", "end", "strand", "order", "attribute"])

        records = []
        for _, row in gff.iterrows():
            sid = row["monomer_id"]
            if sid in fasta:
                records.append({
                    "chr": row["chr"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "monomer_id": sid,
                    "seq": str(fasta[sid].seq)
                })
        df = pd.DataFrame(records)

        # 创建任务列表：tuple形式，保证可被pickle
        tasks = []
        for chrom in df["chr"].unique():
            chrom_df = df[df["chr"] == chrom].sort_values("start").reset_index(drop=True)
            for i in range(len(chrom_df)):
                tasks.append((chrom, i, chrom_df, k, win, dist_thres))

        with mp.Pool() as pool:
            results = list(tqdm(pool.imap_unordered(SatAge.process_window, tasks), total=len(tasks)))

        final = pd.DataFrame([r for sub in results for r in sub])
        final.to_csv(self.monomer_age_file, sep="\t", index=False)
        print(f"✅ 年龄计算结果已保存到 {self.monomer_age_file}")

    @staticmethod
    def process_window(args):
        chrom, center_idx, sub_df, k, win, dist_thres = args

        def kmer_set(s): return set(s[i:i+k] for i in range(len(s)-k+1))
        def jaccard(s1, s2): return 1 - len(kmer_set(s1) & kmer_set(s2)) / len(kmer_set(s1) | kmer_set(s2))
        def vectorize(seq, kmers): return np.array([1 if k in kmer_set(seq) else 0 for k in kmers])

        center_pos = sub_df.iloc[center_idx]["start"]
        window_df = sub_df[(sub_df["start"] >= center_pos - win) & (sub_df["end"] <= center_pos + win)]
        if len(window_df) < 2:
            return []

        dist_mat = np.zeros((len(window_df), len(window_df)))
        for a in range(len(window_df)):
            for b in range(a + 1, len(window_df)):
                dist = jaccard(window_df.iloc[a]["seq"], window_df.iloc[b]["seq"])
                dist_mat[a, b] = dist_mat[b, a] = dist

        labels = DBSCAN(eps=dist_thres, min_samples=2, metric="precomputed").fit_predict(dist_mat)
        window_df = window_df.reset_index(drop=True)
        window_df["clade_id"] = ["clade" + str(c) if c != -1 else "NA" for c in labels]
        result = []
        for clade in window_df["clade_id"].unique():
            if clade == "NA": continue
            clade_df = window_df[window_df["clade_id"] == clade]
            kmers_union = set().union(*[kmer_set(seq) for seq in clade_df["seq"]])
            vecs = np.array([vectorize(seq, kmers_union) for seq in clade_df["seq"]])
            centroid = vecs.mean(axis=0)
            dists = np.linalg.norm(vecs - centroid, axis=1)
            norm_dists = MinMaxScaler().fit_transform(dists.reshape(-1, 1)).flatten()
            for i, (_, row) in enumerate(clade_df.iterrows()):
                result.append({
                    "chr": row["chr"], "start": row["start"], "end": row["end"],
                    "age": round(norm_dists[i], 4), "clade_id": clade, "monomer_id": row["monomer_id"]
                })
        return result


    def plot_age(self,chrom_lengths):
        df = pd.read_csv(self.monomer_age_file, sep="\t").sort_values(["chr", "start"])
        chroms = df["chr"].unique()
        max_len = max(chrom_lengths.get(c, 0) for c in chroms)

        fig = plt.figure(figsize=(15, 3 * len(chroms) + 1))
        gs = GridSpec(len(chroms)+1, 1, figure=fig, height_ratios=[3]*len(chroms)+[1], hspace=0.5)
        axes = [fig.add_subplot(gs[i, 0]) for i in range(len(chroms))]
        cbar_ax = fig.add_subplot(gs[-1, 0])

        cmap = plt.get_cmap("viridis", self.n_bins) if self.discrete_colormap == 'True' else plt.get_cmap("viridis")
        norm = mcolors.BoundaryNorm(np.linspace(0, 1, self.n_bins + 1), self.n_bins) if self.discrete_colormap == 'True' else plt.Normalize(0, 1)

        for ax, chrom in zip(axes, chroms):
            ax.hlines(0.5, 0, chrom_lengths.get(chrom, max_len), color='black', linewidth=3)


            for _, row in tqdm(df[df["chr"] == chrom].iterrows(), total=len(df[df["chr"] == chrom]), desc=f"{chrom}", leave=False):
                ax.add_patch(patches.Rectangle((row["start"], 0.25), row["end"]-row["start"], 0.5, color=cmap(norm(row["age"]))))

            ax.set_xlim(-1000, max_len + 1000)
            ax.set_yticks([])
            ax.text(-0.01, 0.5, chrom, transform=ax.transAxes, ha='right', va='center', fontsize=self.chrom_label)
            ax.tick_params(axis='x', labelsize=self.xtick)

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.spines["bottom"].set_visible(False)

        axes[-1].set_xlabel("Chromosome Position (bp)", fontsize=self.xlabel)
        

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal', ticks=np.linspace(0, 1, self.n_bins+1))
        cbar.set_label("Monomer Relative Age (0=new, 1=old)", fontsize=self.colorbar_label)
        cbar.ax.tick_params(labelsize=self.xtick)

        plt.savefig(self.age_plot_file, dpi=1000, bbox_inches='tight')
        print(f"✅ 图像已保存为 {self.age_plot_file}")

    def run(self):
        if self.blast == 'True':
            self.make_blast_db()
            self.run_blast()
            shutil.rmtree(self.blast_dir, ignore_errors=True)
        chrom_lens = self.get_chrom_lengths(self.genome_fa)
    
        monomer_lens = self.get_chrom_lengths(self.monomers_fa)
        df = pd.read_csv(self.blast_output, sep='\t', header=None,
                         names=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'])
        df = self.filter_blast_results(df, monomer_lens)
        df = df.sort_values(by=['sseqid', 'sstart'], ascending=[True, True]).reset_index(drop=True)
        color_map = self.generate_color_map(monomer_lens.keys())
        self.plot_blast_hits(df, color_map, chrom_lens)
        self.extract_matched_regions(df, chrom_lens)
        if self.age == 'True':
            self.compute_age()
        self.plot_age(chrom_lens)

