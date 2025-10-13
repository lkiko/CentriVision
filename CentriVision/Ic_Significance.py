#!/usr/bin/env python3
# ic_significance.py
# 用途：对多条 monomer 的 per-position IC 文件进行序列级和滑窗显著性检验并可视化
#
# 依赖：pandas, numpy, scipy, statsmodels, matplotlib, seaborn
# pip install pandas numpy scipy statsmodels matplotlib seaborn

import os
import glob
import math
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
import CentriVision.bez as bez

sns.set_style("whitegrid")

class Ic_Significance():
    def __init__(self, options):
        self.pattern = '*.tsv'
        self.outdir = 'ic_results'
        self.min_window = 10
        self.max_window = 10
        self.step = 1
        self.background = 'within'
        self.correction = 'fdr_bh'
        self.wlist = list(range(self.min_window, self.max_window + 1))
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)
        self.min_window = int(self.min_window)
        self.max_window = int(self.max_window)
        self.step = int(self.step)

    # ---------- 工具函数 ----------
    def read_ic_file(self,path: str) -> pd.DataFrame:
        """
        读取单个 IC 文件。自动识别 Position 和 IC 列名（容错）。
        返回 DataFrame，列名为 ['pos','IC']，pos 从 1 起。
        """
        df = pd.read_csv(path, sep=None, engine='python')  # 自动探测分隔符（通常为 tab）
        cols = list(df.columns)
        # 寻找包含Position和IC的列
        pos_col = None
        ic_col = None
        for c in cols:
            cl = c.lower()
            if 'position' in cl or 'pos' == cl or 'bp' in cl:
                pos_col = c
            if 'ic' in cl or 'information' in cl:
                ic_col = c
        # 退化处理：如果只有两列，直接使用
        if pos_col is None or ic_col is None:
            if len(cols) >= 2:
                pos_col = cols[0]
                ic_col = cols[1]
            else:
                raise ValueError(f"无法识别文件 {path} 的列，请确保为两列格式")
        df2 = df[[pos_col, ic_col]].copy()
        df2.columns = ['pos', 'IC']
        # 保证 pos 从 1 开始，IC 为浮点
        df2['pos'] = df2['pos'].astype(int)
        df2['IC'] = pd.to_numeric(df2['IC'], errors='coerce')
        df2 = df2.sort_values('pos').reset_index(drop=True)
        return df2

    def read_ic_dir(self,indir: str, pattern: str = '*.tsv') -> Dict[str, pd.DataFrame]:
        """
        读取目录下所有匹配文件，返回 dict: name -> DataFrame。
        name = 文件名 (不含扩展名)
        """
        files = sorted(glob.glob(os.path.join(indir, pattern)))
        if not files:
            # 也尝试匹配 .txt 或 .csv 或 .tsv
            files = sorted(glob.glob(os.path.join(indir, '*.txt'))) + sorted(glob.glob(os.path.join(indir, '*.csv')))
        ic_dict = {}
        for f in files:
            try:
                df = self.read_ic_file(f)
                name = os.path.splitext(os.path.basename(f))[0]
                ic_dict[name] = df
                print(f"读取 {f} -> {name}, 长度 {len(df)}")
            except Exception as e:
                print(f"跳过 {f}（读取失败）：{e}")
        return ic_dict

    def cohen_d(self,x: np.ndarray, y: np.ndarray) -> float:
        """
        计算 Cohen's d（两组样本的标准化均值差）
        若方差为 0 返回 0
        """
        x = np.asarray(x)
        y = np.asarray(y)
        nx = len(x)
        ny = len(y)
        if nx < 2 or ny < 2:
            return np.nan
        mx = np.mean(x); my = np.mean(y)
        vx = np.var(x, ddof=1); vy = np.var(y, ddof=1)
        pooled_var = ((nx - 1) * vx + (ny - 1) * vy) / (nx + ny - 2) if (nx + ny - 2) > 0 else 0.0
        if pooled_var <= 0:
            return np.nan
        d = (mx - my) / math.sqrt(pooled_var)
        return d

    def significance_label(self,q: float) -> str:
        """
        根据 q 值返回显著性标注
        """
        if pd.isna(q):
            return 'ns'
        if q < 0.001:
            return '***'
        elif q < 0.01:
            return '**'
        elif q < 0.05:
            return '*'
        else:
            return 'ns'

    # ---------- 序列级检验 ----------
    def sequence_level_test(self,ic_dict: Dict[str, pd.DataFrame],
                            correction: str = 'fdr_bh') -> pd.DataFrame:
        """
        对每条序列的 IC 值与其它序列的 IC 值进行 Mann-Whitney U 检验（非参数），
        分别做 one-sided 'greater'（更保守）和 'less'（更不保守）检验。
        对所有序列的 p 值做多重检验校正（分别对 greater 与 less 两组 p 值做校正）。
        返回 DataFrame，包含 mean_IC, median_IC, n, p_greater, q_greater, p_less, q_less, effect sizes, label。
        """
        names = list(ic_dict.keys())
        nseq = len(names)
        results = []
        p_greater = []
        p_less = []
        for name in names:
            arr = ic_dict[name]['IC'].dropna().values
            others = np.concatenate([ic_dict[n]['IC'].dropna().values for n in names if n != name]) if nseq > 1 else np.array([])
            if len(arr) < 1 or len(others) < 1:
                # 不可检验
                p_greater.append(np.nan); p_less.append(np.nan)
                results.append((name, arr.mean() if len(arr) else np.nan, np.median(arr) if len(arr) else np.nan,
                                len(arr), np.nan, np.nan))
                continue
            # one-sided tests
            try:
                stat_g, pgt = mannwhitneyu(arr, others, alternative='greater')
                stat_l, plt_ = mannwhitneyu(arr, others, alternative='less')
            except Exception:
                pgt = np.nan; plt_ = np.nan
            p_greater.append(pgt)
            p_less.append(plt_)
            results.append((name, arr.mean(), np.median(arr), len(arr), pgt, plt_))
        # 多重检验校正（分别对 greater / less）
        p_greater_np = np.array([r[4] for r in results], dtype=float)
        p_less_np = np.array([r[5] for r in results], dtype=float)

        # 校正 greater
        q_greater = np.full_like(p_greater_np, np.nan, dtype=float)
        valid_idx = ~np.isnan(p_greater_np)
        if valid_idx.any():
            reject_g, qvals_g, _, _ = multipletests(p_greater_np[valid_idx], alpha=0.05, method=correction)
            q_greater[valid_idx] = qvals_g

        # 校正 less
        q_less = np.full_like(p_less_np, np.nan, dtype=float)
        valid_idx2 = ~np.isnan(p_less_np)
        if valid_idx2.any():
            reject_l, qvals_l, _, _ = multipletests(p_less_np[valid_idx2], alpha=0.05, method=correction)
            q_less[valid_idx2] = qvals_l

        # 构建最终 DataFrame
        rows = []
        for i, (name, mean_ic, median_ic, n, pgt, plt_) in enumerate(results):
            others_arr = np.concatenate([ic_dict[nm]['IC'].dropna().values for nm in names if nm != name]) if len(names) > 1 else np.array([])
            mdiff = mean_ic - np.mean(others_arr) if len(others_arr) else np.nan
            d = self.cohen_d(ic_dict[name]['IC'].dropna().values, others_arr) if len(others_arr) else np.nan
            qg = q_greater[i] if i < len(q_greater) else np.nan
            ql = q_less[i] if i < len(q_less) else np.nan
            label_g = self.significance_label(qg)
            label_l = self.significance_label(ql)
            # 综合标签：若 greater 显著，标注 'conserved'；若 less 显著，标注 'less_conserved'；二者都不显著为 'ns'
            combined = 'ns'
            if (not pd.isna(qg)) and qg < 0.05:
                combined = 'conserved'
            if (not pd.isna(ql)) and ql < 0.05:
                combined = 'less_conserved'
            # 若两边都显著（极少见），以 q 值小者为准
            if (not pd.isna(qg)) and (not pd.isna(ql)):
                if qg < 0.05 and ql < 0.05:
                    combined = 'conserved' if qg < ql else 'less_conserved'

            rows.append({
                'name': name,
                'mean_IC': mean_ic,
                'median_IC': median_ic,
                'n_pos': n,
                'mean_diff_vs_others': mdiff,
                'cohen_d_vs_others': d,
                'p_greater': pgt,
                'q_greater': qg,
                'sig_greater': label_g,
                'p_less': plt_,
                'q_less': ql,
                'sig_less': label_l,
                'combined_label': combined
            })
        return pd.DataFrame(rows).sort_values('mean_IC', ascending=False)


    # ---------- 滑窗检验 ----------
    def sliding_window_tests(self,ic_dict: Dict[str, pd.DataFrame],
                             window_sizes: List[int] = list(range(10, 21)),
                             step: int = 1,
                             background: str = 'within',  # 'within' 或 'global'
                             correction: str = 'fdr_bh') -> pd.DataFrame:
        """
        对每条序列做滑窗检验。window_sizes 应在 [10,20] 区间内（默认包含端点）。
        background:
          - 'within'：用同一序列的其余位点作为背景
          - 'global'：用整个数据集（所有序列的所有位点，排除当前窗口）作为背景
        返回 DataFrame，包含每个窗口的统计量、p/q 值、显著性等。
        多重检验校正在每一条序列内部（即对该序列的所有窗口 p 值做 FDR 校正），并分别对 greater/less 两类 p 值校正。
        """
        all_seq_names = list(ic_dict.keys())
        # 预先构建全局背景（所有位点）
        global_all_ic = np.concatenate([df['IC'].dropna().values for df in ic_dict.values()]) if len(ic_dict) else np.array([])

        window_rows = []
        for name, df in ic_dict.items():
            arr = df['IC'].values
            L = len(arr)
            # 生成所有窗口（包括不同大小）
            windows = []
            for w in window_sizes:
                if w < 2 or w > L:
                    continue
                for start in range(0, L - w + 1, step):
                    end = start + w  # Python slice [start:end)
                    windows.append((start, end, w))
            if not windows:
                continue

            pvals_g = []
            pvals_l = []
            stats_info = []  # 保存临时信息
            for (start, end, w) in windows:
                window_vals = arr[start:end]
                # background selection
                if background == 'within':
                    bg_vals = np.concatenate([arr[:start], arr[end:]]) if (start > 0 or end < L) else np.array([])
                else:
                    # global background excluding current window values (从 global_all_ic 中减去 window)
                    # 直接排除会比较复杂（若存在相同数值），这里近似使用 global_all_ic 作为背景（不排除当前window）
                    bg_vals = global_all_ic
                # 保证有效样本数
                if len(window_vals) < 2 or len(bg_vals) < 2:
                    p_g = np.nan; p_l = np.nan
                else:
                    try:
                        _, p_g = mannwhitneyu(window_vals, bg_vals, alternative='greater')
                        _, p_l = mannwhitneyu(window_vals, bg_vals, alternative='less')
                    except Exception:
                        p_g = np.nan; p_l = np.nan
                mean_w = np.mean(window_vals) if len(window_vals) else np.nan
                mean_bg = np.mean(bg_vals) if len(bg_vals) else np.nan
                diff = mean_w - mean_bg if (not pd.isna(mean_w) and not pd.isna(mean_bg)) else np.nan
                d = self.cohen_d(window_vals, bg_vals) if (len(window_vals) > 1 and len(bg_vals) > 1) else np.nan
                pvals_g.append(p_g)
                pvals_l.append(p_l)
                stats_info.append((name, start + 1, end, w, mean_w, mean_bg, diff, d, p_g, p_l))  # pos从1开始输出

            # 多重检验校正：对该序列的所有 windows 分别校正 greater 与 less
            pvals_g_np = np.array(pvals_g, dtype=float)
            pvals_l_np = np.array(pvals_l, dtype=float)
            qvals_g = np.full_like(pvals_g_np, np.nan, dtype=float)
            qvals_l = np.full_like(pvals_l_np, np.nan, dtype=float)

            valid_g = ~np.isnan(pvals_g_np)
            valid_l = ~np.isnan(pvals_l_np)
            if valid_g.any():
                _, q_g, _, _ = multipletests(pvals_g_np[valid_g], alpha=0.05, method=correction)
                qvals_g[valid_g] = q_g
            if valid_l.any():
                _, q_l, _, _ = multipletests(pvals_l_np[valid_l], alpha=0.05, method=correction)
                qvals_l[valid_l] = q_l

            # 组织输出行
            for idx, tup in enumerate(stats_info):
                name0, start1, end1, w, mean_w, mean_bg, diff, d, p_g, p_l = tup
                q_g = qvals_g[idx] if idx < len(qvals_g) else np.nan
                q_l = qvals_l[idx] if idx < len(qvals_l) else np.nan
                sig_g = self.significance_label(q_g)
                sig_l = self.significance_label(q_l)
                combined = 'ns'
                if (not pd.isna(q_g)) and q_g < 0.05:
                    combined = 'conserved_window'
                if (not pd.isna(q_l)) and q_l < 0.05:
                    combined = 'less_conserved_window'
                if (not pd.isna(q_g)) and (not pd.isna(q_l)):
                    if q_g < 0.05 and q_l < 0.05:
                        combined = 'conserved_window' if q_g < q_l else 'less_conserved_window'
                window_rows.append({
                    'name': name0,
                    'start_bp': start1,
                    'end_bp': end1,
                    'window_size': w,
                    'mean_window': mean_w,
                    'mean_background': mean_bg,
                    'mean_diff': diff,
                    'cohen_d': d,
                    'p_greater': p_g,
                    'q_greater': q_g,
                    'sig_greater': sig_g,
                    'p_less': p_l,
                    'q_less': q_l,
                    'sig_less': sig_l,
                    'combined_label': combined
                })

        win_df = pd.DataFrame(window_rows)
        # 建议按序列和起始位置排序
        if not win_df.empty:
            win_df = win_df.sort_values(['name', 'start_bp', 'window_size']).reset_index(drop=True)
        return win_df


    def plot_sequence_windows(self,ic_dict: Dict[str, pd.DataFrame],
                              window_df: pd.DataFrame,
                              outdir: str,
                              fontsize_label: int = 15,
                              fontsize_xlabel: int = 18,
                              fontsize_ylabel: int = 18,
                              fontsize_title: int = 20):
        """
        为每条序列绘制 IC 曲线并用颜色条标注显著保守/不保守的窗口。
        保存图片到 outdir/{name}_windows.png
        """
        import os
        import matplotlib.pyplot as plt
        
        os.makedirs(outdir, exist_ok=True)
        for name, df in ic_dict.items():
            L = len(df)
            x = df['pos'].values
            y = df['IC'].values
            fig, ax = plt.subplots(figsize=(12, 2))
            ax.plot(x, y, lw=1, color='black', label='IC')
            ax.set_ylim(0, 2)
            ax.set_xlim(0, L)
            ax.set_xlabel('Position (bp)', fontsize=fontsize_xlabel)
            ax.set_ylabel('IC (bits)', fontsize=fontsize_ylabel)
            ax.set_title(f"{name} IC profile with significant windows", fontsize=fontsize_title)
            ax.tick_params(axis='x', labelsize=fontsize_label)
            ax.tick_params(axis='y', labelsize=fontsize_label)
            # 取该序列所有显著窗口（q < 0.05）
            sub = window_df[window_df['name'] == name]
            if not sub.empty:
                # 绘制保守窗口（sig greater）
                cons = sub[sub['combined_label'] == 'conserved_window']
                for _, r in cons.iterrows():
                    ax.axvspan(r['start_bp'], r['end_bp'] - 1, color='red', alpha=0.25)
                # 绘制低保守窗口（sig less）
                less = sub[sub['combined_label'] == 'less_conserved_window']
                for _, r in less.iterrows():
                    ax.axvspan(r['start_bp'], r['end_bp'] - 1, color='blue', alpha=0.25)
                # 若需要，可在图上标注 stars
                for _, r in sub.iterrows():
                    qg = r['q_greater']; ql = r['q_less']
                    lbl = ''
                    if (not pd.isna(qg)) and qg < 0.05:
                        lbl = self.significance_label(qg)
                    elif (not pd.isna(ql)) and ql < 0.05:
                        lbl = self.significance_label(ql)
                    if lbl != '':
                        ax.text((r['start_bp']+r['end_bp'])/2, 0.5, lbl,
                                ha='center', va='center',
                                fontsize=fontsize_label, rotation=90, fontweight='bold')
            plt.tight_layout()
            outpng = os.path.join(outdir, f"{name}_windows.png")
            plt.savefig(outpng, dpi=1000, bbox_inches='tight')
            plt.close()
            print(f"保存图像 {outpng}")

    # ---------- 主流程（示例） ----------
    def run_pipeline(self,input_dir: str,
                     outdir: str,
                     file_pattern: str = '*.tsv',
                     window_sizes: List[int] = list(range(10, 21)),
                     step: int = 1,
                     background: str = 'within',
                     correction: str = 'fdr_bh'):
        os.makedirs(outdir, exist_ok=True)
        ic_dict = self.read_ic_dir(input_dir, pattern=file_pattern)
        if not ic_dict:
            raise RuntimeError("没有读到任何 IC 文件，请检查目录与文件格式。")

        # 序列级检验
        seq_res = self.sequence_level_test(ic_dict, correction=correction)
        seq_out = os.path.join(outdir, "sequence_level_results.csv")
        seq_res.to_csv(seq_out, index=False)
        print(f"序列级检验结果已保存：{seq_out}")

        # 滑窗检验
        win_res = self.sliding_window_tests(ic_dict, window_sizes=window_sizes, step=step,
                                       background=background, correction=correction)
        win_out = os.path.join(outdir, "window_level_results.csv")
        win_res.to_csv(win_out, index=False)
        print(f"滑窗检验结果已保存：{win_out}")

        # 绘图
        plot_dir = os.path.join(outdir, "plots")
        self.plot_sequence_windows(ic_dict, win_res, plot_dir)
        print("全部完成")

    def run(self):
        self.run_pipeline(self.ic_dir, self.outdir, file_pattern=self.pattern,
                     window_sizes=self.wlist, step=self.step, background=self.background, correction=self.correction)
