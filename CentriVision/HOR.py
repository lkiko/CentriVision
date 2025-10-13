from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import random

class HOR:
    def __init__(self, parameters=None):
        ### 初始化HOR解析器
        # 默认参数
        self.params = {
            'match_score': 1,
            'mismatch_score': -1,
            'open_gap_score': -2,
            'extend_gap_score': -0.5,
            'base_lengths': [150, 190, 320],
            'length_variation': 10,
            'delta_range': range(-10, 11),
            'min_seed_length': 5
        }
        
        # 合并用户自定义参数
        if parameters:
            self.params.update(parameters)
            
        # 初始化比对器
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'global'
        self.aligner.match_score = self.params['match_score']
        self.aligner.mismatch_score = self.params['mismatch_score']
        self.aligner.open_gap_score = self.params['open_gap_score']
        self.aligner.extend_gap_score = self.params['extend_gap_score']
        
        # 预生成候选长度集合
        self.candidate_Ls = self._generate_candidate_lengths()

    def _generate_candidate_lengths(self):
        """生成候选长度集合"""
        return list({
            L + delta
            for L in self.params['base_lengths']
            for delta in self.params['delta_range']
            if (L + delta) >= self.params['min_seed_length']
        })

    def _compute_similarity(self, s1, s2):
        """计算归一化相似度得分"""
        max_len = max(len(s1), len(s2))
        if max_len == 0:
            return 0.0
        return self.aligner.score(s1, s2) / max_len

    def _reverse_complement(self, s):
        """生成反向互补序列"""
        return str(Seq(s).reverse_complement())

    def split(self, sequence):
        """主分割方法
        
        :param sequence: 输入DNA序列字符串
        :return: tuple (最佳seed长度, 分割结果列表)
                 分割结果列表元素为元组：(segment序列, 方向)
        """
        seq = sequence.upper()
        best_L = None
        best_segments = []
        best_avg_score = -1
        var = self.params['length_variation']

        for L in self.candidate_Ls:
            segments = []
            start = 0
            total_score = 0
            comparisons = 0

            while start < len(seq):
                # 计算可能的结束位置范围
                min_end = start + max(L - var, 1)
                max_end = start + L + var
                end_candidates = range(
                    min_end,
                    min(max_end, len(seq)) + 1
                )

                if not segments:  # 处理首段
                    target = start + L
                    best_end = min(end_candidates, key=lambda x: abs(x - target))
                    segments.append((seq[start:best_end], 'forward'))
                    start = best_end
                    continue

                # 评估后续段
                best_score = -1
                best_dir = 'forward'
                best_end = start + L
                
                for end in end_candidates:
                    current = seq[start:end]
                    prev_seg, _ = segments[-1]
                    
                    # 计算两种方向的相似度
                    fw_score = self._compute_similarity(prev_seg, current)
                    rc_score = self._compute_similarity(
                        prev_seg, 
                        self._reverse_complement(current)
                    )
                    
                    # 更新最佳选择
                    if max(fw_score, rc_score) > best_score:
                        best_score = max(fw_score, rc_score)
                        best_dir = 'forward' if fw_score > rc_score else 'reverse'
                        best_end = end

                total_score += best_score
                comparisons += 1
                segments.append((seq[start:best_end], best_dir))
                start = best_end

            # 评估当前候选长度
            if comparisons > 0:
                avg_score = total_score / comparisons
                if avg_score > best_avg_score:
                    best_avg_score = avg_score
                    best_L = L
                    best_segments = segments

        return best_L, best_segments

# # 使用示例
# if __name__ == "__main__":
#     # 初始化解析器（可自定义参数）
#     hor_parser = HOR({
#         'length_variation': 15,  # 调整长度波动范围
#         'delta_range': range(-15, 16)
#     })
    
#     # 测试序列
#     test_seq = "ATCG"*100 + "GCTA"*50  # 示例序列
    
#     # 执行分割
#     seed_length, segments = hor_parser.split(test_seq)
    
#     # 输出结果
#     print(f"最佳Seed长度: {seed_length}")
#     print(f"分割出{len(segments)}个segment")
#     print("前5个segment示例：")
#     for i, (seg, dir) in enumerate(segments[:5]):
#         print(f"[{i+1}] {dir}方向 | 长度{len(seg)} | 前15bp: {seg[:15]}")