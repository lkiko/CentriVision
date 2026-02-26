# CentriVision

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15910908.svg)](https://doi.org/10.5281/zenodo.15910908)
[![PyPI version](https://badge.fury.io/py/centrivision.svg)](https://badge.fury.io/py/centrivision)

## At present there is no detailed manual for this application, you will simply have to play around and see what happens.

📖 用户手册正在紧锣密鼓地编写中！
I'm working hard on a clear and practical guide—stay tuned!
如有疑问或建议，随时提 Issue，我们一起把它打磨得更好。
Questions or ideas? Open an issue and let’s make it better together.

<!-- 
$${\color{green}Success!}$$  
$${\color{orange}\textbf{Warning!}}$$  
$${\color{red}\textsf{Error!}}$$  
$${\color{blue}Use \space \texttt{CentriVision -trf}}$$  
$${\color{red}\textbf{红色加粗}}$$  
$${\color{orange}\textbf{橙色加粗}}$$  
$${\color{blue}\textbf{蓝色加粗}}$$  
$${\color{green}\textbf{绿色加粗}}$$  
$${\color{purple}\textbf{紫色加粗}}$$  
$${\color{red}\textbf{CentriVision}} \space {\color{blue}\textbf{Configuration}} \space {\color{green}\textbf{Completed!}}$$  
这里是 ${\color{red}\textbf{红色加粗文字}}$ 示例  
这是行内颜色示例： ${\color{blue}blue}$   


$${\color{green}Success!}$$  
$${\color{orange}\textbf{Warning!}}$$  
$${\color{red}\textsf{Error!}}$$  
$${\color{blue}Use \space \texttt{CentriVision -trf}}$$  
$${\color{red}\textbf{红色加粗}}$$  
$${\color{orange}\textbf{橙色加粗}}$$  
$${\color{blue}\textbf{蓝色加粗}}$$  
$${\color{green}\textbf{绿色加粗}}$$  
$${\color{purple}\textbf{紫色加粗}}$$  
$${\color{red}\textbf{CentriVision}} \space {\color{blue}\textbf{Configuration}} \space {\color{green}\textbf{Completed!}}$$  
这里是 ${\color{red}\textbf{红色加粗文字}}$ 示例  
这是行内颜色示例： ${\color{blue}blue}$  
-->

<!-- 目录 -->
- [项目介绍](#项目介绍)
- [安装指南](#安装指南)
  - [环境配置](#环境配置)
- [流程介绍](#流程介绍)
  - [设计理念](#设计理念)
  - [研究流程](#研究流程)
  - [模块说明](#模块说明)
    - [着丝粒鉴定模块](#着丝粒鉴定模块)
    - [宏观结构分析](#宏观结构分析)
    - [微观结构解析](#微观结构解析)
    - [重复单体解析](#重复单体解析)
- [使用说明](#使用说明)
- [使用命令](#使用命令)
- [参数详解](#参数详解)
  - [-trf TRF 串联鉴定](#-trf-TRF-串联鉴定)
  - [-ps PALINDROMIC 回文鉴定](#-ps-PALINDROMIC-回文鉴定)
  - [参数详解](#参数详解)
  - [-hm HEATMAP 着丝粒热图](#-hm-HEATMAP-着丝粒热图)
  - [-e EDISTALN 快速比对](#-e-EDISTALN-快速比对)
  - [-ed EDISTDOT 宏观点图](#-ed-EDISTDOT-宏观点图)
  - [-md DOTPLOT 切片分析](#-md-DOTPLOT-切片分析)
  - [-c COUNT_FILE 统计绘图](#-c-COUNT_FILE-统计绘图)
  - [-m MONOMER 重复单体拆分](#-m-MONOMER-重复单体拆分)
  - [-s SEQSIGIL logo图制作](#-s-SEQSIGIL-logo图制作)
  - [-ic IC_SIGNIFICANCE 位点分歧程度](#-ic-IC_SIGNIFICANCE-位点分歧程度)
  - [-sa SATAGE 分歧推断](#-sa-SATAGE-分歧推断)
  - [-cd COMMUNITY_DETECTION 单元聚类](#-cd-COMMUNITY_DETECTION-单元聚类)
  - [-cm REPEAT_COMMUNITY_MAP 类型分布](#-cm-REPEAT_COMMUNITY_MAP-类型分布)
  - [-gc GET_CENTRI 提取区间](#-gc-GET_CENTRI-提取区间)
  - [-gf GET_CENTGFF 提取注释](#-gf-GET_CENTGFF-提取注释)
  - [-gr GET_REPEAT 提取重复](#-gr-GET_REPEAT-提取重复)
- [文献引用](#文献引用)

## 项目介绍

${\color{orange}\textbf{CentriVision}}$ 是一个用于研究 ${\color{orange}\textbf{着丝粒}}$ 结构的软件工具。
---------------------------------

${\color{orange}\textbf{CentriVision}}$ 旨在提供一个简单而强大的工具，用于分析和可视化基因组中着丝粒的结构。它支持从基因组中提取和分析着丝粒相关信息，并提供丰富的可视化功能，以帮助研究人员深入理解着丝粒的组织和功能。同时支持 ${\color{orange}\textbf{植物和动物}}$ 基因组。
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## 安装指南

<p align="center">
  <img src="https://github.com/lkiko/CentriVision/blob/main/video/install.gif?raw=true" width="100%">
</p>

你可以使用 ${\color{green}\textbf{pip}}$ (https://pypi.org/project/CentriVision/) 来安装 CentriVision：

```bash
pip install CentriVision
```

pip三平台统一版本将在下一个版本更新

${\color{green}\textbf{windows}}$ 用户请使用github/dist/centrivision-1.0.1-py3-none-any.whl版本，使用下面方法安装

或者，你也可以从本地安装 CentriVision 的 wheel 文件：

```bash
pip install CentriVision-x.x.x-py3-none-any.whl
```

软件依赖TRF,Mafft,Muscle,clustalw
使用conda或者mamba配置环境命令

```bash
conda create -n centrivision_env -c bioconda -c conda-forge openmpi trf mafft=7.490 muscle=3.8.1551 clustalw=2.1 clustalo=1.2.4 blast bowtie2=2.4.4

```

激活环境

```bash
conda activate centrivision_env
```

---

### 环境配置：

**以下以 ubuntu操作系统 用户名为charles 的miniconda3 python3.13环境为例**

<p align="center">
  <img src="https://github.com/lkiko/CentriVision/blob/main/video/configuration.gif?raw=true" width="100%">
</p>

查找安装路径：

```bash
pip uninstall CentriVision
```

或者

```bash
python -c "import CentriVision; print(CentriVision.__file__)"
```

默认安装路径：

```bash
/home/charles/miniconda3/lib/python3.13/site-packages/CentriVision/
```

查看依赖文件：

```bash
cat /home/charles/miniconda3/lib/python3.13/site-packages/CentriVision/conf.ini
```

```ini
[ini]
# mpirun_path = mpirun 非必需
mpirun_path = /home/charles/miniconda3/bin/mpirun
trf_path = /usr/bin/trf
# MAFFT v7.490 指定版本
mafft_path = /usr/bin/mafft
# MUSCLE v3.8.1551
muscle_path = /usr/bin/muscle
# CLUSTAL 2.1  指定版本
clustalw_path = /usr/bin/clustalw
# 1.2.4  指定版本
clustalo_path = /usr/bin/clustalo
blast_path = /usr/bin/
# Bowtie 2 version 2.4.4  指定版本
bowtie2_path = /usr/bin/

```

使用vim或其它编辑器修改对应依赖软件位置TRF\Mafft\Muscle 并保存
路径查询

```bash
which trf
which mafft
which bowtie2
which blastn
......
```

ubuntu安装命令

```
sudo apt update
sudo apt install -y trf
sudo apt install -y mafft
sudo apt install -y muscle
sudo apt install -y clustalw
sudo apt install -y clustalo
sudo apt install -y ncbi-blast+
sudo apt install -y bowtie2

```

更新命令：

```
pip install --upgrade CentriVision
```

---

## 流程介绍

<img width="1572" height="1916" alt="image" src="https://github.com/user-attachments/assets/009d1140-ebe2-4ec7-b9d8-ae085c6f2622" />

---

### 设计理念总结

CentriVision 支持：

- 串联重复驱动型着丝粒
- 转座子驱动型着丝粒
- 回文结构驱动型着丝粒

并提供：

- 宏观结构可视化
- 局部结构拆分
- 社区聚类分析
- 单体水平进化研究

实现从 **基因组尺度 → 片段尺度 → 单体尺度** 的多层级解析。
--------------------------------------------------------
### 研究流程

```pgsql
着丝粒鉴定  
↓  
宏观结构分析  
↓  
片段拆分与聚类  
↓  
重复单体拆分  
↓  
单体进化与变异分析 
``` 

--- 

```pgsql
Genome Input
    │
    ▼
┌─────────────────────────────┐
│   Centromere Identification │
└─────────────────────────────┘
    │
    ├── -trf  → TRF-based tandem repeat detection
    │
    ├── -cf   → Custom repeat / TE-based prediction
    │            (supports external GFF input)
    │
    └── -ps   → Palindrome-based detection
    │
    ▼
Centromere Coordinates + Sequences
    │
    ▼
┌─────────────────────────────┐
│   Auxiliary Extraction      │
└─────────────────────────────┘
    │
    ├── -gc  → Extract centromere sequences
    ├── -gr  → Extract repeat sequences
    └── -gf  → Re-index annotations (relative to centromere)
    │
    ▼
Centromere Sequences +
Relative Repeat Annotations
    │
    ▼
┌─────────────────────────────┐
│   Macro-scale Analysis      │
└─────────────────────────────┘
    │
    ├── -e   → Multi-centromere alignment
    ├── -ed  → Dotplot between centromeres
    └── -hm  → Triangular heatmap visualization
    │
    ▼
┌─────────────────────────────┐
│   Micro-scale Analysis      │
└─────────────────────────────┘
    │
    ├── -md  → Split into 4–6 kb fragments
    │           • Dotplot per fragment
    │           • Repeat statistics
    │
    ├── -c   → Clean statistics
    │           • Length distribution
    │
    ├── -cd  → Fragment clustering
    │
    ├── -cm  → Community distribution visualization
    │
    ▼
Extract Continuous Community Regions
    │
    ▼
-m  → Repeat Monomer Decomposition
    │
    ├── -s   → Monomer logo
    ├── -ic  → Mutation analysis
    └── -sa  → Divergence estimation
```

---

<details>
<summary>点击展开查看作者的神神叨叨👏</summary>

---

**CentriVision**提供了**两种**鉴定着丝粒的方法，其一是根据**串联重复序列**来鉴定，使用的模块为-trf/-cf。其中-trf会自动调用本地安装的TRF软件，对基因组进行串联重复序列注释，这一步会输出全基因组的串联重复注释结果gff文件，同时还会生产着丝粒预测结果，但是往往不同的物种着丝粒的大小和重复序列聚集程度不同，这个时候我们可以使用-cf模块，输入前面生成的全基因组的串联重复注释结果gff文件，调整参数之后再预测。

同时-cf模块也提供了更多的可能性。在面对大基因组的时候我们本地跑TRF太慢，我们如果有线程的TRF预测结果gff文件可以直接在-cf模块中使用。

此外，这也为**转座子着丝粒**预测提供了入口，将全基因组的转座子预测结果gff文件输入，就可以预测转座子类型的着丝粒。（注：由于转座子在全基因组中相比于串联重复呈现出的更均匀分布特点，这里预测转座子类型的着丝粒的时候最好输入对应物种着丝粒转座子成分的注释结果。）

另外一种鉴定着丝粒的方法是通过**回文序列**来鉴定，也就是CentriVision中的-ps模块。以上这三个模块在完成鉴定的同时会输出着丝粒在基因组上的位置文件以及对应的序列文件。

---

此外在鉴定这一环节还有一些辅助的模块，如果已经有了着丝粒或者指定区域位置文件想要提取对应的序列，可以使用-gc模块；提取重复序列注释中对应序列的-gr模块；以及提取着丝粒或者指定区域的注释结果，并修改索引为相对于目标区域的-gf模块。

---

经过以上步骤我们确定了基因组中的着丝粒区域，并提取出了对应的完整着丝粒序列以及着丝粒区域的重复序列注释（索引是**相对于着丝粒**）。那么我们就可以进行下一步，从宏观的角度去看着丝粒。

---

CentriVision使用-e模块快速的对所有着丝粒序列进行比对，并使用-ed模块绘制所有着丝粒之间的相似性点阵图。然后使用-hm模块，输入着丝粒序列，以及前面通过-gf模块提取的相对于着丝粒的各类重复序列文件，一般就是TRF、Copia、Gypsy和Gene等，绘制每条着丝粒的三角形热图。

---

进一步进行研究小区域的着丝粒，使用-md模块，将所有着丝粒拆分为4000~6000的片段，并为每一个片段生成点阵图，统计每一个片段内的所有重复序列，会输出统计文件，以及每个片段的位置，序列，点图等文件。这个阶段会生成大量的图片。然后使用-c模块可以读取统计文件，会生成清理好的统计文件以及重复序列长度分布图。

---

将输出的统计文件修改为序列文件，每个-md切片使用统计文件中的种子序列作为代表。输入-cd模块中进行聚类。这一模块最好将过短的部分去除，模块会将相似的序列聚集为一个社区。随后可以使用-cm模块将社区在着丝粒上的分部情况可视化。人工将连续的相同社区分部的区间取出，并使用-gc模块提取出对应的序列进入下一步。

---

将同一社区的序列放在一起使用-m模块拆分出重复单体，这一步每一次拆分都会生成一个分值图，以及一个重复单体，最后会出现较多的图片，以及一个单体合集（fasta文件）。然后使用-s模块既可以获得单体的logo图，使用-ic可以获取单体突变的情况，使用-sa可以估计单体的分歧程度。

---

</details>

<details>
<summary>Click to expand to view the author's ramblings 👏</summary>

---

**CentriVision** provides **two** methods for centromere identification. The first is based on **tandem repeats**, using the `-trf` / `-cf` modules.  
`-trf` automatically invokes a local installation of TRF to annotate tandem repeats across the genome. This step outputs a genome-wide tandem repeat annotation GFF file, along with centromere predictions. However, different species often vary in centromere size and the degree of repeat clustering. In such cases, the `-cf` module can be used, taking the previously generated genome-wide tandem repeat GFF file as input and allowing parameter adjustments for refined predictions.

The `-cf` module also offers more flexibility. When dealing with large genomes, running TRF locally can be too slow. If you already have a precomputed TRF annotation GFF file, you can directly input it into the `-cf` module.

Additionally, this provides an entry point for **transposon-based centromere** prediction. By inputting a genome-wide transposon annotation GFF file, you can predict transposon-derived centromeres.  
(Note: Since transposons are generally more evenly distributed across the genome compared to tandem repeats, it is recommended to input a species-specific annotation of transposon components enriched in centromeres when predicting this type of centromere.)

The second method for centromere identification is based on **palindromic sequences**, implemented in the `-ps` module of CentriVision.  
All three modules output the genomic coordinates of identified centromeres along with their corresponding sequence files.

---

In addition, several auxiliary modules are available for the identification step.  
If you already have centromere coordinates or other region files and wish to extract the corresponding sequences, use the `-gc` module.  
To extract specific sequences from repeat annotations, use `-gr`.  
To extract annotation results for centromeres or specified regions and reindex them relative to the target region, use `-gf`.

---

Once the centromeric regions in the genome have been identified and the corresponding full centromere sequences—as well as repeat annotations indexed **relative to the centromere**—have been extracted, we can proceed to the next step: viewing the centromere from a macro perspective.

---

CentriVision uses the `-e` module to rapidly align all centromere sequences, and the `-ed` module to generate dot plot similarity matrices between them.  
The `-hm` module then takes the centromere sequences along with the various repeat annotation files (typically TRF, Copia, Gypsy, Gene, etc.) extracted via `-gf` and generates a triangular heatmap for each centromere.

---

For further investigation at a finer scale, the `-md` module splits each centromere into 4000–6000 bp fragments, generates a dot plot for each fragment, and tallies all repeat sequences within each fragment. This step outputs summary statistics, along with positional, sequence, and dot plot files for each fragment. This stage produces a large number of images.  
The `-c` module can then read the summary statistics file to generate a cleaned statistical file and a length distribution plot of the repeat sequences.

---

The output statistical file can be reformatted into a sequence file, using the seed sequence from each `-md` fragment as a representative. These are fed into the `-cd` module for clustering. It is recommended to remove overly short fragments before clustering. The module groups similar sequences into communities.  
Subsequently, the `-cm` module visualizes the distribution of these communities along the centromeres. Intervals containing the same community in succession can be manually extracted, and the `-gc` module can be used to retrieve the corresponding sequences for the next step.

---

Sequences from the same community are pooled and processed using the `-m` module to decompose them into repeat monomers. Each run generates a score plot and a repeat monomer, producing many images and a monomer collection (FASTA file).  
The `-s` module then generates a sequence logo of the monomer, `-ic` produces mutation profiles, and `-sa` estimates the divergence level among monomers.

---

</details>

<details>
<summary>クリックして作者の独り言を表示👏</summary>

---

**CentriVision**は、セントロメアを特定する**2つ**の方法を提供します。1つ目は**タンデムリピート**に基づく方法で、`-trf` / `-cf` モジュールを使用します。  
`-trf`はローカルにインストールされたTRFソフトウェアを自動的に呼び出し、ゲノム全体のタンデムリピートアノテーションを実行します。このステップでは、全ゲノムのタンデムリピートアノテーションGFFファイルとともに、セントロメア予測結果も出力されます。しかし、種によってセントロメアのサイズやリピートの集積度は異なるため、そのような場合は`-cf`モジュールを使用し、先ほど生成した全ゲノムタンデムリピートアノテーションGFFファイルを入力として、パラメータを調整した上で再予測することができます。

また、`-cf`モジュールはさらなる柔軟性も提供します。大規模ゲノムを扱う場合、ローカルでTRFを実行するには時間がかかりすぎます。既存のTRF予測結果GFFファイルがあれば、それを直接`-cf`モジュールに入力することが可能です。

さらに、このモジュールは**トランスポゾン由来セントロメア**の予測への入り口にもなります。全ゲノムのトランスポゾン予測結果GFFファイルを入力することで、トランスポゾンタイプのセントロメアを予測できます。  
（注：トランスポゾンはタンデムリピートに比べてゲノム全体により均等に分布する傾向があるため、このタイプのセントロメアを予測する際には、対象種のセントロメアに富むトランスポゾン成分のアノテーション結果を入力することを推奨します。）

2つ目のセントロメア特定方法は**パリンドローム配列**に基づくもので、CentriVisionの`-ps`モジュールで実行します。  
以上の3つのモジュールは、特定したセントロメアのゲノム上の位置情報ファイルと、対応する配列ファイルを出力します。

---

また、特定作業には補助モジュールも用意されています。  
既にセントロメアや目的領域の位置情報ファイルがあり、対応する配列を抽出したい場合は`-gc`モジュールを使用します。  
リピートアノテーションから特定の配列を抽出するには`-gr`モジュールを、セントロメアや目的領域のアノテーション結果を抽出し、インデックスを対象領域に対して相対的に修正するには`-gf`モジュールを使用します。

---

以上のステップを経て、ゲノム中のセントロメア領域を特定し、対応する完全長セントロメア配列、および**セントロメア相対**でインデックスされたセントロメア領域のリピートアノテーションを抽出します。これにより、次のステップである巨視的な視点からのセントロメア解析に進むことができます。

---

CentriVisionの`-e`モジュールは、すべてのセントロメア配列を高速にアラインメントし、`-ed`モジュールはセントロメア間の類似性を示すドットプロットを描画します。  
次に`-hm`モジュールは、セントロメア配列と、先ほど`-gf`モジュールで抽出したセントロメア相対の各種リピートファイル（一般的にはTRF、Copia、Gypsy、Geneなど）を入力として、各セントロメアの三角形ヒートマップを描画します。

---

さらに詳細な局所領域のセントロメア解析を行うには、`-md`モジュールを使用します。このモジュールはすべてのセントロメアを4000～6000 bpの断片に分割し、各断片のドットプロットを作成するとともに、各断片内の全リピート配列を集計します。出力されるのは統計ファイル、各断片の位置情報、配列、ドットプロットファイルなどです。この段階では多数の画像が生成されます。  
次に`-c`モジュールで統計ファイルを読み込むと、クリーニングされた統計ファイルとリピート配列の長さ分布図が生成されます。

---

出力された統計ファイルをシーケンスファイルに変換し、各`-md`断片のシード配列を代表配列として使用します。これを`-cd`モジュールに入力してクラスタリングを実行します。このモジュールでは、短すぎる断片は事前に除去することが推奨されます。類似した配列はコミュニティとしてグループ化されます。  
その後、`-cm`モジュールを使用して、セントロメア上におけるこれらのコミュニティの分布を可視化します。同じコミュニティが連続して分布する区間を手動で抽出し、`-gc`モジュールで対応する配列を取得して次のステップに進みます。

---

同一コミュニティの配列をまとめて`-m`モジュールに入力し、リピートモノマーに分割します。この分割処理のたびにスコアプロットとリピートモノマーが生成され、最終的に多数の画像とモノマーコレクション（FASTAファイル）が得られます。  
その後、`-s`モジュールでモノマーのロゴ図を、`-ic`モジュールでモノマーの変異プロファイルを取得でき、`-sa`モジュールではモノマー間の分岐度を推定することができます。

---

</details>

---

### 模块说明

CentriVision 提供了 **着丝粒鉴定 → 宏观结构分析 → 微观结构分析 → 重复结构解析** 的完整分析流程。整体流程可分为四个阶段：

1. 着丝粒区域鉴定
2. 宏观结构分析
3. 微观结构分析
4. 重复单体解析

---

#### 着丝粒鉴定模块

CentriVision 提供三种核心鉴定策略：

| 方法         | 模块               | 原理                 |
| ------------ | ------------------ | -------------------- |
| 串联重复驱动 | `-trf` / `-cf` | 基于串联重复序列富集 |
| 转座子驱动   | `-cf`            | 基于转座子富集       |
| 回文结构驱动 | `-ps`            | 基于回文序列特征     |

所有鉴定模块均会输出：

- 着丝粒在基因组中的 **位置文件**
- 对应的 **着丝粒序列文件**

---

##### 1. 基于串联重复的鉴定

###### 1.1 `-trf` 模块（自动流程）

该模块会：

1. 自动调用本地安装的 **TRF（Tandem Repeats Finder）**
2. 对全基因组进行串联重复注释
3. 输出：
   - 全基因组串联重复注释 `GFF` 文件
   - 初步的着丝粒预测结果

适用于：

- 中小型基因组
- 串联重复高度聚集的物种

---

###### 1.2 `-cf` 模块（可调节预测）

当不同物种的着丝粒大小或重复聚集程度差异较大时，推荐使用 `-cf`：

输入：

- 全基因组串联重复注释 `GFF` 文件（来自 `-trf` 或外部TRF）

功能：

- 自定义参数重新预测着丝粒
- 支持更灵活的阈值调整

---

###### 1.3 `-cf` 的扩展用途

##### （1）大基因组加速分析

如果本地运行 TRF 过慢：

- 可直接输入已有的 **多线程TRF预测结果GFF文件**
- 无需再次运行TRF

##### （2）转座子型着丝粒预测

`-cf` 同样可用于基于转座子富集的预测：

输入：

- 全基因组转座子注释 GFF 文件

⚠️ 注意：

由于转座子通常在全基因组中分布较为均匀，建议：

- 输入已筛选的 **着丝粒相关转座子成分注释结果**
- 而不是全部转座子注释

---

##### 2. 基于回文序列的鉴定

##### `-ps` 模块

该模块通过识别回文序列结构进行着丝粒预测。

适用于：

- 具有明显回文结构特征的物种
- 串联重复不显著的情况

输出内容同样包括：

- 着丝粒位置文件
- 着丝粒序列文件

---

##### 附、辅助提取模块

在鉴定完成后，CentriVision 提供多个辅助模块用于区域和注释提取。

| 模块    | 功能                                           |
| ------- | ---------------------------------------------- |
| `-gc` | 根据位置文件提取对应序列                       |
| `-gr` | 从重复注释中提取对应序列                       |
| `-gf` | 提取目标区域的注释，并重建为相对于该区域的索引 |

---

##### 典型使用流程

在完成着丝粒鉴定后：

1. 使用 `-gc` 提取完整着丝粒序列
2. 使用 `-gf` 提取着丝粒区域内的重复注释
   （索引转换为相对于着丝粒）

此时我们获得：

- 完整着丝粒序列
- 相对坐标体系下的重复注释

接下来进入结构分析阶段。

---

#### 宏观结构分析

##### 1. 着丝粒之间的整体相似性分析

###### `-e` 模块

- 对所有着丝粒序列进行快速比对

###### `-ed` 模块

- 绘制着丝粒间相似性点阵图

---

##### 2. 单条着丝粒结构可视化

###### `-hm` 模块

输入：

- 着丝粒序列
- 相对于着丝粒的重复注释文件（来自 `-gf`）

常见注释类型：

- TRF（串联重复）
- Copia
- Gypsy
- Gene

输出：

- 每条着丝粒的三角形热图

用于展示重复结构的空间分布特征。

---

#### 微观结构解析（片段级分析）

用于研究局部重复单元组织结构。

---

##### 1. 片段拆分与统计

###### `-md` 模块

功能：

- 将所有着丝粒拆分为 4000–6000 bp 片段
- 为每个片段：
  - 生成点阵图
  - 统计重复序列
  - 输出统计文件
  - 输出片段位置与序列文件

⚠️ 该阶段会生成大量图片文件。

---

##### 2. 统计结果整理

###### `-c` 模块

输入：

- `-md` 生成的统计文件

输出：

- 清理后的统计文件
- 重复序列长度分布图

---

##### 3. 片段聚类分析

步骤：

1. 将统计文件转换为序列文件
2. 每个片段以统计文件中的“种子序列”作为代表

###### `-cd` 模块

功能：

- 对片段进行聚类
- 将相似序列划分为不同社区

建议：

- 过滤过短序列后再进行聚类

---

###### `-cm` 模块

- 可视化不同社区在着丝粒上的分布情况

随后：

- 人工提取连续的相同社区区间
- 使用 `-gc` 提取对应序列

进入重复单体解析阶段。

---

#### 重复单体解析

##### 1. 单体拆分

###### `-m` 模块

功能：

- 将同一社区序列拆分为重复单体

输出：

- 每轮拆分的分值图
- 单个重复单体文件
- 单体合集 FASTA 文件

该步骤会产生较多图像文件。

---

##### 2. 单体特征分析

| 模块    | 功能             |
| ------- | ---------------- |
| `-s`  | 生成单体 Logo 图 |
| `-ic` | 分析单体突变情况 |
| `-sa` | 估计单体分歧程度 |



## 使用命令

```bash
CentriVision -h
```

```bash
CentriVision options ?/xx.conf
```

运行命令参数解释
${\color{red}\textbf{?}}$ 表示询问模块参数
${\color{red}\textbf{xx.conf}}$ 配置文件内是模块需要的参数
${\color{red}\textbf{? > xx.conf}}$ 询问模块参数并将输出的内容覆盖式输入到后续的xx.conf配置文件中
${\color{red}\textbf{? >> xx.conf}}$ 询问模块参数并将输出的内容追加输入到后续的xx.conf配置文件中

```
usage: CentriVision [options]
runing CentriVision
options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -ps PALINDROMIC       Palindromic sequence 查询基因组中的回文序列
  -trf TRF              run TRF(Tandem Repeat Finder) 通过TRF查找串联重复序列；
  -cf CENTRIFINDER      Centrifinder 着丝粒预测；
  -md DOTPLOT           mini Dotplot 重复序列点图；
  -hm HEATMAP           Heatmap 区域相似度热图；
  -m MONOMER            Monomer scanning 重复单体扫描；
  -s SEQSIGIL           SeqSigil scanning 重复单体logo；
  -ic IC_SIGNIFICANCE   Ic Significance 单体保守性IC检验；
  -sa SATAGE            SatAge Monomer 重复时间推断（拟分子钟）；
  -gc GET_CENTRI        Get_centri 提取基因组的指定区域；
  -gf GET_CENTGFF       Get_centgff 提取基因组的指定区域gff,index修改为相对着丝粒；
  -gr GET_REPEAT        Get_repeat 根据gff3提取基因组的重复序列；
  -c COUNT_FILE         Count_file 统计dotplot文件；
  -r HOR                HOR HOR搜索；
  -ed EDISTDOT          EdistDot EdistDot 点阵图；
  -e EDISTALN           EdistAln EdistAln 快速比对；
  -cd COMMUNITY_DETECTION
                        Community_detection 重复序列社区发现；
  -cm REPEAT_COMMUNITY_MAP
                        Repeat_community_map 重复序列社区映射；

```

---

## 参数详解

### -trf TRF 串联鉴定

调用TRF(Tandem Repeat Finder)扫描重复序列，根据重复序列判断着丝粒。

<p align="center">
  <img src="https://github.com/lkiko/CentriVision/blob/main/video/trf.gif?raw=true" width="100%">
</p>

查看参数：

```bash
CentriVision -trf ?
```

![参数](https://github.com/user-attachments/assets/a4480953-31ee-461d-a35e-cd71fd3e5dbc)

参数重定向到配置文件total.conf

覆盖式命令：

```bash
CentriVision -trf ? > total.conf
```

追加式命令：

```bash
CentriVision -trf ? >> total.conf
```

![参数](https://github.com/user-attachments/assets/7be1829e-65e6-4b83-b05a-6bca0711766c)

配置文件：

```
[TRF]
genome_file = genome file
lens = lens file
chip_seq = chip_seq map file or None
colors = hish,centri,chip or hish,centri,None or #38b48b,#1e50a2,#d7003a
trfgff = out gff
trffasta = out fasta
windows = 10000
step = 5000
gap = 40
centrigff = centri gff
centrifasta = centri fasta
```

genome_file = genome file 基因组fasta文件
lens = lens file 染色体文件
chip_seq = chip_seq map file or None ChIP-seq或其它数据的先验着丝粒位置文件
colors = hish,centri,chip or hish,centri,None or #38b48b,#1e50a2,#d7003a 颜色配置
trfgff = out gff TRF输出gff3结果
trffasta = out fasta TRF输出fasta文件
windows = 10000 重复序列密度窗口跨度 根据基因组大小和重复序列密度自行调整
step = 5000 重复序列密度窗口滑动步长 根据基因组大小和重复序列密度自行调整
gap = 40 重复区域连续性容错宽度 gap\*windows 根据基因组大小和重复序列密度自行调整
centrigff = centri gff 鉴定候选区结果
centrifasta = centri fasta 候选区fasta文件

lens文件：

```
#染色体号 染色体长度 基因数
osa2  36447916  6689
```

chip_seq文件：

```
#染色体号 起始  终止
osa2  13619000  14176000
```

trfgff文件：

```
#染色体号 TRF TandemRepeat  起始  终止  - - - 注释
osa2  TRF TandemRepeat  1 6615        ID=TRF00001;PeriodSize=7;CopyNumber=944.9;Consensus=CCCTAAA
osa2  TRF TandemRepeat  9594  9661        ID=TRF00002;PeriodSize=34;CopyNumber=2.0;Consensus=CTCCAAAACCATGGAGGAAGTCAAATTACACCGA
osa2  TRF TandemRepeat  19795 19826       ID=TRF00003;PeriodSize=3;CopyNumber=10.7;Consensus=CGG
osa2  TRF TandemRepeat  20033 20061       ID=TRF00004;PeriodSize=6;CopyNumber=5.0;Consensus=GGGGCG
osa2  TRF TandemRepeat  22330 22433       ID=TRF00005;PeriodSize=27;CopyNumber=3.9;Consensus=TTCCCAGGAGGGATGCCTGGTGGAGGT
osa2  TRF TandemRepeat  22326 22457       ID=TRF00006;PeriodSize=54;CopyNumber=2.5;Consensus=GGGCTTCCCAGGTGCTATGCCTGGTGGAGGTTTCCCAGGAGGAATGCCTGGT
osa2  TRF TandemRepeat  25152 25199       ID=TRF00007;PeriodSize=24;CopyNumber=2.0;Consensus=GCATGCAAAGCAAGTAATAATAGG
osa2  TRF TandemRepeat  25446 25500       ID=TRF00008;PeriodSize=6;CopyNumber=9.2;Consensus=ATATAG
......
```

修改配置文件如下：

```
[TRF]
genome_file = NIP-T2T-osa2.fa
lens = osa.lens
chip_seq = None
colors = #38b48b,#1e50a2,#d7003a
trfgff = out.gff
trffasta = out.fasta
windows = 10000
step = 5000
gap = 40
centrigff = centri.gff
centrifasta = centri.fasta
```

功能执行
命令：

```bash
CentriVision -trf total.conf
```

<p align="center">
  <img src="https://github.com/lkiko/CentriVision/blob/main/video/trf-run.gif?raw=true" width="100%">
</p>
输出结果:

![图示](https://github.com/user-attachments/assets/3e2e19d7-7256-4cd9-9495-047d7f29a4be)
淡蓝色为重复序列分布图，蓝色区域为着丝粒候选区域

![结果](https://github.com/user-attachments/assets/83d37692-95db-4380-b960-89230164c18a)

$$
{\color{red}\textbf{TRF在面对大区域重复的时候扫描特别慢，可以单独切片运行TRF注释}}$$  

### -cf CENTRIFINDER 快速鉴定

模块，输入文件兼容 ${\color{orange}\textbf{串联重复注释文件、转座子注释文件或者回文序列注释}}$  

覆盖式命令：  
```bash
CentriVision -cf ? > total.conf
```

追加式命令：  
```bash
CentriVision -cf ? >> total.conf
```
配置文件：  
```
[Centrifinder]
genome_file = genome file
lens = lens file
chip_seq = chip_seq map file or None
colors = hish,centri,chip or hish,centri,None or #38b48b,#1e50a2,#d7003a
trfgff = out gff
windows = 10000
step = 5000
gap = 40
centrigff = centri gff
centrifasta = centri fasta
```
配置文件和trf模块类似  
genome_file = genome file 基因组fasta文件  
lens = lens file 染色体文件  
chip_seq = chip_seq map file or None ChIP-seq或其它数据的先验着丝粒位置文件  
colors = hish,centri,chip or hish,centri,None or #38b48b,#1e50a2,#d7003a 颜色配置  
trfgff = out gff ${\color{orange}\textbf{串联重复注释文件、转座子注释文件或者回文序列注释文件}}$  
trffasta = out fasta TRF输出fasta文件  
windows = 10000 重复序列密度窗口跨度  
step = 5000 重复序列密度窗口滑动步长  
gap = 40 重复区域连续性容错宽度 gap\*windows  
centrigff = centri gff 鉴定候选区结果  
centrifasta = centri fasta 候选区fasta文件  

运行方式同上
命令：  
```bash
CentriVision -cf total.conf
```

### -ps PALINDROMIC 回文鉴定

PALINDROMIC 模块，通过染色体 ${\color{orange}\textbf{回文序列}}$ 密度来鉴定着丝粒

覆盖式命令：  
```bash
CentriVision -ps ? > total.conf
```

追加式命令：  
```bash
CentriVision -ps ? >> total.conf
```
配置文件：  
```
[Palindromic]
genome_file = genome file
length = 10
reach = 2000
windows = 10000
step = 5000
coln = 3
width = 15
height = 10
gff_file = Palindromic gff
savefile = save file (*.png, *.pdf, *.svg)
```
配置文件  
genome_file = genome file 基因组fasta文件  
length = 10 回文序列长度  
reach = 2000 回文最大距离  
windows = 10000 密度窗口  
step = 5000 密度窗口滑动距离  
coln = 3 绘图列数  
width = 15 绘图宽度  
height = 10 绘图高度  
gff_file = Palindromic gff 输出gff文件  
savefile = save file (\*.png, \*.pdf, \*.svg) 可视化输出  

运行方式同上
命令：  
```bash
CentriVision -ps total.conf
```
![ps](https://github.com/user-attachments/assets/d8917a8f-a888-49d7-ab8a-17f99b71ee47)

增加ChIP-seq先验数据(红色区域)  

![ps](https://github.com/user-attachments/assets/fbe6a8f5-a0f1-421b-8318-64f9b109fe1c)


---

### -hm HEATMAP 着丝粒热图

着丝粒内差异化热图  

查看参数：  
```bash
CentriVision -hm ?
```

![md参数](https://github.com/user-attachments/assets/47e4323a-eef5-44c9-817d-bfb597a676eb)


参数重定向到配置文件total.conf

覆盖式命令：  
```bash
CentriVision -hm ? > total.conf
```

追加式命令：  
```bash
CentriVision -hm ? >> total.conf
```

配置文件：  
```
# 图片生成可以中断再运行
[Heatmap]
centromere_file = genome file
align_software = ialign or muscle or mafft or hamming
# 比对模式，global考虑全局相似性，全部序列比对；local序列两两比对; ialign c模块快速比对
model = global or local
# 负值表示反向互补/反向
reverse_complement = True/False
# 软件支持离散着色和渐变着色 Discrete or Gradient
color_mode = Discrete or Gradient
# 选择是否需要绘制注释信息，GC,TRF,TE
annotation = True/False
trf_gff = None or gff:#e2041b
te_gff = None or gff:#19448e
gene_gff = None or gff:#b44c97
chip_file = None or txt:#3eb370
gcc = #0d0015
# 确保着丝粒拆分的不能太短，协调计算机的内存、算力和着丝粒长度/global模式下，是最短着丝粒的切割份数，local模式下，是每条着丝粒单独拆分
split = 1000
# segment_length为0时split生效，若segment_length不为零为拆分长度
segment_length = 0
out_path = out path
```
参数详解 ${\color{orange}\textbf{可中断接续运行}}$ ：  
centromere_file = genome file 着丝粒文件  
align_software = ialign/muscle/mafft/hamming 比对软件，可以调用现有软件，也可以使用 ${\color{orange}\textbf{ialign}}$ c模块快速比对  
model = global/local 比对模式，global考虑全局相似性，全部序列同时比对；local序列两两比对  
reverse_complement = True/False 负值表示 ${\color{orange}\textbf{反向互补/反向}}$  
color_mode = Discrete/Gradient 软件支持 ${\color{orange}\textbf{离散着色/渐变着色}}$  Discrete or Gradient  
annotation = True/False 选择是否需要绘制注释信息，GC,TRF,TE  
trf_gff = None or gff:#e2041b ${\color{orange}\textbf{None}}$ 表示无注释，若有注释则使用 ${\color{orange}\textbf{:}}$ 分割文件名和展示颜色  
te_gff = None or gff:#19448e  
gene_gff = None or gff:#b44c97  
chip_file = None or txt:#3eb370  
gcc = #0d0015 GC含量展示颜色  
split = 1000 着丝粒拆分的不能太短，协调计算机的内存、算力和着丝粒长度/global模式下，是最短着丝粒的切割份数，local模式下，是每条着丝粒单独拆分  
segment_length = 0 segment_length为0时split生效，若segment_length不为零为拆分长度  
out_path = out path 输出路径  

修改配置文件如下：  
```
[Heatmap]
centromere_file = centri.fasta
align_software = ialign
model = local
reverse_complement = True
color_mode = Gradient
annotation = False
trf_gff = None
te_gff = None
gene_gff = None
chip_file = None
gcc = #0d0015
split = 1000
segment_length = 0
out_path = hmap

```

![md参数修改](https://github.com/user-attachments/assets/121aaaa2-f386-46f6-a6c7-6341faf6fde0)

功能执行  
命令：  
```bash
CentriVision -hm total.conf
```
![hm-run](https://github.com/user-attachments/assets/30b08149-a273-47fe-a866-48ce43e3fbae)

ialign比对结果
![ialign](https://github.com/user-attachments/assets/abd6f57f-8c5f-460e-b713-4b5984969402)

重复序列相似性热图：
![s02_1](https://github.com/lkiko/CentriVision/assets/57740432/d0b95ae5-d83f-4997-9410-2768ddc296bf)

${\color{red}\textbf{切片大小}}$ 与分辨率和计算机内存大小挂钩  

#### 测试样本 osa2着丝粒   ${\color{orange}\textbf{585000bp}}$  
拆分为 ${\color{orange}\textbf{1000}}$ 份，需要将1000份子序列进行两两比对，平台 ${\color{orange}\textbf{Intel(R) Core(TM) i7-8565U CPU @ 1.80GHz 8进程}}$ 运算  
${\color{orange}\textbf{汉明距离}}$ （速度最快，但是准确性太差，主要问题是对于低相似性序列之间效果太差）  
![汉明距离时间](https://github.com/user-attachments/assets/cd55c48d-3b6f-46a8-ad52-b635ef6504a2)  
运行时长 ${\color{orange}\textbf{118}}$ 秒  
![汉明距离](https://github.com/user-attachments/assets/6582650e-196e-4575-85e6-560c17ca2d6c)  

${\color{orange}\textbf{C语言脚本}}$ （速度快，准确率也不错） ${\color{red}\textbf{推荐使用ialign}}$  
![ialign时间](https://github.com/user-attachments/assets/f505a2e9-c1f1-4a13-a441-01fefa0af36c)  
运行时长 ${\color{orange}\textbf{2626}}$ 秒（ ${\color{orange}\textbf{43}}$ 分钟）  
![ialign结果](https://github.com/user-attachments/assets/5c04c2d1-0b2e-47ba-8ae2-a246d37d0815)  

调用 ${\color{orange}\textbf{mafft}}$ （速度慢）  
![mafft时间](https://github.com/user-attachments/assets/962d4dcf-a978-4ac7-853a-259e26fee501)  
运行时长 ${\color{orange}\textbf{69619}}$ 秒（ ${\color{orange}\textbf{19}}$ 小时）  
![mafft](https://github.com/user-attachments/assets/275b1e3d-c4ff-432c-ab18-2a3e487d4f08)  

调用 ${\color{orange}\textbf{muscle}}$ （速度比mafft快，比ialign慢，对精确度没有ialign和mafft好） 
![muscle时间](https://github.com/user-attachments/assets/171561ba-a06d-4d28-9529-3a59b970903e)  
运行时长 ${\color{orange}\textbf{8867}}$ 秒（ ${\color{orange}\textbf{147}}$ 分钟） 
![muscle](https://github.com/user-attachments/assets/daee0171-ee6b-4dc4-9d10-c38ab85ead1c)  


---

### -e EDISTALN 快速比对

着丝粒并快速比对。  

查看参数：  
```bash
CentriVision -e ?
```

参数重定向到配置文件total.conf

覆盖式命令：  
```bash
CentriVision -e ? > total.conf
```

追加式命令：  
```bash
CentriVision -e ? >> total.conf
```

配置文件：  
```
[EdistAln]
centri_sequence = centri file
window = 2000
cpu = 8
out_file = out file (\*.tsv)
```
参数详解：  
centri_sequence = centri file 着丝粒fasta文件  
window = 2000 切片比对宽度  
cpu = 8 多进程  
out_file = out file (\*.tsv) 比对结果输出  

功能执行  
命令：  
```bash
CentriVision -e total.conf
```
---


### -ed EDISTDOT 宏观点图

基因组着丝粒比对结果绘制。  

![ed](https://github.com/user-attachments/assets/aa34eb40-d1f5-47bf-8439-ffcc208ad7ed)  

查看参数：  
```bash
CentriVision -ed ?
```

参数重定向到配置文件total.conf

覆盖式命令：  
```bash
CentriVision -ed ? > total.conf
```

追加式命令：  
```bash
CentriVision -ed ? >> total.conf
```

配置文件：  
```
[EdistDot]
genepairs = colinearity file
genepairsfile_type = EdistAln/BLAST/MCScanX/ColinearScan
gff1 =  gff1 file
gff2 =  gff2 file
lens1 = lens1 file
lens2 = lens2 file
genome1_name =  Genome1 name
genome2_name =  Genome2 name
position = order
blast_reverse = false
block = 0
markersize = 0.5
figsize = 10,10
savefig = savefile(.png, .pdf, .svg)

# 其他参数
genome_name_size = 30
chr_name_size = 20
tandem = True
levels = 1:1:0
q_s = 1:1
```
参数详解：  
genepairs = colinearity file 比对结果文件  
genepairsfile_type = EdistAln/BLAST/MCScanX/ColinearScan 比对结果格式  
gff1 =  gff1 file gff文件  
gff2 =  gff2 file gff文件  
lens1 = lens1 file lens文件  
lens2 = lens2 file lens文件  
genome1_name =  Genome1 name  
genome2_name =  Genome2 name  
position = order/end order使用相对位置，end使用绝对位置  
blast_reverse = false 是否需要交换顺序  
block = 0 最断共线性  
markersize = 0.5 点大小  
figsize = 10,10 图片比例  
savefig = savefile(.png, .pdf, .svg) 保存格式  
genome_name_size = 30 基因组名字体大小  
chr_name_size = 20 染色体名字体大小  
tandem = True 是否去除串联重复序列  
levels = 1:1:0 blast结果显示比例  
q_s = 1:1  

简化gff文件：  
```
#染色体号 切片id  起始  终止  - - - 注释
chr01_1 chr01_1_s0  0 6000
chr01_1 chr01_1_s1  6000  12000
chr01_1 chr01_1_s2  12000 18000
chr01_1 chr01_1_s3  18000 24000
chr01_1 chr01_1_s4  24000 30000
chr01_1 chr01_1_s5  30000 36000
chr01_1 chr01_1_s6  36000 42000
chr01_1 chr01_1_s7  42000 48000
chr01_1 chr01_1_s8  48000 54000
chr01_1 chr01_1_s9  54000 60000
......
```

功能执行  
命令：  
```bash
CentriVision -ed total.conf
```
---

### -md DOTPLOT 切片分析

拆分着丝粒并扫描重复序列。  

查看参数：  
```bash
CentriVision -md ?
```

![md参数](https://github.com/user-attachments/assets/47e4323a-eef5-44c9-817d-bfb597a676eb)


参数重定向到配置文件total.conf

覆盖式命令：  
```bash
CentriVision -md ? > total.conf
```

追加式命令：  
```bash
CentriVision -md ? >> total.conf
```

配置文件：  
```
[Dotplot]
# 窗口宽度根据内存大设定默认4000
genome_file = genome file
minlength = 10
windows = 4000
poly = False
plot = False
temp = False
cpu = 16
outfile = out dotplot
```
参数详解：  
genome_file = genome file 着丝粒fasta文件  
minlength = 10 最小重复单元  
windows = 4000 切片大小，根据不同物种的重复单元大小，计算机内存大小等合理设置，一般保持 ${\color{orange}\textbf{20个重复单元}}$ 左右最为清晰  
poly = False 去除序列中的单碱基重复区域，默认不开启  
plot = False 输出单独的自相似矩阵图  
temp = False 是否保留比对矩阵  
cpu = 16 多进程  
outfile = out dotplot 输出文件，包含每个切片的特征矩阵(tab隔开)  

修改配置文件如下：  
```
[Dotplot]
# 窗口宽度根据内存大设定默认4000
genome_file = centri.fasta
minlength = 8
windows = 4000
poly = False
plot = True
temp = False
cpu = 8
outfile = out.dotplot
```

![md参数修改](https://github.com/user-attachments/assets/121aaaa2-f386-46f6-a6c7-6341faf6fde0)

功能执行  
命令：  
```bash
CentriVision -md total.conf
```

结果：  
点阵图：展示序列重复规律的点图  
![md](https://github.com/user-attachments/assets/af6059fe-f27f-49c2-8fe7-8c90d39b4207) 

查找重复单元的相位纠正图  
![md](https://github.com/user-attachments/assets/58e31459-35d8-4d9d-aace-2dbc4bba9d55)  

${\color{red}\textbf{切片大小}}$ 与分辨率和计算机内存大小挂钩，大型矩阵极其消耗内存；对于具有 ${\color{red}\textbf{超大着丝粒}}$ 的物种，切片数量非常多，是否需要输出所有自相似矩阵图以及比对矩阵需要适当选择，可利用输出文件可选的生成对应切片的自相似矩阵图和比对矩阵

### -c COUNT_FILE 统计绘图

统计着丝粒重复单元信息。  

查看参数：  
```bash
CentriVision -c ?
```

参数重定向到配置文件total.conf

覆盖式命令：  
```bash
CentriVision -c ? > total.conf
```

追加式命令：  
```bash
CentriVision -c ? >> total.conf
```

配置文件：  
```
[Count_file]
dot_file = dotplot file
# 计入统计的重复单元长度范围
lmmin = 8
lmmax = 1000
# 分箱宽度
bin_size = 10
# 关注前几个柱体
peak_index = 1
# 关注其它柱体
peak_indices = None or Other bars example: 1,2,3
# y_break_min不为0时绘制断轴图，设置省略范围
y_break_min = 0
y_break_max = 0
out_file = new dotplot file
savefile = save file (\*.png, \*.pdf, \*.svg)
```
参数详解：  
dot_file = dotplot file md模块输出文件  
lmmin = 8 计入统计的最小重复单元长度   
lmmax = 1000 计入统计的最大重复单元长度  
bin_size = 10 分箱宽度  
peak_index = 1 关注前几个柱体   
peak_indices = None or Other bars example: 1,2,3 关注其它柱体  
y_break_min = 0 y_break_min不为0时绘制断轴图，设置省略范围下限  
y_break_max = 0 设置省略范围上限  
out_file = new dotplot file 过滤md模块输出文件  
savefile = save file (\*.png, \*.pdf, \*.svg) 可视化输出  

功能执行  
命令：  
```bash
CentriVision -c total.conf
```
---


### -m MONOMER 重复单体拆分

将较为均匀的重复区域拆分为单体。  

![md](https://github.com/user-attachments/assets/226ab6d2-9188-4fd7-88d6-df5a66cd9116)  

查看参数：  
```bash
CentriVision -m ?
```

参数重定向到配置文件total.conf

覆盖式命令：  
```bash
CentriVision -m ? > total.conf
```

追加式命令：  
```bash
CentriVision -m ? >> total.conf
```

配置文件：  
```
# 单体扫描
[Monomer]
centri_sequence = centri file
seed = 320
window = 20
```
参数详解：  
centri_sequence = centri file 着丝粒fasta文件  
seed = 320 提示重复单元长度  
window = 20 允许的差异范围  

功能执行  
命令：  
```bash
CentriVision -m total.conf
```
---

### -s SEQSIGIL logo图制作

拆分着丝粒并扫描重复序列。  

![md](https://github.com/user-attachments/assets/2ca1ccab-b426-4329-bb6b-22f85ede5c5f)  
![md](https://github.com/user-attachments/assets/eccfd242-0337-4e32-a310-47f2a00a711f)  

查看参数：  
```bash
CentriVision -s ?
```

参数重定向到配置文件total.conf

覆盖式命令：  
```bash
CentriVision -s ? > total.conf
```

追加式命令：  
```bash
CentriVision -s ? >> total.conf
```

配置文件：  
```
[SeqSigil]
monomer_seq = monomer file
align_software = muscle or mafft or clustalw or clustalo
missing_threshold = 0.5
split_position = 150
savefig = savefile(.png, .pdf, .svg)

# 可调参数：标题和坐标轴字体大小
title_fontsize = 20
axis_fontsize = 18
```
参数详解：  
monomer_seq = monomer file m模块拆分的重复单体fasta文件  
align_software = muscle or mafft or clustalw or clustalo 比对方法  
missing_threshold = 0.5 最低比对  
split_position = 150 最大宽度  
savefig = savefile(.png, .pdf, .svg) 可视化输出  
title_fontsize = 20 标题字体大小  
axis_fontsize = 18 坐标轴字体大小  

功能执行  
命令：  
```bash
CentriVision -s total.conf
```
---

### -ic IC_SIGNIFICANCE 位点分歧程度

扫描重复单体不同区域的保守性。  

查看参数：  
```bash
CentriVision -ic ?
```

参数重定向到配置文件total.conf

覆盖式命令：  
```bash
CentriVision -ic ? > total.conf
```

追加式命令：  
```bash
CentriVision -ic ? >> total.conf
```

配置文件：  
```
[Ic_Significance]
ic_dir = IC file idr
pattern = *.tsv
min_window = 10
max_window = 10
step = 1
# 背景选择：within 使用同序列其余位点；global 使用全体位点作为背景（默认 within）
background = within
# 多重检验校正方法（statsmodels 支持的方法），默认 fdr_bh
# bonferroni：Bonferroni 校正
# sidak：Sidak 校正
# holm-sidak：Holm-Sidak 校正
# holm：Holm 校正
# simes-hochberg：Simes-Hochberg 校正
# hommel：Hommel 校正
# fdr_bh：Benjamini-Hochberg FDR 校正（默认）
# fdr_by：Benjamini-Yekutieli FDR 校正
# fdr_tsbh：Two-stage Benjamini-Hochberg FDR 校正
# fdr_tsbky：Two-stage Benjamini-Krieger-Yekutieli FDR 校正
correction = fdr_bh
```
参数详解：  
ic_dir = IC file idr  
pattern = \*.tsv  
min_window = 10  
max_window = 10  
step = 1  
background = within/global within 使用同序列其余位点;global 使用全体位点作为背景(默认 within)  
correction = fdr_bh  
#多重检验校正方法(statsmodels 支持的方法),默认 fdr_bh  
#bonferroni:Bonferroni 校正  
#sidak:Sidak 校正  
#holm-sidak:Holm-Sidak 校正  
#holm:Holm 校正  
#simes-hochberg:Simes-Hochberg 校正  
#hommel:Hommel 校正  
#fdr_bh:Benjamini-Hochberg FDR 校正(默认)  
#fdr_by:Benjamini-Yekutieli FDR 校正  
#fdr_tsbh:Two-stage Benjamini-Hochberg FDR 校正  
#fdr_tsbky:Two-stage Benjamini-Krieger-Yekutieli FDR 校正  

功能执行  
命令：  
```bash
CentriVision -ic total.conf
```
---

### -sa SATAGE 分歧推断

计算Monomer之间的相对距离，不指定祖先序列，只展示序列相对差异大小，理论上数值小可以代表最近扩增。  

查看参数：  
```bash
CentriVision -sa ?
```

参数重定向到配置文件total.conf

覆盖式命令：  
```bash
CentriVision -sa ? > total.conf
```

追加式命令：  
```bash
CentriVision -sa ? >> total.conf
```

配置文件：  
```
[SatAge]
genome_fa = centromere fasta
blast = True
monomers_fa = monomer fasta
blast_output = blast_results.tsv
blast_hits_file = monomer_blast_hits.png
age = True
monomer_age_file = monomer_age.tsv
age_plot_file = monomer_age.png
# monomer 环境窗口
window_size = 10000
# kmer 长度
k = 5

# 可调参数：标题和坐标轴字体大小
distance_threshold = 0.6
chrom_label = 20
xlabel = 18
xtick = 16
colorbar_label = 18
legend_fontsize = 18
discrete_colormap = False
n_bins = 8
mismatches = 5
```
参数详解：  
genome_fa = centromere fasta 着丝粒fasta序列  
blast = True 需要blast，若已有blast结果则设置为False  
monomers_fa = monomer fasta 重复单元种子文件  
blast_output = blast_results.tsv blast比对结果  
blast_hits_file = monomer_blast_hits.png 重复单元分布图  
age = True True表示需要推断相对突变距离，若有其它计算方法得到则设置为False  
monomer_age_file = monomer_age.tsv 相对突变距离文件  
age_plot_file = monomer_age.png 带年龄的重复单元分布图  
window_size = 10000 monomer 环境窗口  
k = 5 kmer 长度  
distance_threshold = 0.6 序列相似性或距离的阈值  
chrom_label = 20 染色体标签（或其它分类标签）字体大小  
xlabel = 18 X 轴标题的字体大小  
xtick = 16 X 轴刻度字体大小  
colorbar_label = 18 色条标签的字体大小  
legend_fontsize = 18 图例文字的字体大小  
discrete_colormap = False 是否使用离散配色  
n_bins = 8 离散配色时，颜色条分成多少个颜色块  
mismatches = 5 允许的最大不匹配数  

功能执行  
命令：  
```bash
CentriVision -sa total.conf
```
---

### -cd COMMUNITY_DETECTION 单元聚类

根据md模块得到的每个切片的采样进行聚类，划分不同的重复单元组。  

查看参数：  
```bash
CentriVision -cd ?
```

参数重定向到配置文件total.conf

覆盖式命令：  
```bash
CentriVision -cd ? > total.conf
```

追加式命令：  
```bash
CentriVision -cd ? >> total.conf
```

配置文件：  
```
[Community_detection]
fasta_file = fatsa file
gap = 10
identity = 75
alignment = 75
out_file = out community
```
参数详解：  
fasta_file = fatsa file md模块得到的种子文件  
gap = 10 比对时允许存在的最大gap  
identity = 75 分组最低相似度  
alignment = 75 分组最低匹配长度比例  
out_file = out community 社区文件  

功能执行  
命令：  
```bash
CentriVision -cd total.conf
```
---

### -cm REPEAT_COMMUNITY_MAP 类型分布

cd聚类得到的社区在所有着丝粒上的分布。  

![md](https://github.com/user-attachments/assets/7f54e0f2-a817-45f2-bf61-b4d5e531add0)  

查看参数：  
```bash
CentriVision -cm ?
```

参数重定向到配置文件total.conf

覆盖式命令：  
```bash
CentriVision -cm ? > total.conf
```

追加式命令：  
```bash
CentriVision -cm ? >> total.conf
```

配置文件：  
```
[Repeat_community_map]
community_file = community file
lens = lens file
repeat_gff = CentriVision repeat gff
focus_areas = areas_file1,areas_file2 or None
focus_name = areas_name1:red,areas_name2:green or None
model = global or local
top = 10
windows = 50000
step = 5000
savefile = save file (\*.png, \*.pdf, \*.svg)
```
参数详解：  
community_file = community file cd模块的到的社区文件  
lens = lens file 着丝粒的lens文件  
repeat_gff = CentriVision repeat gff cd模块输出的切片gff文件  
focus_areas = areas_file1,areas_file2 or None 关注区域(一般时ChIP先验着丝粒区域)可以输入多个文件  
focus_name = areas_name1:red,areas_name2:green or None 配色设置，每个文件对应一个名字一个颜色  
model = global or local 关注全局前top个社区还是单条着丝粒的前top个社区  
top = 10 只可视化数量靠前的社区分布  
windows = 50000 滑窗大小  
step = 5000 滑窗每次移动距离  
savefile = save file (\*.png, \*.pdf, \*.svg) 可视化输出  

功能执行  
命令：  
```bash
CentriVision -cm total.conf
```
---

### -gc GET_CENTRI 提取区间

根据gff提取对应的区域。  

查看参数：  
```bash
CentriVision -gc ?
```

参数重定向到配置文件total.conf

覆盖式命令：  
```bash
CentriVision -gc ? > total.conf
```

追加式命令：  
```bash
CentriVision -gc ? >> total.conf
```

配置文件：  
```
[Get_centri]
genome_file = genome file
# 确保着丝粒gff中的染色体和基因组染色体一致
gff_file = centri gff
out_fasta = out fasta
```
参数详解：  
genome_file = genome file 基因组fasta文件  
gff_file = centri gff 着丝粒位置文件  
out_fasta = out fasta 提取结果  

功能执行  
命令：  
```bash
CentriVision -gc total.conf
```
---

### -gf GET_CENTGFF 提取注释

从全基因组注释文件中提取着丝粒区域的注释信息。  

查看参数：  
```bash
CentriVision -gf ?
```

参数重定向到配置文件total.conf

覆盖式命令：  
```bash
CentriVision -gf ? > total.conf
```

追加式命令：  
```bash
CentriVision -gf ? >> total.conf
```

配置文件：  
```
[Get_centgff]
centromere_file = centri file
gff_file = genome gff
start = 3
end = 4
locmin = 9
output_file = centri gff
```
参数详解：  
centromere_file = centri file 着丝粒位置文件  
gff_file = genome gff 基因组注释文件gff/gff3格式  
start = 3 起始位置列号  
end = 4 终止位置列号  
locmin = 9 最大列数  
output_file = centri gff 输出结果  

功能执行  
命令：  
```bash
CentriVision -gf total.conf
```
---

### -gr GET_REPEAT 提取重复

从常用注释软件的结果中提取指定信息，例如提取Gypsy注释/Copia注释。  

查看参数：  
```bash
CentriVision -gr ?
```

参数重定向到配置文件total.conf

覆盖式命令：  
```bash
CentriVision -gr ? > total.conf
```

追加式命令：  
```bash
CentriVision -gr ? >> total.conf
```

配置文件：  
```
[Get_repeat]
genome = genome file
repeat_gff3 = repeat file
idtag = ID
classtag = Classification
out_path = out path
```
参数详解：  
genome = genome file 基因组文件  
repeat_gff3 = repeat file 重复序列注释文件  
idtag = ID 提取目标标签  
classtag = Classification 目标分类标记Gypsy/Copia等  
out_path = out path 输出文件  

功能执行  
命令：  
```bash
CentriVision -gr total.conf
```
---

## 文献引用
### Citating CentriVision
If you use CentriVision in your work, please cite:

Mei-Fang Lan, Xi-Yin Wang, Xian-Chun Zhang.2026.CentriVision: An integrated platform for multiscale centromere analysis in plants,Plant Communications,7(2):101689.    https://doi.org/10.1016/j.xplc.2025.101689.

***
$$
