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
-->

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


${\color{orange}\textbf{CentriVision}}$ 是一个用于研究 ${\color{orange}\textbf{着丝粒}}$ 结构的软件工具。
---

## 简介

${\color{orange}\textbf{CentriVision}}$ 旨在提供一个简单而强大的工具，用于分析和可视化基因组中着丝粒的结构。它支持从基因组中提取和分析着丝粒相关信息，并提供丰富的可视化功能，以帮助研究人员深入理解着丝粒的组织和功能。同时支持 ${\color{orange}\textbf{植物和动物}}$ 基因组。  
---

## 安装方法  
<p align="center">
  <img src="https://github.com/lkiko/CentriVision/blob/main/video/install.gif?raw=true" width="100%">
</p>

你可以使用  ${\color{green}\textbf{pip}}$ (https://pypi.org/project/CentriVision/) 来安装 CentriVision：  
```bash
pip install CentriVision
```

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
## 配置：  
**以下以 ubuntu操作系统 用户名为charles 的miniconda3 python3.13环境为例**  

<p align="center">
  <img src="https://github.com/lkiko/CentriVision/blob/main/video/configuration.gif?raw=true" width="100%">
</p>

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
which mafft
which bowtie2
which mpirun
```

---

## 使用方法
命令：  
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
***
### 着丝粒鉴定 -trf TRF ${\color{orange}\textbf{TRF(Tandem Repeat Finder)}}$  
调用TRF(Tandem Repeat Finder)扫描重复序列，根据重复序列判断着丝粒，也可以输入现有的重复注释gff3文件，同时适用于 ${\color{orange}\textbf{串联重复类型和转座子类型}}$ 的着丝粒。  
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
windows = 10000 重复序列密度窗口跨度  
step = 5000 重复序列密度窗口滑动步长  
gap = 40 重复区域连续性容错宽度 gap\*windows  
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

$${\color{red}\textbf{TRF在面对大区域重复的时候扫描特别慢，可以单独切片运行TRF注释}}$$  

#### 已有注释文件时使用 -cf CENTRIFINDER 模块，输入文件兼容 ${\color{orange}\textbf{串联重复注释文件、转座子注释文件或者回文序列注释}}$  

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

#### 回文序列注释 -ps PALINDROMIC 模块，通过染色体 ${\color{orange}\textbf{回文序列}}$ 密度来鉴定着丝粒，注释结果输入 -cf CENTRIFINDER 模块

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
step = 5000 密度窗口华东距离  
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

---

### 重复模式 -md DOTPLOT ${\color{orange}\textbf{分块点图}}$  
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

![md参数修改](https://github.com/user-attachments/assets/121aaaa2-f386-46f6-a6c7-6341faf6fde0)

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

genome_file = genome file 基因组fasta文件  
lens = lens file 染色体文件  
chip_seq = chip_seq map file or None ChIP-seq或其它数据的先验着丝粒位置文件  
colors = hish,centri,chip or hish,centri,None or #38b48b,#1e50a2,#d7003a 颜色配置  
trfgff = out gff TRF输出gff3结果  
trffasta = out fasta TRF输出fasta文件  
windows = 10000 重复序列密度窗口跨度  
step = 5000 重复序列密度窗口滑动步长  
gap = 40 重复区域连续性容错宽度 gap\*windows  
centrigff = centri gff 鉴定候选区结果  
centrifasta = centri fasta 候选区fasta文件  

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
功能执行  
命令：  
```bash
CentriVision -md total.conf
```

结果：  
点阵图：展示序列重复规律的点图
![chr02_1_s98-dotplot](https://github.com/lkiko/CentriVision/assets/57740432/47a3700d-c49d-4000-8f67-3ee4c7c1ffaa)

查找重复单元的相位纠正图
![chr02_1_s98](https://github.com/lkiko/CentriVision/assets/57740432/0de54b96-94bd-4c05-a60a-b60e64dd8020)

$${\color{red}\textbf{TRF在面对大区域重复的时候扫描特别慢，可以单独切片运行TRF注释}}$$  



***

### 着丝粒热图 -hmap


查看对应参数

命令：  
```bash
CentriVision -hmap ?
```

将参数重定向到配置文件total.conf

覆盖式命令：  
```bash
CentriVision -hmap ? > total.conf
``` 

追加式命令：  
```bash
CentriVision -hmap ? >> total.conf
``` 

配置文件：

[Heatmap] 

centromere_file = genome file 

align_software = muscle or mafft 

color_mode = Discrete or Gradient 

split = 1000 

out_path = out path 


参数解释： 

centromere_file 着丝粒fasta文件 

align_software 多序列比对软件支持 muscle or mafft 

color_mode 颜色模式支持 Discrete or Gradient 

split 每着丝粒拆分为1000等份 

out_path 输出路径 


功能执行：

命令：  
```bash
CentriVision -hmap total.conf
```

结果

重复序列相似性热图：
![s02_1](https://github.com/lkiko/CentriVision/assets/57740432/d0b95ae5-d83f-4997-9410-2768ddc296bf)


***
