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

### 着丝粒热图 -hm HEATMAP ${\color{orange}\textbf{热图}}$  
着丝粒热图  

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

### 快速比对 -e EDISTALN ${\color{orange}\textbf{整体点图比对}}$  
拆分着丝粒并快速比对。  

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


### 大型点图绘制 -ed EDISTDOT ${\color{orange}\textbf{整体点图}}$  
比对结果绘制。  

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
![chr02_1_s98-dotplot](https://github.com/lkiko/CentriVision/assets/57740432/47a3700d-c49d-4000-8f67-3ee4c7c1ffaa)

查找重复单元的相位纠正图
![chr02_1_s98](https://github.com/lkiko/CentriVision/assets/57740432/0de54b96-94bd-4c05-a60a-b60e64dd8020)

${\color{red}\textbf{切片大小}}$ 与分辨率和计算机内存大小挂钩，大型矩阵极其消耗内存；对于具有 ${\color{red}\textbf{超大着丝粒}}$ 的物种，切片数量非常多，是否需要输出所有自相似矩阵图以及比对矩阵需要适当选择，可利用输出文件可选的生成对应切片的自相似矩阵图和比对矩阵

### 统计 -c COUNT_FILE ${\color{orange}\textbf{统计绘图}}$  
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


### 重复单体拆分 -m MONOMER ${\color{orange}\textbf{拆分单体}}$  
将较为均匀的重复区域拆分为单体。  

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

### 单体模式查询 -s SEQSIGIL ${\color{orange}\textbf{重复单体logo图}}$  
拆分着丝粒并扫描重复序列。  

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

### 重复单体局部保守性 -ic IC_SIGNIFICANCE ${\color{orange}\textbf{香农信息熵}}$  
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

### 相对突变距离 -sa SATAGE ${\color{orange}\textbf{Monomer相对年龄}}$  
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

### 切片种子序列聚类 -cd COMMUNITY_DETECTION ${\color{orange}\textbf{聚类}}$  
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

### 社区在所有着丝粒上的分布 -cm REPEAT_COMMUNITY_MAP ${\color{orange}\textbf{社区分布}}$  
cd聚类得到的社区在所有着丝粒上的分布。  

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

### 提取着丝粒 -gc GET_CENTRI ${\color{orange}\textbf{提取序列}}$  
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

### 提取着丝粒gff -gf GET_CENTGFF ${\color{orange}\textbf{提取着丝粒区域gff}}$  
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

### 提取重复信息 -gr GET_REPEAT ${\color{orange}\textbf{提取repeatmasker/EDTA文件}}$  
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

***
