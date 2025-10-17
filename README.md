[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15910908.svg)](https://doi.org/10.5281/zenodo.15910908)

# CentriVision

## At present there is no detailed manual for this application, you will simply have to play around and see what happens.
## 📖 用户手册正在紧锣密鼓地编写中！  
## I'm working hard on a clear and practical guide—stay tuned!
## 如有疑问或建议，随时提 Issue，我们一起把它打磨得更好。  
## Questions or ideas? Open an issue and let’s make it better together.

CentriVision 是一个用于研究基因组着丝粒结构的软件工具。
---

## 简介

CentriVision 旨在提供一个简单而强大的工具，用于分析和可视化基因组中着丝粒的结构。它支持从基因组中提取和分析着丝粒相关信息，并提供丰富的可视化功能，以帮助研究人员深入理解着丝粒的组织和功能。
---

## 安装
<p align="center">
  <img src="https://github.com/lkiko/CentriVision/blob/main/video/install.gif?raw=true" width="100%">
</p>

你可以使用 pip 来安装 CentriVision：

`pip install CentriVision`

或者，你也可以从本地安装 CentriVision 的 wheel 文件：

`pip install CentriVision-x.x.x-py3-none-any.whl`

软件依赖TRF,Mafft,Muscle,clustalw

---
## 配置：

**以下以 ubuntu操作系统 用户名为charles 的miniconda3 python3.13环境为例**

<p align="center">
  <img src="https://github.com/lkiko/CentriVision/blob/main/video/configuration.gif?raw=true" width="100%">
</p>

默认安装路径：

`/home/charles/miniconda3/lib/python3.13/site-packages/CentriVision/`

查看依赖文件：

`cat /home/charles/miniconda3/lib/python3.13/site-packages/CentriVision/conf.ini`

```
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

---

## 使用
命令： `CentriVision -h`

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
### 着丝粒鉴定 -trf
<p align="center">
  <img src="https://github.com/lkiko/CentriVision/blob/main/video/trf.gif?raw=true" width="100%">
</p>

查看参数：`CentriVision -trf ?`

![参数](https://github.com/user-attachments/assets/a4480953-31ee-461d-a35e-cd71fd3e5dbc)


参数重定向到配置文件total.conf

覆盖式命令：  `CentriVision -trf ? > total.conf`

追加式命令：  `CentriVision -trf ? >> total.conf`

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
功能执行：

命令： `CentriVision -trf total.conf`

<p align="center">
  <img src="https://github.com/lkiko/CentriVision/blob/main/video/trf-run.gif?raw=true" width="100%">
</p>
输出结果:

![图示](https://github.com/user-attachments/assets/3e2e19d7-7256-4cd9-9495-047d7f29a4be)
淡蓝色为重复序列分布图，蓝色区域为着丝粒候选区域

![结果](https://github.com/user-attachments/assets/83d37692-95db-4380-b960-89230164c18a)
鉴定结果判断


### 着丝粒点阵图 -d

查看对应参数

命令： `CentriVision -d ?`

将参数重定向到配置文件total.conf

覆盖式命令：  `CentriVision -d ? > total.conf`

追加式命令：  `CentriVision -d ? >> total.conf`

配置文件：

[Dotplot] 

genome_file = genome file 

windows = 4000 

minlength = 8 

poly = False 

cpu = 16 

outfile = out dotplot 



参数解释：

genome_file 着丝粒fasta文件 

windows 窗口宽度 

minlength 最短重复片段 

poly 单碱基重复去除 

cpu 线程数 

outfile 着丝粒特征文件 


功能执行：

命令： `CentriVision -d total.conf`

结果：
点阵图：展示序列重复规律的点图
![chr02_1_s98-dotplot](https://github.com/lkiko/CentriVision/assets/57740432/47a3700d-c49d-4000-8f67-3ee4c7c1ffaa)

查找重复单元的相位纠正图
![chr02_1_s98](https://github.com/lkiko/CentriVision/assets/57740432/0de54b96-94bd-4c05-a60a-b60e64dd8020)


***

### 着丝粒热图 -hmap


查看对应参数

命令： `CentriVision -hmap ?`

将参数重定向到配置文件total.conf

覆盖式命令：  `CentriVision -hmap ? > total.conf` 

追加式命令：  `CentriVision -hmap ? >> total.conf` 

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

命令： `CentriVision -hmap total.conf`

结果

重复序列相似性热图：
![s02_1](https://github.com/lkiko/CentriVision/assets/57740432/d0b95ae5-d83f-4997-9410-2768ddc296bf)


***
