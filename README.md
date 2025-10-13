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

你可以使用 pip 来安装 CentriVision：

`pip install CentriVision`

或者，你也可以从本地安装 CentriVision 的 wheel 文件：

`pip install CentriVision-x.x.x-py3-none-any.whl`

软件依赖TRF,Mafft,Muscle

---
## 配置：

**以下以 ubuntu操作系统 用户名为charles 的conda3 python3.11环境为例**


默认安装路径：

`/home/charles/anaconda3/lib/python3.11/site-packages/CentriVision/`

查看依赖文件：

`cat /home/charles/anaconda3/lib/python3.11/site-packages/CentriVision/conf.ini`

[ini] 

trf_path = /usr/bin/trf 

mafft_path = /usr/bin/mafft 

muscle_path = /usr/bin/muscle 


使用vim或其它编辑器修改对应依赖软件位置TRF\Mafft\Muscle 并保存

---

## 使用
命令： `CentriVision -h`

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


***

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
