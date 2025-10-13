# -*- encoding: utf-8 -*-
'''
@File        :bez.py
@Time        :2021/09/28 11:26:28
@Author        :charles kiko
@Version        :1.0
@Contact        :charles_kiko@163.com
@Desc        :None
'''

import configparser
import os
import re
import math
from math import *
import numpy as np
import pandas as pd
import CentriVision
from Bio import Seq, SeqIO, SeqRecord
import codecs
from tqdm import trange
import gc
import matplotlib.pyplot as plt
from matplotlib.patches import *
from matplotlib.patches import Circle, Ellipse
from pylab import *
from collections import Counter
import time

def TRF2gff(file,gff3,fasta,path):
    dat = open(file,'r').read()
    dat = dat.split('Sequence: ')
    # print(len(dat))
    gff = open(path+gff3,'w')
    TRF = open(path+fasta,'w')
    gff.write('')
    chro = ''
    num = 0
    for i in trange(1,len(dat)):
        lt = dat[i].split('\n')
        #print(len(lt))
        for line in lt:
            lt0 = line.split()
            if len(lt0) == 1:
                chro = lt0[0]
            elif len(lt0) == 0 or lt0[0] == 'Parameters:':
                continue
            else:
                # print(len(lt0))
                num += 1
                ID = 'TRF'+str(num).zfill(5)
                out_l = [chro,'TRF','TandemRepeat',lt0[0],lt0[1],'.','.','.',"ID="+ID+";PeriodSize="+str(lt0[2])+';CopyNumber='+str(lt0[3])+';Consensus='+lt0[-2]]
                # print(out_l)
                ID = '>'+ID+'#TRF'
                seq = lt0[-1]
                gff.write('\t'.join(out_l)+"\n")
                TRF.write(ID+'\n'+seq+'\n')
    gff.close()
    TRF.close()

def read_zsl(file):
    zsl = {}
    for line in open(file,'r'):
        lt = line.strip('\n').split('\t')
        # lt[0] = 'chr'+lt[0][3:]
        if lt[0] == 'Chr_ID':
            continue
        if lt[0] not in zsl.keys():
            zsl[lt[0]] = []
            zsl[lt[0]].append([int(lt[1]),int(lt[2]),abs(int(lt[1])-int(lt[2]))+1])
        else:
            zsl[lt[0]].append([int(lt[1]),int(lt[2]),abs(int(lt[1])-int(lt[2]))+1])
    return zsl

def readblast_large(genepairs_,blast_reverse,score_levels_,position_,levels_,tandem_,gff,lens,chr1_start,chr2_start,step1,step2):
    print("Count the number of lines in the BLAST file")
    start_time = time.time()  # 记录开始时间
    total_lines = sum(1 for _ in open(genepairs_,'r'))
    end_time = time.time()  # 记录结束时间
    print(f"Function execution time: {end_time - start_time:.2f} seconds")
    x1,x2,x3,y1,y2,y3 = [],[],[],[],[],[]
    # 初始化进度条

    with trange(total=math.ceil(total_lines/100000), desc="Reading BLAST", unit="block") as pbar:
        # 分块读取文件
        chunk_iter = pd.read_csv(genepairs_,header = None, sep='\t', comment='#', chunksize=100000)
        processed_lines = 0
        tail = pd.DataFrame()
        for chunk in chunk_iter:
            processed_lines += len(chunk)
            chunk.loc[:, 12] = chunk.apply(lambda x : cat_gene(x[0],x[1]), axis=1)
            chunk = chunk.drop_duplicates(subset=12 , keep='first')
            del chunk[12]
            if blast_reverse == 'True':
                chunk[[0, 1]] = chunk[[1, 0]]
            chunk = chunk.drop(chunk[(chunk[0] == chunk[1])].index)#.sort_values(by=[0,11],ascending= [True,False])
            last_key = list(chunk[0].to_list())[-1]
            if not tail.empty:
                chunk = pd.concat([tail, chunk], ignore_index=True)
                chunk.reset_index(drop=True,inplace=True)
            dic = chunk.groupby(0).groups# 按照第几列分组
            # 判断是否为最后一块
            if processed_lines == total_lines:
                # print("This is the last chunk.")
                tail = pd.DataFrame()
            else:
                tail = chunk.loc[dic[last_key]]
                tail.reset_index(drop=True,inplace=True)
                del dic[last_key]

            for key in dic.keys():
                # print(" ******** ",key," ******** ")
                local = chunk.loc[dic[key]]
                local.reset_index(drop=True,inplace=True)

                if key not in gff.keys():
                    continue
                if score_levels_ == 'False':
                    if tandem_ == 'True':
                        object0 = local[1].to_list()
                        object0 = list(dict.fromkeys(object0))# 去除tandem保留顺序
                        for i in range(levels_[0]):
                            if i < len(object0):
                                if object0[i] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[object0[i]]['chr'] not in lens.keys():
                                    continue
                                if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                    x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                                else:
                                    x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                                if 'rev' in lens[gff[object0[i]]['chr']].keys() and lens[gff[object0[i]]['chr']]['rev'] == '-':
                                    y = -(chr2_start[gff[object0[i]]['chr']] + lens[gff[object0[i]]['chr']][position_] - gff[object0[i]][position_] + 1)*step2
                                else:
                                    y = -(chr2_start[gff[object0[i]]['chr']] + gff[object0[i]][position_])*step2
                                x1.append(x)
                                y1.append(y)
                        for i in range(levels_[0],levels_[0]+levels_[1]):
                            if i < len(object0):
                                if object0[i] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[object0[i]]['chr'] not in lens.keys():
                                    continue
                                if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                    x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                                else:
                                    x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                                if 'rev' in lens[gff[object0[i]]['chr']].keys() and lens[gff[object0[i]]['chr']]['rev'] == '-':
                                    y = -(chr2_start[gff[object0[i]]['chr']] + lens[gff[object0[i]]['chr']][position_] - gff[object0[i]][position_] + 1)*step2
                                else:
                                    y = -(chr2_start[gff[object0[i]]['chr']] + gff[object0[i]][position_])*step2
                                x2.append(x)
                                y2.append(y)
                        for i in range(levels_[0]+levels_[1],levels_[0]+levels_[1]+levels_[2]):
                            if i < len(object0):
                                if object0[i] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[object0[i]]['chr'] not in lens.keys():
                                    continue
                                if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                    x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                                else:
                                    x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                                if 'rev' in lens[gff[object0[i]]['chr']].keys() and lens[gff[object0[i]]['chr']]['rev'] == '-':
                                    y = -(chr2_start[gff[object0[i]]['chr']] + lens[gff[object0[i]]['chr']][position_] - gff[object0[i]][position_] + 1)*step2
                                else:
                                    y = -(chr2_start[gff[object0[i]]['chr']] + gff[object0[i]][position_])*step2
                                x3.append(x)
                                y3.append(y)
                    else:
                        object0 = local[1].to_list()
                        for i in range(levels_[0]):
                            if i < len(object0):
                                if object0[i] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[object0[i]]['chr'] not in lens.keys():
                                    continue
                                if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                    x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                                else:
                                    x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                                if 'rev' in lens[gff[object0[i]]['chr']].keys() and lens[gff[object0[i]]['chr']]['rev'] == '-':
                                    y = -(chr2_start[gff[object0[i]]['chr']] + lens[gff[object0[i]]['chr']][position_] - gff[object0[i]][position_] + 1)*step2
                                else:
                                    y = -(chr2_start[gff[object0[i]]['chr']] + gff[object0[i]][position_])*step2
                                x1.append(x)
                                y1.append(y)
                        for i in range(levels_[0],levels_[0]+levels_[1]):
                            if i < len(object0):
                                if object0[i] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[object0[i]]['chr'] not in lens.keys():
                                    continue
                                if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                    x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                                else:
                                    x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                                if 'rev' in lens[gff[object0[i]]['chr']].keys() and lens[gff[object0[i]]['chr']]['rev'] == '-':
                                    y = -(chr2_start[gff[object0[i]]['chr']] + lens[gff[object0[i]]['chr']][position_] - gff[object0[i]][position_] + 1)*step2
                                else:
                                    y = -(chr2_start[gff[object0[i]]['chr']] + gff[object0[i]][position_])*step2
                                x2.append(x)
                                y2.append(y)
                        for i in range(levels_[0]+levels_[1],levels_[0]+levels_[1]+levels_[2]):
                            if i < len(object0):
                                if object0[i] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[object0[i]]['chr'] not in lens.keys():
                                    continue
                                if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                    x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                                else:
                                    x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                                if 'rev' in lens[gff[object0[i]]['chr']].keys() and lens[gff[object0[i]]['chr']]['rev'] == '-':
                                    y = -(chr2_start[gff[object0[i]]['chr']] + lens[gff[object0[i]]['chr']][position_] - gff[object0[i]][position_] + 1)*step2
                                else:
                                    y = -(chr2_start[gff[object0[i]]['chr']] + gff[object0[i]][position_])*step2
                                x3.append(x)
                                y3.append(y)
                else:
                    score_levels = [float(i) for i in score_levels_.split(',')]
                    for index,row in local.iterrows():
                        if row[2] >= score_levels[0]:
                            if row[1] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[row[1]]['chr'] not in lens.keys():
                                continue
                            if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                            else:
                                x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                            if 'rev' in lens[gff[row[1]]['chr']].keys() and lens[gff[row[1]]['chr']]['rev'] == '-':
                                y = -(chr2_start[gff[row[1]]['chr']] + lens[gff[row[1]]['chr']][position_] - gff[row[1]][position_] + 1)*step2
                            else:
                                y = -(chr2_start[gff[row[1]]['chr']] + gff[row[1]][position_])*step2
                            x1.append(x)
                            y1.append(y)
                        elif row[2] >= score_levels[1]:
                            if row[1] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[row[1]]['chr'] not in lens.keys():
                                continue
                            if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                            else:
                                x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                            if 'rev' in lens[gff[row[1]]['chr']].keys() and lens[gff[row[1]]['chr']]['rev'] == '-':
                                y = -(chr2_start[gff[row[1]]['chr']] + lens[gff[row[1]]['chr']][position_] - gff[row[1]][position_] + 1)*step2
                            else:
                                y = -(chr2_start[gff[row[1]]['chr']] + gff[row[1]][position_])*step2
                            x2.append(x)
                            y2.append(y)
                        elif row[2] >= score_levels[2]:
                            if row[1] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[row[1]]['chr'] not in lens.keys():
                                continue
                            if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                            else:
                                x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                            if 'rev' in lens[gff[row[1]]['chr']].keys() and lens[gff[row[1]]['chr']]['rev'] == '-':
                                y = -(chr2_start[gff[row[1]]['chr']] + lens[gff[row[1]]['chr']][position_] - gff[row[1]][position_] + 1)*step2
                            else:
                                y = -(chr2_start[gff[row[1]]['chr']] + gff[row[1]][position_])*step2
                            x3.append(x)
                            y3.append(y)
            # pbar.update(len(chunk))
            pbar.update()
            
    return x1,x2,x3,y1,y2,y3

def readcolineartly_dotplot(file,blockmin,step,genepairsfile_type,gff,chr_list):# 块首尾
    block = []
    if genepairsfile_type == 'MCScanX':
        alphagenepairs = open(file, 'r', encoding='utf-8').read()
        file = alphagenepairs.split('## Alignment ')[1:]
        for local in file:
            local = local.strip('\n').split('\n')[1:]
            if len(local) < blockmin:
                continue
            else:
                data = [[i.split('\t')[1],i.split('\t')[2]] for i in local]
                if data[0][0] not in gff.keys() or data[0][1] not in gff.keys():
                    continue
                chr1,chr2 = gff[data[0][0]]['chr'],gff[data[0][1]]['chr']
                if chr1 not in chr_list or chr2 not in chr_list:
                    continue
                if chr1 == chr2 and min(abs(gff[data[0][0]]['order'] - gff[data[0][1]]['order']),abs(gff[data[0][0]]['order'] - gff[data[-1][1]]['order']),abs(gff[data[-1][0]]['order'] - gff[data[0][1]]['order']),abs(gff[data[-1][0]]['order'] - gff[data[-1][1]]['order'])) > step:
                    block.append(data)
                elif chr1 != chr2:
                    block.append(data)
    elif genepairsfile_type == 'WGDI':
        alphagenepairs = open(file, 'r', encoding='utf-8').read()
        file = alphagenepairs.split('# Alignment ')[1:]
        for local in file:
            local = local.strip('\n').split('\n')[1:]
            if len(local) < blockmin:
                continue
            else:
                data = [[i.split()[0],i.split()[2]] for i in local]
                if data[0][0] not in gff.keys() or data[0][1] not in gff.keys():
                    continue
                chr1,chr2 = gff[data[0][0]]['chr'],gff[data[0][1]]['chr']
                if chr1 not in chr_list or chr2 not in chr_list:
                    continue
                if chr1 == chr2 and min(abs(gff[data[0][0]]['order'] - gff[data[0][1]]['order']),abs(gff[data[0][0]]['order'] - gff[data[-1][1]]['order']),abs(gff[data[-1][0]]['order'] - gff[data[0][1]]['order']),abs(gff[data[-1][0]]['order'] - gff[data[-1][1]]['order'])) > step:
                    block.append(data)
                elif chr1 != chr2:
                    block.append(data)

    else:
        print('genepairsfile_type error: File Format not recognized! 共线性文件读取方式为完整block 2023.8.21修改')
        exit()
    return block


def get_scorefile():
    score = CentriVision.__path__[0]
    return score+'/example/score.xlsx'

def get_path():
    from pathlib import Path
    # 获取当前脚本的路径
    script_path = Path(__file__).resolve()
    # 获取当前脚本所在的目录
    script_dir = script_path.parent
    # print(f"Script path: {script_path}")
    # print(f"Script directory: {script_dir}")
    return script_dir

# def config():
#     conf = configparser.ConfigParser()
#     conf.read(os.path.join(CentriVision.__path__[0], 'conf.ini'))
#     return conf.items('ini')

# def load_conf(file, section):
#     conf = configparser.ConfigParser()
#     conf.read(file)
#     return conf.items(section)

def config():
    conf = configparser.ConfigParser()
    conf.read(os.path.join(CentriVision.__path__[0], 'conf.ini'), encoding='utf-8')
    return conf.items('ini')

def load_conf(file, section):
    conf = configparser.ConfigParser()
    conf.read(file, encoding='utf-8')
    return conf.items(section)


    # 多进程错误打印
def print_error(value):
    print("Process pool error, the cause of the error is :", value)

def read_lens(lens):
    chr_dic = {}
    chr_list = []
    lens = pd.read_csv(lens,header = None, sep='\t', comment='#')
    chr_list = lens[0].to_list()
    for index,row in lens.iterrows():
        dic0 = {}
        dic0['end'] = row[1]
        try:
            dic0['order'] = row[2]
        except:
            pass
        # dic0['order'] = row[2]
        try:
            dic0['rev'] = row[3]
        except:
            pass
        chr_dic[row[0]] = dic0
    return chr_dic,chr_list

def read_gff(gff,chr_list):
    gene = {}
    gff = pd.read_csv(gff,header = None, sep='\t', comment='#')
    gff[0] = gff[0].astype(str)
    # print(gff)
    for index,row in gff.iterrows():
        if row[0] not in chr_list:
            continue
        dic0 = {}
        dic0['chr'] = (row[0])
        dic0['end'] = row[2]
        dic0['order'] = row[5]
        gene[row[1]] = dic0
    return gene