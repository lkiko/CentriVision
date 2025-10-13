# -*- encoding: utf-8 -*-
'''
@File        :run.py
@Time        :2021/09/28 09:04:51
@Author        :charles kiko
@Version        :1.0
@Contact        :charles_kiko@163.com
@Desc        :主程序
'''


import argparse
import os
import sys
import configparser
import pandas as pd
import CentriVision
import CentriVision.bez as bez


from CentriVision.Palindromic import Palindromic
from CentriVision.Centrifinder import Centrifinder
from CentriVision.TRF import TRF
from CentriVision.Heatmap import Heatmap
from CentriVision.Monomer import Monomer
from CentriVision.Dotplot import Dotplot
from CentriVision.Get_centri import Get_centri
from CentriVision.Get_repeat import Get_repeat
from CentriVision.Get_centgff import Get_centgff
from CentriVision.Count_file import Count_file
from CentriVision.HOR import HOR
from CentriVision.SeqSigil import SeqSigil
from CentriVision.Ic_Significance import Ic_Significance
from CentriVision.SatAge import SatAge
from CentriVision.EdistDot import EdistDot
from CentriVision.EdistAln import EdistAln
from CentriVision.Community_detection import Community_detection
from CentriVision.Repeat_community_map import Repeat_community_map



parser = argparse.ArgumentParser(
    prog = 'CentriVision', usage = '%(prog)s [options]', epilog = "", formatter_class = argparse.RawDescriptionHelpFormatter,)
parser.description = '''\
runing CentriVision
    -------------------------------------- '''
parser.add_argument("-v", "--version", action = 'version', version='1.0.1')



parser.add_argument("-ps", dest = "Palindromic",
                    help = "Palindromic sequence 查询基因组中的回文序列")
parser.add_argument("-trf", dest = "TRF",
                    help = "run TRF(Tandem Repeat Finder) 通过TRF查找串联重复序列；")
parser.add_argument("-cf", dest = "Centrifinder",
                    help = "Centrifinder 着丝粒预测；")
parser.add_argument("-md", dest = "Dotplot",
                    help = "mini Dotplot 重复序列点图；")
parser.add_argument("-hm", dest = "Heatmap",
                    help = "Heatmap 区域相似度热图；")
parser.add_argument("-m", dest = "Monomer",
                    help = "Monomer scanning 重复单体扫描；")
parser.add_argument("-s", dest = "SeqSigil",
                    help = "SeqSigil scanning 重复单体logo；")
parser.add_argument("-ic", dest = "Ic_Significance",
                    help = "Ic Significance  单体保守性IC检验；")
parser.add_argument("-sa", dest = "SatAge",
                    help = "SatAge Monomer 重复时间推断（拟分子钟）；")
parser.add_argument("-gc", dest = "Get_centri",
                    help = "Get_centri 提取基因组的指定区域；")
parser.add_argument("-gf", dest = "Get_centgff",
                    help = "Get_centgff 提取基因组的指定区域gff,index修改为相对着丝粒；")
parser.add_argument("-gr", dest = "Get_repeat",
                    help = "Get_repeat 根据gff3提取基因组的重复序列；")
parser.add_argument("-c", dest = "Count_file",
                    help = "Count_file 统计dotplot文件；")
parser.add_argument("-r", dest = "HOR",
                    help = "HOR HOR搜索；")
parser.add_argument("-ed", dest = "EdistDot",
                    help = "EdistDot EdistDot 点阵图；")
parser.add_argument("-e", dest = "EdistAln",
                    help = "EdistAln EdistAln 快速比对；")
parser.add_argument("-cd", dest = "Community_detection",
                    help = "Community_detection 重复序列社区发现；")
parser.add_argument("-cm", dest = "Repeat_community_map",
                    help = "Repeat_community_map 重复序列社区映射；")
args = parser.parse_args()



def run_Palindromic():
    options = bez.load_conf(args.Palindromic, 'Palindromic')
    Palindromic1 = Palindromic(options)
    Palindromic1.run()

def run_Centrifinder():
    options = bez.load_conf(args.Centrifinder, 'Centrifinder')
    Centrifinder1 = Centrifinder(options)
    Centrifinder1.run()

def run_TRF():
    options = bez.load_conf(args.TRF, 'TRF')
    TRF1 = TRF(options)
    TRF1.run()

def run_Heatmap():
    options = bez.load_conf(args.Heatmap, 'Heatmap')
    Heatmap1 = Heatmap(options)
    Heatmap1.run()

def run_Monomer():
    options = bez.load_conf(args.Monomer, 'Monomer')
    Monomer1 = Monomer(options)
    Monomer1.run()

def run_SeqSigil():
    options = bez.load_conf(args.SeqSigil, 'SeqSigil')
    SeqSigil1 = SeqSigil(options)
    SeqSigil1.run()

def run_Ic_Significance():
    options = bez.load_conf(args.Ic_Significance, 'Ic_Significance')
    Ic_Significance1 = Ic_Significance(options)
    Ic_Significance1.run()

def run_SatAge():
    options = bez.load_conf(args.SatAge, 'SatAge')
    SatAge1 = SatAge(options)
    SatAge1.run()

def run_EdistDot():
    options = bez.load_conf(args.EdistDot, 'EdistDot')
    EdistDot1 = EdistDot(options)
    EdistDot1.run()

def run_EdistAln():
    options = bez.load_conf(args.EdistAln, 'EdistAln')
    EdistAln1 = EdistAln(options)
    EdistAln1.run()

def run_Get_centri():
    options = bez.load_conf(args.Get_centri, 'Get_centri')
    Get_centri1 = Get_centri(options)
    Get_centri1.run()

def run_Get_centgff():
    options = bez.load_conf(args.Get_centgff, 'Get_centgff')
    Get_centgff1 = Get_centgff(options)
    Get_centgff1.run()

def run_Get_repeat():
    options = bez.load_conf(args.Get_repeat, 'Get_repeat')
    Get_repeat1 = Get_repeat(options)
    Get_repeat1.run()

def run_Count_file():
    options = bez.load_conf(args.Count_file, 'Count_file')
    Count_file1 = Count_file(options)
    Count_file1.run()

def run_HOR():
    options = bez.load_conf(args.HOR, 'HOR')
    HOR1 = HOR(options)
    HOR1.run()

def run_Community_detection():
    options = bez.load_conf(args.Community_detection, 'Community_detection')
    Community_detection1 = Community_detection(options)
    Community_detection1.run()

def run_Repeat_community_map():
    options = bez.load_conf(args.Repeat_community_map, 'Repeat_community_map')
    Repeat_community_map1 = Repeat_community_map(options)
    Repeat_community_map1.run()

def run_Dotplot():
    options = bez.load_conf(args.Dotplot, 'Dotplot')
    Dotplot1 = Dotplot(options)
    Dotplot1.run()


def module_to_run(argument):
    switcher = {


        'Palindromic': run_Palindromic,
        'Centrifinder': run_Centrifinder,
        'TRF': run_TRF,
        'Heatmap': run_Heatmap,
        'Monomer': run_Monomer,
        'SeqSigil': run_SeqSigil,
        'Ic_Significance': run_Ic_Significance,
        'SatAge': run_SatAge,
        'Get_centri': run_Get_centri,
        'Get_centgff': run_Get_centgff,
        'Get_repeat': run_Get_repeat,
        'Count_file': run_Count_file,
        'HOR': run_HOR,
        'EdistDot': run_EdistDot,
        'EdistAln': run_EdistAln,
        'Community_detection': run_Community_detection,
        'Repeat_community_map': run_Repeat_community_map,
        'Dotplot': run_Dotplot,

    }
    return switcher.get(argument)()

def main():
    path = CentriVision.__path__[0]
    options = {

               'Palindromic': 'Palindromic.conf',
               'Centrifinder': 'Centrifinder.conf',
               'TRF': 'TRF.conf',
               'Heatmap': 'Heatmap.conf',
               'Monomer': 'Monomer.conf',
               'SeqSigil': 'SeqSigil.conf',
               'Ic_Significance': 'Ic_Significance.conf',
               'SatAge': 'SatAge.conf',
               'Get_centri': 'Get_centri.conf',
               'Get_centgff': 'Get_centgff.conf',
               'Get_repeat': 'Get_repeat.conf',
               'Count_file': 'Count_file.conf',
               'EdistAln': 'EdistAln.conf',
               'HOR': 'HOR.conf',
               'EdistDot': 'EdistDot.conf',
               'Community_detection': 'Community_detection.conf',
               'Repeat_community_map': 'Repeat_community_map.conf',
               'Dotplot': 'Dotplot.conf',

               }
    for arg in vars(args):
        value = getattr(args, arg)
        # print(value)
        if value is not None:
            if value in ['?', 'help', 'example']:
                f = open(os.path.join(path, 'example', options[arg]), encoding='utf-8')
                print(f.read())
            elif value == 'e':
                out = '''\
        File example
        [fpchrolen]
        chromosomes number_of_bases
        *   *
        *   *
        *   *
        [fpgff]
        chromosomes gene    start   end
        *   *   *   *
        *   *   *   *
        *   *   *   *
        [fpgenefamilyinf]
        gene1   gene2   Ka  Ks
        *   *   *   *
        *   *   *   *
        *   *   *   *
        [alphagenepairs]
        gene1   gene2
        *   *   *
        *   *   *
        *   *   *

        The file columns are separated by Tab
        -----------------------------------------------------------    '''
                print(out)
            elif not os.path.exists(value):
                print(value+' not exits')
                sys.exit(0)
            else:
                module_to_run(arg)

