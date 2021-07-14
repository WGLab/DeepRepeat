#!/usr/bin/env python

import re
import os
import sys
import string
import math
import random

import numpy
import argparse

import logging

from argparse import RawTextHelpFormatter

from scripts import Repeat_CNN_Prediction_motif_gen2
from scripts import Repeat_CNN_Prediction_model
from scripts import myHeader



def printOptions(op):
    print('The following options are used (included default):')
    opkeys = sorted(list(op.keys()))
    for opk in opkeys:
        print(''.join([('%20s' % opk), '\t(', str(op[opk]), ');']))
    print('')

def detect(margs):
    prg_options = {}
    p_options = {};
    prg_options['feature_num'] = 50;
    prg_options['nb_size'] = 3; # 3 adjacent repeat copies--->RGB channel;
    prg_options['label_size'] = 4;
    prg_options['outlog'] = myHeader.M_DEBUG
    
    errorStr = "";
    prg_options["bam"] = margs.bam;
    if prg_options["bam"]==None or prg_options["bam"]=="":
       errorStr += "\tNo bam file (specified by --bam).\n"
    prg_options["outputfolder"] = margs.o;
    if prg_options["outputfolder"]==None or prg_options["outputfolder"]=="":
       errorStr += "\tNo output folder (specified by --o).\n"
    else:
       if not os.path.isdir(prg_options["outputfolder"]):
          os.system('mkdir -p {}'.format(prg_options["outputfolder"]))
    prg_options["f5folder"] = margs.f5folder
    if prg_options["f5folder"]==None or prg_options["f5folder"]=="":
       errorStr += "\tNo f5folder (specified by --f5folder).\n"
    prg_options["f5i"] = margs.f5i;
    if prg_options["f5i"]==None or prg_options["f5i"]=="":
       if not (prg_options["f5folder"]==None or prg_options["f5folder"]==""):
          prg_options["f5i"] = prg_options["f5folder"]+("" if len(prg_options["f5folder"])>0 and prg_options["f5folder"][-1] in ['/'] else "/")+"f5.f5index"
       else:
          prg_options["f5i"] = "f5.f5index"
       prg_options["f5i_basefile"] = "f5"
    else:
       prg_options["f5i_basefile"] = prg_options["f5i"].split("/")[-1].split('.f5index')[0]
       if prg_options["f5i_basefile"]=="":
          prg_options["f5i_basefile"] = "f5"
    prg_options["repeat"] = margs.repeat; # chr4:3074876-3074933:CAG:3
    if prg_options["repeat"] ==None or prg_options["repeat"] =="":
       errorStr += "\tNo repeat information (specified by --repeat).\n"
    else:
       info_sp = prg_options["repeat"].split(":")
       if len(info_sp)<4: errorStr += "\tRepeat format `"+prg_options["repeat"]+"` is not correct. Example: `chr4:3074876-3074933:CAG:3`.\n"
       else:
           p_options['repeat_chr'] = info_sp[0]
           #p_options['repeat_unit_size'] = int(info_sp[3]);
           p_options['repeat_unit_size'] = len(info_sp[2]);
           prg_options['repeat_pat'] = info_sp[2]
           prg_options["merge_gap"] = p_options['repeat_unit_size']*1.5
           pos_2 = info_sp[1].split('-')
           if len(pos_2)==2:
               p_options['repeat_start_pos'] = int(pos_2[0]);
               p_options['repeat_end_pos'] = int(pos_2[1]);
           else:
               errorStr += "\tRepeat format `"+prg_options["repeat"]+"` is not correct for positions. Example: `chr4:3074876-3074933:CAG:3`.\n"
    prg_options["repeat_name"] = margs.repeatName;
    if prg_options["repeat_name"] ==None or prg_options["repeat_name"] =="":
       prg_options["repeat_name"] = "m_repeat"
    prg_options["rpg"] = margs.rpg;
    if prg_options["rpg"]==None or prg_options["rpg"]=="":
       prg_options["rpg"] = "data/trf.v0.bed"
    if not os.path.isfile(prg_options["rpg"]):
       if prg_options["rpg"] == "data/trf.v0.bed":
          prg_options["rpg"] = '/'.join(os.path.abspath(sys.argv[0]).split("/")[:-1]) +"/data/trf.v0.bed"
    if not os.path.isfile( prg_options["rpg"] ):
       errorStr += "\t "+prg_options["rpg"]+" does not exist. Should specified by --rpg.\n"

    prg_options["f5config"] = margs.f5config;
    if prg_options["f5config"] ==None or prg_options["f5config"] =="":
       prg_options["f5config"] = "data/config/fast5_path.config"
    if not os.path.isfile( prg_options["f5config"] ):
       if prg_options["f5config"] == "data/config/fast5_path.config":
          prg_options["f5config"] = '/'.join(os.path.abspath(sys.argv[0]).split("/")[:-1]) +"/data/config/fast5_path.config"
    if not os.path.isfile( prg_options["f5config"] ):
       errorStr += "\t "+prg_options["f5config"]+" does not exist. Should specified by --f5config.\n"
    prg_options["nbsize"] = margs.nbsize;
    prg_options["pcr"] = margs.is_pcr;
    prg_options["UniqueID"] = margs.UniqueID;
    if prg_options["UniqueID"]==None or prg_options["UniqueID"]=="":
       prg_options["UniqueID"] = "_mpred_"

    prg_options["basecalled_path"] = margs.basecalled_path
    prg_options["summary_file"] = margs.summary_file
    if prg_options["summary_file"]==None or prg_options["summary_file"]=="":
       if not (prg_options["f5folder"]==None or prg_options["f5folder"]==""):
          prg_options["summary_file"] = prg_options["f5folder"] + ("" if len(prg_options["f5folder"])>0 and prg_options["f5folder"][-1] in ['/'] else "/")+"sequencing_summary.txt"
       else: prg_options["summary_file"] = "sequencing_summary.txt"
    if (not os.path.isfile( prg_options["summary_file"])) and (not os.path.isfile(prg_options["f5i"])):
       errorStr += "\t "+prg_options["summary_file"]+" does not exist. Should specified by --summary_file.\n"
    prg_options["multif5"] = margs.multif5
    if not prg_options["multif5"]: prg_options["multif5"] = 0;

    prg_options['mod_version'] = 2
    prg_options['mod_path'] = margs.mod_path

    printOptions(prg_options);
    if not errorStr=="":
        print("Error happen!")
        print(errorStr)
        parser.print_help()
        parser.parse_args(['Detect', '--help'])
        #sys.exit(140)

    #print ("{} {}".format(prg_options["f5i"], os.path.isfile(prg_options["f5i"]) ))
    if (not os.path.isfile(prg_options["f5i"])) or os.path.getsize(prg_options["f5i"])<1000:
       _index_cmd = "{}/{} {} {} {} {} {}".format('/'.join(os.path.abspath(sys.argv[0]).split('/')[:-1]), "scripts/IndexF5files", prg_options["f5folder"], prg_options["basecalled_path"],  prg_options["f5i_basefile"], prg_options["summary_file"], prg_options["multif5"])
       print("Generating f5index: {}".format(_index_cmd))
       sys.stdout.flush()
       os.system(_index_cmd)
    if not os.path.isfile(prg_options["f5i"]):
       print("Error!! Cannot generate index file for {}.".format(prg_options["f5i"]))
    else:
       _fs_list_op = ["{}/{} ".format('/'.join(os.path.abspath(sys.argv[0]).split('/')[:-1]), "scripts/genomic1FE")]
       _fs_list_op.append(prg_options["bam"])
       _fs_list_op.append(prg_options["outputfolder"]+"/"+prg_options["UniqueID"]+".fs");
       _fs_list_op.append(prg_options["f5folder"]+"/")
       _fs_list_op.append(prg_options["f5i"])
       _fs_list_op.append(prg_options["repeat"])
       _fs_list_op.append(str(prg_options["rpg"]))
       _fs_list_op.append(prg_options["f5config"])
       _fs_list_op.append("{:.1f}".format(prg_options["nbsize"]))
       _fs_list_op.append("{}".format(margs.algLen))
       _fs_cmd = ' '.join(_fs_list_op);
       print("Generating features: {}".format(_fs_cmd))
       sys.stdout.flush()
       os.system(_fs_cmd)
       
       if not os.path.isfile(prg_options["outputfolder"]+"/"+prg_options["UniqueID"]+".fs"):
          print("Error! Cannot generate fs file: {}".format( prg_options["outputfolder"]+"/"+prg_options["UniqueID"]+".fs" ))
       else: ##  Repeat_CNN_Prediction_motif.gen2.py
          p_options['valid_repeat'] = 1;
          p_options['c_fs_file_list'] = [prg_options["outputfolder"]+"/"+prg_options["UniqueID"]+".fs"]
          p_options['repeat_relax_bp'] = margs.repeat_relax_bp;

          prg_options['filtern'] = 32;
          prg_options['epchon'] = margs.epchon; #200;
          _m_rep = Repeat_CNN_Prediction_model.get_rep_group(prg_options['repeat_pat'])
          prg_options['gn'] = margs.gn; 'hx1-na12878'
          prg_options['mod_version'] = '0.2'
          if prg_options['mod_path']==None:
             #prg_options['mod_path_base'] = "trainedmod_{}_{}/hg_group/{}/{}_{}/".format(prg_options['filtern'], prg_options['mod_version'] , prg_options['gn'], _m_rep, len(_m_rep));
             prg_options['mod_path_base'] = "trainedmod_{}_{}/{}/{}_{}/".format(prg_options['filtern'], prg_options['mod_version'] , prg_options['gn'], _m_rep if len(_m_rep)<6 else _m_rep[:5], len(_m_rep));
             prg_options['mod_path'] = prg_options['mod_path_base']+str(prg_options['epchon'])+'_v'+prg_options['mod_version']+'motif'+str(prg_options['filtern'])+'/np_dl_'+str(p_options['repeat_unit_size']);
          else:
             prg_options['mod_path_base'] = '/'.join(prg_options['mod_path'].split("/")[:-2])

          print("Modused: {}".format(prg_options['mod_path']))
          _m_rp_c = [int((p_options['repeat_end_pos'] - p_options['repeat_start_pos'])/p_options['repeat_unit_size'] +0.5), int((p_options['repeat_end_pos'] - p_options['repeat_start_pos'])/p_options['repeat_unit_size'] +0.5)]
          Repeat_CNN_Prediction_motif_gen2.m_pred(p_options, prg_options, [_m_rp_c])

          if margs.TempRem:
             print("{}".format(margs.TempRem))
             os.system("rm "+prg_options["outputfolder"]+"/"+prg_options["UniqueID"]+".fs")
          


######################################################

parser = argparse.ArgumentParser(description="Determine microsatellite repeat of interests using Nanopore signals.", epilog="For example, \n \
\tpython %(prog)s Detect: detect repeat counts of locus of interest", formatter_class=RawTextHelpFormatter);

subparsers = parser.add_subparsers()
parent_parser = argparse.ArgumentParser(add_help=False)

com_group_arg = parent_parser.add_argument_group('Common options for DeepRepeat')
com_group_arg.add_argument("--f5folder", type=str, 
                                 help="The folder path where is the f5 files for Nanopore data. ");

parser_det = subparsers.add_parser('Detect', parents=[parent_parser], help="Detect repeat counts of locus of interest", description="Detect repeat counts of locus of interest", epilog="For example, \n \
python %(prog)s --f5folder f5f/ ; \n \
", formatter_class=RawTextHelpFormatter)
parser_det.add_argument("--bam", type=str,
                         help="The bam file. Compulsory. ");
parser_det.add_argument("--o", type=str,
                         help="The output folder. Compulsory. ");
parser_det.add_argument("--f5i", type=str,
                         help="The index of f5 files. If empty, `f5.f5index` under `f5folder` will be used. If not exist, `f5.f5index` will be created. Must end with \"f5index\".");
parser_det.add_argument("--repeat", type=str,
                         help="The locus of repeat of interest. For example: \"chr4:3074876-3074933:CAG:3\". Compulsory.");
parser_det.add_argument("--repeatName", type=str,
                         help="The name of the repeat. Default:m_repeat.");
parser_det.add_argument("--rpg", type=str,
                         help="The bed file for neighborhood repeat. If not provided, `data/trf.v0.bed` will be used.");
parser_det.add_argument("--f5config", type=str,
                         help="The config file of basecalling. Default: data/config/fast5_path.config");
parser_det.add_argument("--nbsize", type=float, default=-1.5,
                         help="The size for neighborhood(for feature generation). Minus values for percentage, while plus values for absolute size.");
parser_det.add_argument("--is_pcr", type=bool, default=False,
                         help="Is the data is PCR-based. Default: False.");
parser_det.add_argument("--UniqueID", type=str, default="_mpred_",
                         help="Unique name for the results. Default: \"_mpred_\".");
parser_det.add_argument("--basecalled_path", type=str, default="workspace/pass/",
                         help="The basecalled f5 files under the data folder. Default: \"workspace/pass/\".");
parser_det.add_argument("--summary_file", type=str, #Default="sequencing_summary.txt",
                         help="The summary file of basecalling under the data folder. Default: \"sequencing_summary.txt\" under data folder.");
parser_det.add_argument("--multif5", type=bool, default=False,
                         help="Is multi-fast5 format. Default: False (single-fast5 format).");
parser_det.add_argument("--mod_path", type=str,
                         help="The path of well-trained models. If not provided, it will be automatically generated according to repeat patterns. You need provide the path if you have your own trained models.");
parser_det.add_argument("--epchon", type=int, default=200,
                         help="The well-trained models with epcho. Default: 200.");
parser_det.add_argument("--repeat_relax_bp", type=int, default=15,
                         help="Relaxed nb size (for prediction). Default: 15.");
parser_det.add_argument("--TempRem", type=bool, default=True,
                         help="Will Temp files be removed? Default: True (yes).");
parser_det.add_argument("--gn", type=str, default='hx1-na12878',
                         help="Well-trained models for a genome('hx1', 'na12878' or 'hx1-na12878'). Default: 'hx1-na12878'.");
parser_det.add_argument("--algLen", type=int, default=100,
                         help="Minimal alignment length. Default: 100.");


parser_det.set_defaults(func=detect)



if len(sys.argv) < 2:
   parser.print_help()
else:
   args = parser.parse_args()
   args.func(args)


