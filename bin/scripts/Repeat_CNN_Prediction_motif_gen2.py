
import os, sys

from collections import defaultdict;
import numpy as np
import time;
import copy
import glob

from . import myHeader
from . import myGaussianMixtureModel

from . import Repeat_CNN_Train_genome_motif2 as Repeat_CNN_Train

import tensorflow as tf

from . import Repeat_CNN_Prediction_model

def m_pred(p_options, prg_options, for_test_gt=None):
   np.random.seed(7)

   m_pmd = Repeat_CNN_Prediction_model.Pred_Mod(prg_options['mod_path'])

   #for_test = [[18,41], [26,46]]
   for c_fs_file_id in range(len(p_options['c_fs_file_list'])):
      c_fs_file = p_options['c_fs_file_list'][c_fs_file_id]
      _this_d = Repeat_CNN_Train.read_fs(c_fs_file, p_options, prg_options, 1)
      if _this_d == None: continue;
      c_ty, c_tx, c_pt2 = _this_d
      id_str_info = Repeat_CNN_Prediction_model.read_strand_read_info(c_fs_file)
      c_pt = c_pt2[0]
      c_l_pt = c_pt2[1]
      m_pmd.pred(c_tx);
      m_pmd.sum_for_count(c_pt, c_l_pt, c_ty, id_str_info, p_options, prg_options, for_test_gt[c_fs_file_id] if not for_test_gt==None else None );
   
   print("Testing Finished!")

if __name__=='__main__':
   gn = sys.argv[2];
   filtern = int(sys.argv[3]);
   epchon = int(sys.argv[4]);

   prg_options = {}
   prg_options['feature_num'] = 50;
   prg_options['nb_size'] = 3;
   prg_options['label_size'] = 4;
   prg_options['UniqueID'] = '_pred_'+gn+"Mpred"+sys.argv[11]
   prg_options['merge_gap'] = len(sys.argv[6])*1.5;
   prg_options['repeat_name'] = sys.argv[11]
   prg_options['repeat_pat'] = sys.argv[6]
   prg_options['mod_version'] = sys.argv[5]
   _m_rep = Repeat_CNN_Prediction_model.get_rep_group(sys.argv[6])
   if gn.lower() not in ['hx1', 'na12878', 'hx1_na12878', 'na12878_hx1', 'hx1-na12878', 'na12878-hx1']:
      prg_options['mod_path_base'] = 'sep_mod_output_v'+prg_options['mod_version']+'motif'+str(filtern)+"/sep/"+gn+"/"
   else:
      if not sys.argv[12] == 'general':
         prg_options['mod_path_base'] = 'repeat_mod_output_v'+prg_options['mod_version']+'motif'+str(filtern)+"/"+gn+"/"+_m_rep+"_"+str(sys.argv[7])+"/"
      else: 
         prg_options['mod_path_base'] = 'repeat_mod_output_v'+prg_options['mod_version']+'motif'+str(filtern)+"/"+gn+"/"+str(sys.argv[7])+"/"
   prg_options['mod_path'] = prg_options['mod_path_base']+str(epchon)+'_v'+prg_options['mod_version']+'motif'+str(filtern)+'/np_dl_'+str(sys.argv[7])

   if sys.argv[12] == 'pcr':
      prg_options['pcr'] = 1
   prg_options['outlog'] = myHeader.M_DEBUG

   p_options = {}
   p_options['repeat_unit_size'] = int(sys.argv[7]);
   p_options['repeat_chr'] = sys.argv[8]
   p_options['repeat_start_pos'] = int(sys.argv[9])
   p_options['repeat_end_pos'] = int(sys.argv[10])
   p_options['repeat_relax_bp'] = int(sys.argv[14]) ### can be used when repeat motif is trained separately.
   p_options['valid_repeat'] = int(sys.argv[15]) # 1: the repeat above is used, otherwise: not.
   p_options['c_fs_file_list'] = ['nanopore_data/fs/'+sys.argv[12]+'/motif_test/'+sys.argv[6]+'_'+sys.argv[11]+'_'+str(sys.argv[7])+'.fs']
   p_options['c_fs_file_list'] = [sys.argv[16]];
   print ("Test: {}".format(p_options['c_fs_file_list']))
   _m_rp_c = sys.argv[13].split(',')
   _m_rp_c[0] = int(_m_rp_c[0])
   _m_rp_c[1] = int(_m_rp_c[1])
   m_pred(p_options, prg_options, [_m_rp_c])


