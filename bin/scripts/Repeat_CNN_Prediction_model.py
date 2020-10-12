
import os, sys

from collections import defaultdict;
import numpy as np
import time;
import copy
import glob

from . import myHeader
from . import myGaussianMixtureModel

import tensorflow as tf

class Pred_Mod:
   def __init__(self, _mod_path, _batch_size = 1024):
      self.batch_size = _batch_size;
      self.mod_path = _mod_path;

      config = tf.ConfigProto()
      config.gpu_options.allow_growth = True
      
      self.sess = tf.Session(config=config)
      self.new_saver = tf.train.import_meta_graph(_mod_path+'.meta')
      self.new_saver.restore(self.sess, tf.train.latest_checkpoint(_mod_path[:(_mod_path.rfind('/')+1)]));

      #for gv in tf.global_variables():
      #   print(gv)
      #print (tf.global_variables())
    
      self.mgraph = tf.get_default_graph()
      ######
      #for n in self.mgraph.as_graph_def().node:
      #   print (n.name)
      #for op in self.mgraph.get_operations():
      #   print(op)
      ##########
      self.X = self.mgraph.get_tensor_by_name('X:0')
      self.dropout = self.mgraph.get_tensor_by_name('dropout:0')
      self.pred_labels = self.mgraph.get_tensor_by_name('pred_labels:0')
      self.learning_rate = self.mgraph.get_tensor_by_name('learning_rate:0') 
 
   def pred(self, _p_tx):
      self.pred_list = []
      c_group_size = int(len(_p_tx)/self.batch_size)
      if c_group_size<1: c_group_size = 1
      c_x_list = np.array_split(_p_tx, c_group_size)
      for c_g_i in range(c_group_size):
          moutp = self.sess.run([self.pred_labels], feed_dict={self.X:c_x_list[c_g_i],self.dropout:0,self.learning_rate:0});
          if((c_g_i+1)%100==0): print("moutp = {}/{} {} {}".format(len(moutp), len(moutp[0]), moutp[:10], len(self.pred_list) ))
          if len(self.pred_list)==0: self.pred_list = moutp[0]
          else: self.pred_list = np.concatenate((self.pred_list, moutp[0]));

   def sum_for_count(self, p_pt, p_l_pt, p_ty, id_str_info, p_options, prg_options, gt=None):
      print ("sum_for_count={}/{}/{} {} {}; {}".format(len(self.pred_list), len(p_ty), len(p_l_pt), self.pred_list[:10], self.pred_list.shape, len(id_str_info)));
      print ("{} {} {} {}".format(p_pt[:5], p_pt[-5:], len(id_str_info), len(p_pt)))
      self.repeat_count_list = defaultdict()
      c_pl_i = 0;
      c_pt_i = 1;

      _tp, _fp, _fn, _tn = 0, 0, 0, 0
      _c_tp, _c_fp, _c_fn, _c_tn = 0, 0, 0, 0
      _perf_strand = [[], []]
      for _i_ in range(len(p_l_pt)):
         if (c_pt_i<len(p_pt) and _i_<p_pt[c_pt_i]) or (c_pt_i>=len(p_pt) and _i_<len(p_l_pt)):
            if p_l_pt[_i_][1] in [1,2]:
               if p_ty[_i_][0]==1:
                  if self.pred_list[_i_]==0: _c_tn += 1;
                  else: _c_fp += 1;
               else:
                  if self.pred_list[_i_]==0: _c_fn += 1;
                  else: _c_tp += 1;
         else:
            _c_p = _c_tp/(float(_c_tp + _c_fp) if _c_tp + _c_fp>0 else 1);
            _c_r = _c_tp/(float(_c_tp + _c_fn) if _c_tp + _c_fn>0 else 1);
            _c_f1= 2*_c_p*_c_r/((_c_p+_c_r) if _c_p+_c_r>0 else 1)
            print ("{}.{} tp={} fp={} fn={} tn={} p={:.3f} r={:.3f} f1={:.3f}".format(c_pt_i,'/'.join(id_str_info[c_pt_i-1]) if c_pt_i-1<len(id_str_info) else "None", _c_tp, _c_fp, _c_fn, _c_tn, _c_p,_c_r,_c_f1))
            if c_pt_i-1<len(id_str_info):
               _perf_strand[int(id_str_info[c_pt_i-1][1])].append(_c_f1)
            _c_tp, _c_fp, _c_fn, _c_tn = 0, 0, 0, 0
            c_pt_i += 1;
         if p_l_pt[_i_][1] in [1,2]:
            if p_ty[_i_][0]==1:
               if self.pred_list[_i_]==0: _tn += 1;
               else: _fp += 1;
            else:
               if self.pred_list[_i_]==0: _fn += 1;
               else: _tp += 1;
            continue;
            #if p_ty[_i_][self.pred_list[_i_]] >0 :
            #   if self.pred_list[_i_]==0: _tn += 1;
            #   else: _tp += 1;
            #else:
            #   if self.pred_list[_i_]==0: _fp += 1;
            #   else: _fn += 1;
      _c_p = _c_tp/(float(_c_tp + _c_fp) if _c_tp + _c_fp>0 else 1);
      _c_r = _c_tp/(float(_c_tp + _c_fn) if _c_tp + _c_fn>0 else 1);
      _c_f1= 2*_c_p*_c_r/((_c_p+_c_r) if _c_p+_c_r>0 else 1)
      print ("{}.{} tp={} fp={} fn={} tn={} p={:.3f} r={:.3f} f1={:.3f}\n".format(c_pt_i, '/'.join(id_str_info[c_pt_i-1]) if c_pt_i-1<len(id_str_info) else "None", _c_tp, _c_fp, _c_fn, _c_tn, _c_p,_c_r,_c_f1))

      print("Strand-perf: +:{:.3f}/{:.3f}/{:,}; -:{:.3f}/{:.3f}/{:,}\n".format(np.average(_perf_strand[0]),np.std(_perf_strand[0]), len(_perf_strand[0]), np.average(_perf_strand[1]),np.std(_perf_strand[1]), len(_perf_strand[1])))

      _p = _tp/(float(_tp + _fp) if _tp + _fp>0 else 1);
      _r = _tp/(float(_tp + _fn) if _tp + _fn>0 else 1);
      _f1= 2*_p*_r/((_p+_r) if _p+_r>0 else 1)
      print ("tp={} fp={} fn={} tn={} p={:.3f} r={:.3f} f1={:.3f}\n".format(_tp, _fp, _fn, _tn, _p,_r,_f1))

      if p_options['valid_repeat']==0: return;

      c_pl_i = 0;
      c_pt_i = 1;
      fwd_rev_num = [[0,0], [0,0]] 
      while c_pl_i < len(self.pred_list):
         #if c_pt_i > len(p_pt)-5:
         #   print ('p_pt= {} {} {}'.format(c_pl_i, c_pt_i, (p_pt[c_pt_i] if c_pt_i<len(p_pt) else -1)))
         if c_pt_i<len(p_pt):
            c_read_pred = self.pred_list[c_pl_i:p_pt[c_pt_i]]
         else:
            c_read_pred = self.pred_list[c_pl_i:]

         try:
            fwd_rev_num[0][int(id_str_info[c_pt_i-1][1])] += 1;
         except:
            print ("Error of prediction: {} {}".format(c_pt_i, len(id_str_info)));
            print ("Error of prediction: {} {} {}".format(c_pt_i, id_str_info[c_pt_i-1][1], id_str_info[c_pt_i-1]))

         rep_list = [[]]
         rep_list_pos = [[]]
         add_del = 0;
         for _rp_i in range(len(c_read_pred)):
            #if c_read_pred[_rp_i] in [1,2,3]:
            if c_read_pred[_rp_i] in [1,2,3] and ((p_options['repeat_start_pos']-p_options['repeat_relax_bp'])<p_l_pt[c_pl_i+_rp_i][0]<(p_options['repeat_end_pos']+p_options['repeat_relax_bp'])):
               if (len(rep_list[-1])>0 and (_rp_i-rep_list[-1][-1])>prg_options['merge_gap']):
                  rep_list.append([])
                  rep_list_pos.append([])
               if (c_read_pred[_rp_i]==3):
                  rep_list[-1].append(_rp_i)
                  rep_list_pos[-1].append(p_l_pt[c_pl_i+_rp_i][0])
               elif (c_read_pred[_rp_i]==2): pass;
               else: 
                  add_del += 1;
                  rep_list[-1].append(_rp_i)
                  rep_list[-1].append(_rp_i)
                  rep_list_pos[-1].append(p_l_pt[c_pl_i+_rp_i][0])
                  rep_list_pos[-1].append(p_l_pt[c_pl_i+_rp_i][0])

         #rep_list = rep_list_pos
         #if int(id_str_info[c_pt_i-1][1]) == 1:
         #   for _rev_i in range(len(rep_list)):
         #      rep_list[_rev_i] = rep_list[_rev_i][::-1]
         #   rep_list = rep_list[::-1]
                     
         #print (rep_list) 
         old_rep_list_pos = copy.deepcopy(rep_list_pos)
         while True:
            #if len(rep_list)>0 and len(rep_list[0])<p_options['repeat_unit_size']*2: # htt;
            if len(rep_list)>0 and len(rep_list[0])<p_options['repeat_unit_size']*1.5: # 2htt;
            #if len(rep_list)>0 and len(rep_list[0])<p_options['repeat_unit_size']*2:
               rep_list = rep_list[1:]
               rep_list_pos = rep_list_pos[1:]
            else: break;
         while True:
            #if len(rep_list)>0 and len(rep_list[-1])<p_options['repeat_unit_size']*2: # htt;
            if len(rep_list)>0 and len(rep_list[-1])<p_options['repeat_unit_size']*1.5: # 2htt;
            #if len(rep_list)>0 and len(rep_list[-1])<p_options['repeat_unit_size']*2:
               rep_list = rep_list[:-1]
               rep_list_pos = rep_list_pos[:-1]
            else: break;

         old_rep_list = rep_list;
         if len(rep_list)>0:
            non_e_i = 0;
            while non_e_i<len(rep_list) and len(rep_list[non_e_i])==0:
                non_e_i += 1
            if non_e_i<len(rep_list): rep_n_list_ = [copy.deepcopy(rep_list[non_e_i])]
            else: rep_n_list_ = []
            for _rl_i in range(non_e_i+1, len(rep_list)):
               if len(rep_list[_rl_i])==0: continue;
               if ( (rep_list[_rl_i][0] - rep_n_list_[-1][-1])*3 <= ( len(rep_list[_rl_i]) if len(rep_list[_rl_i]) < len(rep_n_list_[-1]) else len(rep_n_list_[-1]) )):
                  rep_n_list_[-1].extend(rep_list[_rl_i]) 
               else: rep_n_list_.append(copy.deepcopy(rep_list[_rl_i]))
            rep_list = rep_n_list_;
        
         _reg_len = 0;
         if len(rep_list)>0:
            _reg_len = int(( rep_list[-1][-1] - rep_list[0][0] )/p_options['repeat_unit_size']+0.5 )
 
         rep_c = [0]; 
         len_sum = 0;
         for _rl in rep_list:
            rep_c.append(int(len(_rl)/p_options['repeat_unit_size']+0.5))
            ##### sum all repeats
            rep_c[0] += rep_c[-1]; 
            if rep_c[-1]>2:
               len_sum += rep_c[-1]

         rep_c = sorted(rep_c)
         if True: #rep_c[-1]>0:
            _reg_len = 0;
            if len(rep_list)>0:
               _reg_len = int(( rep_list[-1][-1] - rep_list[0][0] )/p_options['repeat_unit_size']+0.5 )
            rep_c.append(_reg_len+2)

            if rep_c[-1] not in self.repeat_count_list:
               self.repeat_count_list[rep_c[-1]] = 1
            else: self.repeat_count_list[rep_c[-1]] += 1
            if not gt==None:
               if (rep_c[-1] not in range(gt[0]-4, gt[0]+5)) and (rep_c[-1] not in range(gt[1]-4, gt[1]+5)):
                  try:
                     fwd_rev_num[1][int(id_str_info[c_pt_i-1][1])] += 1;
                  except:
                     print ("Error of prediction: {} {} {}".format(c_pt_i, id_str_info[c_pt_i-1][1], id_str_info[c_pt_i-1])) 
                  if rep_c[-1] < (gt[0]+gt[1])/2:
                     print ("RX1:{}/{}/{}/{} {} {} >>>> {} <<<< {}".format(rep_c[-1], rep_c[-2] if len(rep_c)>1 else -1, len_sum, add_del, id_str_info[c_pt_i-1], old_rep_list, rep_list_pos, old_rep_list_pos))
                  else:
                     print ("RX2:{}/{}/{}/{} {} {} >>>> {} <<<< {}".format(rep_c[-1], rep_c[-2] if len(rep_c)>1 else -1, len_sum, add_del, id_str_info[c_pt_i-1], old_rep_list, rep_list_pos, old_rep_list_pos))
               else:
                  pass; #print ("RC:{}/{}/{} {} {}".format(rep_c[-1], rep_c[-2] if len(rep_c)>1 else -1, len_sum, id_str_info[c_pt_i-1], old_rep_list))
            else:
               if rep_c[-1] < 12:
                  print ("R1:{}/{}/{} {} {} >>>> {}".format(rep_c[-1], rep_c[-2] if len(rep_c)>1 else -1, len_sum, id_str_info[c_pt_i-1], old_rep_list, rep_list_pos))
               elif 31 <rep_c[-1] < 36:
                  print ("R2:{}/{}/{} {} {} >>>> {}".format(rep_c[-1], rep_c[-2] if len(rep_c)>1 else -1, len_sum, id_str_info[c_pt_i-1], old_rep_list, rep_list_pos))
               else: print ("Rx:{}/{}/{} {} {} >>>> {}".format(rep_c[-1], rep_c[-2] if len(rep_c)>1 else -1, len_sum, id_str_info[c_pt_i-1], old_rep_list, rep_list_pos))

         if c_pt_i<len(p_pt):
            c_pl_i = p_pt[c_pt_i]
            c_pt_i += 1
         else: c_pl_i = len(self.pred_list)

      c_rckeys = sorted(list(self.repeat_count_list.keys()))
      self.c_rk_list = []
      for _rck in c_rckeys:
         self.c_rk_list.append("{}:{}".format(_rck, self.repeat_count_list[_rck]))
      print("Repeat_count distribution= "+' '.join(self.c_rk_list))
      print("fwd_rev_num(0:fwd;1:rev)={}".format(fwd_rev_num))

      self._gmm_res = myGaussianMixtureModel.myGMM(self.repeat_count_list, 5, commonoptions=prg_options)
      print("GMM: {}-{} {} {}".format(prg_options['repeat_name'], prg_options['repeat_pat'], self._gmm_res[0], self._gmm_res[1]))
 
   def __del__(self):
      self.sess.close();

def read_strand_read_info(fmore_file):
   _gz_str = '.gz'
   m_id_str_info = np.loadtxt(fmore_file+'.more' if not fmore_file[-len(_gz_str):]==_gz_str else fmore_file[:-len(_gz_str)]+'.more', dtype={'names': ('chr', 'strand', 'read_id'), 'formats': ('U100', 'U1', 'U300')});
   return m_id_str_info;

def getComplementary(bp, na):
    if na.upper() not in bp:
        print(na, na.upper(), bp)
    return bp[na.upper()]

def getComplementary3(bp, na3):
    c3 = ''
    for li in range(len(na3)):
        c3 += getComplementary(bp, na3[li])
    return (c3[::-1])

def get_rep_group(rep_u_k):
   rep_u_k_list = [rep_u_k];
   for _ap_i in range(1, len(rep_u_k)):
      rep_u_k_list.append(rep_u_k_list[-1][1:]+rep_u_k_list[-1][0])
   for _ruk in range(len(rep_u_k)):
      rep_u_k_list.append(getComplementary3(myHeader.NA_bp_, rep_u_k_list[_ruk]))
   rep_u_k_list = sorted(rep_u_k_list)
   return rep_u_k_list[0]

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
   _m_rep = get_rep_group(sys.argv[6])
   if sys.argv[12] == 'general':
      prg_options['mod_path_base'] = 'repeat_mod_output_v'+prg_options['mod_version']+'motif'+str(filtern)+"/"+gn+"/"+((str(len(_m_rep))+"/") if not len(_m_rep)==3 else "")
   else: 
      prg_options['mod_path_base'] = 'repeat_mod_output_v'+prg_options['mod_version']+'motif'+str(filtern)+"/"+gn+"/"+_m_rep+"/"
   prg_options['mod_path'] = prg_options['mod_path_base']+str(epchon)+'_v'+prg_options['mod_version']+'motif'+str(filtern)+'/np_dl_'+str(len(_m_rep))

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
   #m_pred(p_options, prg_options, [_m_rp_c])


