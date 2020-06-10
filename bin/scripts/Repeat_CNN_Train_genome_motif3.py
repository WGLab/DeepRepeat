
import os, sys
import Repeat_CNN_model_motif3 as Repeat_CNN_model
from collections import defaultdict;
import numpy as np
import time;

import glob
import gc

def m_reshape(prg_options, p_options, traing_or_testing):
   n50_size = prg_options['nb_size'] * p_options['repeat_unit_size']
   is_odd = n50_size % 2
   lr_size = int(n50_size/2);
   
   n0_list = [];
   for n50_i in range(len(p_options['ty'])):
      if p_options['p_t'][n50_i][0]==0:
          n0_list.append(n50_i);

   formatted_x = []
   formatted_y = []
   if traing_or_testing==1: 
      formatted_pt = [0];
      formatted_l_pt = []
   n50_i = 1;
   next_n0 = 1;
   # print ("In m_reshape: {} {}".format(n50_size, lr_size))
   # print (n0_list[:20])
   while n50_i < len(p_options['ty'])-n50_size:
      if n50_i+n50_size > n0_list[next_n0]:
         # print ("B: {} {} {}/{}".format(n50_i, n50_size, next_n0, n0_list[next_n0]))
         n50_i = n0_list[next_n0]
         next_n0 += 1
         # print ("B: {} {} {}/{}".format(n50_i, n50_size, next_n0, n0_list[next_n0]))
         # if (next_n0)>2: break;
         if traing_or_testing==1:
            formatted_pt.append(len(formatted_y))
      else:
         if (p_options['p_t'][n50_i+lr_size][1] in [1,2]) or traing_or_testing==1:
            formatted_x.append(np.transpose([p_options['tx'][(n50_i+nbi*p_options['repeat_unit_size']):(n50_i+(nbi+1)*p_options['repeat_unit_size'])] for nbi in range(prg_options['nb_size'])] ) )
            #print ("F: {} {} {}".format(p_options['p_t'][n50_i+lr_size][0], p_options['p_t'][n50_i+lr_size][1], p_options['ty'][n50_i+lr_size]))
            #for nbi in range(prg_options['nb_size']):
            #   for rui in range(p_options['repeat_unit_size']):
            #      print (' '.join([str(fv) for fv in formatted_x[-1][:,rui,nbi].astype(int)]))
            #   print ()
            ######################################
            #formatted_x.append([p_options['tx'][(n50_i+nbi*p_options['repeat_unit_size']):(n50_i+(nbi+1)*p_options['repeat_unit_size'])] for nbi in range(prg_options['nb_size'])] )
            ######################################
            #formatted_x.append([np.transpose(p_options['tx'][(n50_i+nbi*p_options['repeat_unit_size']):(n50_i+(nbi+1)*p_options['repeat_unit_size'])]) for nbi in range(prg_options['nb_size'])] )
            #print ("F: {} {} {}".format(p_options['p_t'][n50_i+lr_size][0], p_options['p_t'][n50_i+lr_size][1], p_options['ty'][n50_i+lr_size]))
            #for nbi in range(prg_options['nb_size']):
            #   for rui in range(p_options['repeat_unit_size']):
            #      print (' '.join([str(fv) for fv in formatted_x[-1][nbi][:,rui].astype(int)]))
            #   print ()
            #
            formatted_y.append(p_options['ty'][n50_i+lr_size])
            if traing_or_testing==1: 
               formatted_l_pt.append(p_options['p_t'][n50_i+lr_size])
      n50_i += 1;

   print("Size: {} {} {} after {} {}".format(len(p_options['p_t']),len(p_options['ty']),len(p_options['tx']), len(formatted_y), len(formatted_x)  ))

   return (formatted_y, formatted_x, ((formatted_pt, formatted_l_pt) if traing_or_testing==1 else None)) 
      

def read_fs(fs_file, p_options, prg_options, traing_or_testing=0):
   print("\nRead file = {}".format(fs_file));
   sys.stdout.flush()
   c_data = np.loadtxt(fs_file, dtype=np.float32);

   if len(c_data)<1: return None;

   p_t, ty, ts, tx = np.split(c_data, [2, 6, 9], axis = 1);

   p_options['p_t'] = p_t.astype(np.int32);
   p_options['ty'] = ty.astype(int);
   p_options['tx'] = tx

   return m_reshape(prg_options, p_options, traing_or_testing)

def _re_group_data(_x, _y):
   _y1_p = (_y[:,0]==1);
   _y0_p = (_y[:,0]==0);
   _pos_size = np.sum(_y1_p);
   _neg_size = np.sum(_y0_p);

   _third4 = int(_pos_size*0.75+0.5)
   if _third4<64: _third4 = 64;
   #print _pos_size, _third4, _third4, _third4*1.4
   if _neg_size < _third4*1.4:
      return [(_x, _y)]
   else:
      _x_1 = _x[_y1_p, :]
      _y_1 = _y[_y1_p, :];

      _x_0 = _x[_y0_p, :];
      _y_0 = _y[_y0_p, :];

      _n_grp_size = int(_neg_size/_third4+0.7)
      _x_0_list = np.array_split(_x_0, _n_grp_size)
      _y_0_list = np.array_split(_y_0, _n_grp_size)

      _m_ret_data = []
      for _m_i in range(len(_x_0_list)):
         _x_1_c = np.concatenate((_x_1, _x_0_list[_m_i]))
         _y_1_c = np.concatenate((_y_1, _y_0_list[_m_i]))
         _m_ret_data.append((_x_1_c, _y_1_c))
      return _m_ret_data

def m_train(p_options, prg_options):
   np.random.seed(7)

   c_r3cnn = Repeat_CNN_model.Repeat_CNN_model(prg_options['label_size'], prg_options['feature_num'], p_options['repeat_unit_size'], prg_options['nb_size'], prg_options['gpu_cpu'], prg_options['num_filter_in_first']);
   c_r3cnn.prepare_run();
   
   #train_epoches = 100;
   #batch_size = 1024;
   ##batch_size = 512;
   #batch_size = 2048
   #batch_size = 512;
   if prg_options['num_filter_in_first']>=256: batch_size = 256
   elif prg_options['num_filter_in_first']>=128: batch_size = 512;
   else: 
      batch_size = 2048;
      batch_size = 4096;
   batch_size = 2048
   batch_size = 1024
   batch_size = 512;
   if prg_options['num_filter_in_first']>=256: batch_size = 512
   elif prg_options['num_filter_in_first']>=128: batch_size = 1024;
   else: batch_size = 2048;

   if prg_options['gpu_cpu'] in  ['CPU','cpu']: batch_size = 4096;

   batch_size = 512
   if "batch_size" in p_options and p_options["batch_size"]>32:
      batch_size = p_options["batch_size"]

   save_ep = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
   ## save_ep = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
   #if ('motif_of_interest' in prg_options and prg_options['motif_of_interest'] in ['ACG', 'ACT']):
   #   save_ep_n = [_se*10 for _se in save_ep];
   #   save_ep = save_ep_n
   file_id_list = np.arange(len(p_options['c_fs_file_list']))
   ##############################################
   read_10_time = time.time()
   data_10 = []; #[read_fs(c_fs_file, p_options, prg_options) for c_fs_file in p_options['c_fs_file_list']]
   less_sample_n = 0
   for c_fs_file in p_options['c_fs_file_list']:
      _this_d = read_fs(c_fs_file, p_options, prg_options)
      if _this_d==None: continue;
      data_10.append(_this_d)
      if len(data_10[-1][0])<1000000: less_sample_n += 1;
      gc.collect()
   if len(data_10)==0:  
      print("None data for {}".format(p_options['c_fs_file_list']))
      return;
   del p_options['p_t']
   del p_options['ty']
   del p_options['tx']
   gc.collect()
   print("Total IO time: {:.3f} {}".format(time.time()-read_10_time, less_sample_n))
   more_times = 1; #10;
   save_ep_n = [_se*more_times for _se in save_ep];
   save_ep = save_ep_n
   #
   #if less_sample_n >= len(data_10)/2:
   #   save_ep_n = [_se*10 for _se in save_ep];
   #   save_ep = save_ep_n
   #else:
   #   save_ep.append(300)
   #   save_ep.append(500)
   #   #save_ep.append(700)
   #   #save_ep.append(1000)
   train_epoches = save_ep[-1]
   if "max_epoch" in p_options and p_options["max_epoch"]>1:
      train_epoches = p_options["max_epoch"]
      if p_options["max_epoch"] not in save_ep:
         save_ep.append(p_options["max_epoch"])

   format_data = []; m_group_size = []; max_g_id = -1;
   if len(data_10)>0:
      for c_fs_file_ind in range(len(p_options['c_fs_file_list'])):
         m_group_size.append(int(len(data_10[c_fs_file_ind][0])/batch_size))
         if  m_group_size[-1] >  m_group_size[max_g_id] or max_g_id==-1:
             max_g_id = c_fs_file_ind
         format_data.append((np.array_split(data_10[c_fs_file_ind][0], m_group_size[-1]), np.array_split(data_10[c_fs_file_ind][1], m_group_size[-1])))
         _re_shuffle = []
         for _o_i in range(len(format_data[-1][0])):
            _re_shuffle.append((format_data[-1][0][_o_i], format_data[-1][1][_o_i]))
         np.random.shuffle(_re_shuffle);
         format_data[-1] = _re_shuffle

      #print("{} {} {}".format(m_group_size, len(format_data[0][0]), len(format_data[0][1])))
      #for _g_i in range(m_group_size[0]):
      #   print("{} {}".format(len(format_data[0][0][_g_i]), len(format_data[0][1][_g_i])))
      oplist = ["Test group:"]
      for c_fs_file_ind in range(len(p_options['c_fs_file_list'])):
         oplist.append("{}/{}/{}/{};{}/{}\n".format(len(format_data[c_fs_file_ind]), len(format_data[c_fs_file_ind][0][1]), len(format_data[c_fs_file_ind][0][1][0]), len(format_data[c_fs_file_ind][0][1][0][0]), len(format_data[c_fs_file_ind][0][0]), len(format_data[c_fs_file_ind][0][0][0])))
      print(' '.join(oplist))
      sys.stdout.flush()
   del data_10
   gc.collect()
   ################################################

   start_time = time.time();
   for c_e in range(1, train_epoches+1):
      #if c_e<2: this_learning_rate = 0.01;
      #elif c_e<10: this_learning_rate = 0.001;
      #elif c_e<20: this_learning_rate = 0.0001;
      #elif c_e<30: this_learning_rate = 0.00005;
      #else: this_learning_rate = 0.00001;
      #if c_e<2: this_learning_rate = 0.001;
      #elif c_e<5: this_learning_rate = 0.0005;
      #elif c_e<10: this_learning_rate = 0.0001;
      #elif c_e<20: this_learning_rate = 0.00005;
      #else: this_learning_rate = 0.00001;
      if c_e<(p_options['lr']):                 this_learning_rate = 0.005
      elif c_e<(p_options['lr'])*10*more_times:   this_learning_rate = 0.001;
      elif c_e<(p_options['lr'])*100*more_times:  this_learning_rate = 0.0005;
      #elif c_e<100*more_times: this_learning_rate = 0.0001;
      else:                     this_learning_rate = 0.0001;
      
      _this_id_list = np.array([-1 for _ in file_id_list])
      for _gs_i in range(m_group_size[max_g_id]):
         #if c_e==1 and _gs_i<1000: this_learning_rate = 0.005
         #elif c_e==1 and _gs_i>=1000: this_learning_rate = 0.001;
         #if c_e==1 and _gs_i<m_group_size[max_g_id]/2: this_learning_rate = 0.005
         #elif c_e==1 and _gs_i>=m_group_size[max_g_id]/2: this_learning_rate = 0.001;

         for _c_id_ in range(len(file_id_list)):
            _this_id_list[_c_id_] += 1
            if _this_id_list[_c_id_] >= m_group_size[_c_id_]:
               _this_id_list[_c_id_] = 0;
         c_start_time = time.time();
         np.random.shuffle(file_id_list);

         if (_gs_i)%1000==0: 
            print ("==={}=========step========={}/{}/{}/{}/{}".format(c_e, _gs_i, _this_id_list[max_g_id], m_group_size[max_g_id], batch_size, this_learning_rate));
            for c_fs_file_ind in range(len(file_id_list)):
               c_r3cnn.sess.run(c_r3cnn.init_l);
               try:
                  loss, aucm, acc, p, r = c_r3cnn.sess.run([c_r3cnn.loss_op,c_r3cnn.auc_op[1],c_r3cnn.accuracy,c_r3cnn.precision[1],c_r3cnn.recall[1]], feed_dict={c_r3cnn.X:format_data[c_fs_file_ind][_this_id_list[c_fs_file_ind]][1], c_r3cnn.Y:format_data[c_fs_file_ind][_this_id_list[c_fs_file_ind]][0], c_r3cnn.dropout:0, c_r3cnn.learning_rate:0});
                  print(">>>Tratin#files {:5d}/g={:5d}/{}, loss={:.7f} AUC={:.3f} acc={:.3f} prec={:.3f} rec={:.3f}. Time consuming={:.3f}(c={:.3f})".format(_this_id_list[c_fs_file_ind], m_group_size[c_fs_file_ind], c_fs_file_ind, loss,aucm,acc,p,r, time.time()-start_time, time.time()-c_start_time  ))
               except:
                  print(">>>Tratin#filesError {} for {}. Time consuming={:.3f}(c={:.3f}) g={} X-dim={}/{}/{}/{} Y-dim={}/{}".format(_this_id_list[c_fs_file_ind], c_fs_file_ind, time.time()-start_time, time.time()-c_start_time, len(format_data[c_fs_file_ind]), len(format_data[c_fs_file_ind][_this_id_list[c_fs_file_ind]][1]),len(format_data[c_fs_file_ind][_this_id_list[c_fs_file_ind]][1][0]),len(format_data[c_fs_file_ind][_this_id_list[c_fs_file_ind]][1][0][0]),len(format_data[c_fs_file_ind][_this_id_list[c_fs_file_ind]][1][0][0][0]), len(format_data[c_fs_file_ind][_this_id_list[c_fs_file_ind]][0]),len(format_data[c_fs_file_ind][_this_id_list[c_fs_file_ind]][0][0]) )) 
               sys.stdout.flush()

         if len(file_id_list)==1:
            _this_x = format_data[file_id_list[0]][_this_id_list[file_id_list[0]]][1]
            _this_y = format_data[file_id_list[0]][_this_id_list[file_id_list[0]]][0]
         else:
            _this_x = np.concatenate((format_data[file_id_list[0]][_this_id_list[file_id_list[0]]][1], format_data[file_id_list[1]][_this_id_list[file_id_list[1]]][1]))
            _this_y = np.concatenate((format_data[file_id_list[0]][_this_id_list[file_id_list[0]]][0], format_data[file_id_list[1]][_this_id_list[file_id_list[1]]][0]))
            if len(file_id_list)>2:
               for c_fs_file_ind in file_id_list[2:]:
                  _this_x = np.concatenate((_this_x, format_data[file_id_list[c_fs_file_ind]][_this_id_list[file_id_list[c_fs_file_ind]]][1]))
                  _this_y = np.concatenate((_this_y, format_data[file_id_list[c_fs_file_ind]][_this_id_list[file_id_list[c_fs_file_ind]]][0]))
         if c_e==1 and _gs_i<5:
            print("_this: X-dim={}/{}/{}/{} Y-dim={}/{} x={} y={}".format(len(_this_x),len(_this_x[0]),len(_this_x[0][0]),len(_this_x[0][0][0]),  len(_this_y),len(_this_y[0]), c_r3cnn.X.get_shape(), c_r3cnn.Y.get_shape()))
            sys.stdout.flush()
         #trn_op = c_r3cnn.sess.run([c_r3cnn.train_op,c_r3cnn.loss_op], feed_dict={c_r3cnn.X:_this_x, c_r3cnn.Y:_this_y, c_r3cnn.dropout:0.5, c_r3cnn.learning_rate:this_learning_rate})
         _n_g = _re_group_data(_this_x, _this_y)
         for _n_i in range(len(_n_g)):
            trn_op = c_r3cnn.sess.run([c_r3cnn.train_op,c_r3cnn.loss_op], feed_dict={c_r3cnn.X:_n_g[_n_i][0], c_r3cnn.Y:_n_g[_n_i][1], c_r3cnn.dropout:0.5, c_r3cnn.learning_rate:this_learning_rate})
      if c_e in save_ep:
         if (not os.path.isdir(prg_options['save_mod_path']+'/'+str(c_e))):
            os.system('mkdir -p '+prg_options['save_mod_path']+'/'+str(c_e)); time.sleep(3)
         c_r3cnn.saver.save(c_r3cnn.sess, prg_options['save_mod_path']+'/'+str(c_e)+'/'+prg_options['UniqueID']);
         if (not os.path.isdir(''+prg_options['save_mod_path']+'/'+str(c_e)+'_v'+prg_options['mod_version']+'')):
            os.system('mkdir -p '+prg_options['save_mod_path']+'/'+str(c_e)+'_v'+prg_options['mod_version']+'');
         time.sleep(30)
         os.system('cp -r '+prg_options['save_mod_path']+'/'+str(c_e)+' '+prg_options['save_mod_path']+'/'+str(c_e)+'_v'+prg_options['mod_version']+'/')
         

   print("Training Finished!")

if __name__=='__main__':
   prg_options = {}
   prg_options['gpu_cpu'] = 'gpu'
   #prg_options['gpu_cpu'] = 'cpu'
   if len(sys.argv)>5:
      prg_options['gpu_cpu'] = sys.argv[5]
   prg_options['feature_num'] = 50;
   prg_options['nb_size'] = 3;
   prg_options['label_size'] = 4;
   prg_options['num_filter_in_first'] = int(sys.argv[2])
   if len(sys.argv)>4:
      if sys.argv[4].lower() not in ['all']:
         prg_options['motif_of_interest'] = sys.argv[4]
   prg_options['mod_version'] = '0.3motif'+str(prg_options['num_filter_in_first']); # flexible filters for last cnn layers
   prg_options['g_name'] = ['hx1']
   prg_options['g_name'] = sys.argv[1].split(',')
   #prg_options['save_mod_path'] = 'repeat_mod_output_v'+prg_options['mod_version']+'/'+'-'.join(prg_options['g_name'])+'/'+(prg_options['motif_of_interest']+"/" if 'motif_of_interest' in prg_options else "")
   _len_str_ = str(int(sys.argv[3]))+"/"; #(str(int(sys.argv[3]))+"/") if ((len(sys.argv)>3) and (not int(sys.argv[3])==3)) else ""
   prg_options['save_mod_path'] = 'repeat_mod_output_v'+prg_options['mod_version']+'/'+'-'.join(prg_options['g_name'])+'/'+(prg_options['motif_of_interest']+"_"+_len_str_+"/" if 'motif_of_interest' in prg_options else _len_str_)
   if (not os.path.isdir(prg_options['save_mod_path'])):
      os.system('mkdir -p '+prg_options['save_mod_path'])
 
   p_options = {}
   p_options['repeat_unit_size'] = 3;
   if len(sys.argv)>3:
      p_options['repeat_unit_size'] = int(sys.argv[3]);

   prg_options['UniqueID'] = 'np_dl_'+str(p_options['repeat_unit_size'])

   if len(sys.argv)>6:
      p_options['lr'] = int(sys.argv[6])
   else: p_options['lr'] = 5;
   if len(sys.argv)>7:
      p_options["max_epoch"] = int(sys.argv[7])
   if len(sys.argv)>8:
      p_options["batch_size"] = int(sys.argv[8])

   p_options['c_fs_file_list'] = [];
   for c_fs_folder in prg_options['g_name']:
      _this_fs_file_list = glob.glob(os.path.join('nanopore_data/fs/'+c_fs_folder+'/motif/', ('*_'+str(p_options['repeat_unit_size'])+'.fs.gz' if 'motif_of_interest' not in prg_options else prg_options['motif_of_interest']+'*_'+str(p_options['repeat_unit_size'])+'.fs.gz')))
      for _m_t_ff in _this_fs_file_list:
         if (_m_t_ff.split('/')[-1] not in ['ACG_3.fs.gz', 'ACG_3.fs']) or ('motif_of_interest' in prg_options and prg_options['motif_of_interest'] in ['ACG']):
            print(_m_t_ff)
            p_options['c_fs_file_list'].append(_m_t_ff)
      #p_options['c_fs_file_list'].extend(glob.glob(os.path.join('nanopore_data/fs/'+c_fs_folder+'/motif/', ('*.fs.gz' if 'motif_of_interest' not in prg_options else prg_options['motif_of_interest']+'*.fs.gz'))));
   p_options['c_fs_file_list'] = sorted(p_options['c_fs_file_list'])

   m_train(p_options, prg_options)

