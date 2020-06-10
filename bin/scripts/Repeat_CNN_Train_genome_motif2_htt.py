
import glob
import os, sys;

import Repeat_CNN_Train_genome_motif2

if __name__=='__main__':
   prg_options = {}
   prg_options['gpu_cpu'] = 'gpu'
   #prg_options['gpu_cpu'] = 'cpu'
   prg_options['gpu_cpu'] = sys.argv[5]
   prg_options['feature_num'] = 50;
   prg_options['nb_size'] = 3;
   prg_options['label_size'] = 4;
   prg_options['num_filter_in_first'] = int(sys.argv[2])
   prg_options['motif_of_interest'] = sys.argv[4]
   prg_options['mod_version'] = '0.2motif'+str(prg_options['num_filter_in_first']); # flexible filters for last cnn layers
   prg_options['g_name'] = ['hx1']
   prg_options['g_name'] = sys.argv[1].split(',')
   _len_str_ = str(int(sys.argv[3]))+"/"
   prg_options['save_mod_path'] = 'sep_mod_output_v'+prg_options['mod_version']+'/'+'-'.join(prg_options['g_name'])+'/'+prg_options['motif_of_interest']+"_"+sys.argv[6]+"_"+str(int(sys.argv[7]))+"/"
   if (not os.path.isdir(prg_options['save_mod_path'])):
      os.system('mkdir -p '+prg_options['save_mod_path'])

   p_options = {}
   p_options['repeat_unit_size'] = 3;
   p_options['repeat_unit_size'] = int(sys.argv[3]);

   #p_options["max_epoch"] = int(sys.argv[8])
   if len(sys.argv)>8:
      p_options['lr'] = int(sys.argv[8])
   else: p_options['lr'] = 5;
   if len(sys.argv)>9:
      p_options["max_epoch"] = int(sys.argv[9])
   if len(sys.argv)>10:
      p_options["batch_size"] = int(sys.argv[10])

   prg_options['UniqueID'] = 'np_dl_'+str(p_options['repeat_unit_size'])

   fs_file_list = glob.glob(os.path.join(sys.argv[6]+'_data/fs_results/', '*.fs'));
   fs_file_list = sorted(fs_file_list)
   if int(sys.argv[7]) >= len(fs_file_list):
      p_options['c_fs_file_list'] = fs_file_list
   else:
      p_options['c_fs_file_list'] = fs_file_list[:int(sys.argv[7])]
      print ("Training: {}".format(fs_file_list[:int(sys.argv[7])]))
      print ("Remained: {}".format(fs_file_list[int(sys.argv[7]):]))

   Repeat_CNN_Train_genome_motif2.m_train(p_options, prg_options)


