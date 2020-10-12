
import os,sys

barcode_folder = sys.argv[1]
fq_save_folder = sys.argv[2];

_bc_folders = os.listdir(barcode_folder);

_bc_str = 'barcode'

for _bcf in _bc_folders:
   #print (_bcf)
   if os.path.isdir(barcode_folder+"/"+_bcf) and _bcf[:len(_bc_str)]==_bc_str:
      _s_fq_f = fq_save_folder+"/"+_bcf+".fq"
      _cmd = "cat "+barcode_folder+"/"+_bcf+"/*.fastq > "+_s_fq_f
      print(_cmd)
      os.system(_cmd)

