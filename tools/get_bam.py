
import os, sys
import glob

if len(sys.argv)<5:
   print("Usage: "+sys.argv[0]+" fq_folder bam_folder suf_str log_subf")
   print("\tpython "+sys.argv[0]+" strique_g4c2_data/m_fq/ bam_results/ \".fq\" g4c2 ref"); 
   sys.exit(1);

fq_folder = sys.argv[1]
bam_folder = sys.argv[2];

suf_str = sys.argv[3]
if not suf_str[0]=='.': suf_str = '.'+suf_str;
fq_files = glob.glob(os.path.join(fq_folder, "*"+suf_str));
#suf_str = '.minimap2.sorted.fq'
#suf_str = sys.argv[3]

ref_file = 'hg_ref/hg38/hg38.fa'
if len(sys.argv)>5:
   ref_file = sys.argv[5]

for fqf in fq_files:
   print (fqf);
   pref_str = fqf.split('/')[-1][:-len(suf_str)]
   
   #cmd = 'echo "minimap2 -ax map-ont %s %s | samtools sort | samtools view -bS > %s/%s.bam ; samtools index %s/%s.bam " | qsub -V -cwd -N mp2_%s -o log/%s/mp2_%s.log -e log/%s/mp2_%s.log -l h_vmem=10G' % (ref_file,fqf, bam_folder,pref_str, bam_folder,pref_str, pref_str,pref_str,sys.argv[4],pref_str,sys.argv[4])
   cmd = "minimap2 -ax map-ont %s %s | samtools sort | samtools view -bS > %s/%s.bam ; samtools index %s/%s.bam ".format(ref_file,fqf, bam_folder,pref_str, bam_folder,pref_st)
   print(cmd)
   sys.stdout.flush()
   os.system(cmd);

