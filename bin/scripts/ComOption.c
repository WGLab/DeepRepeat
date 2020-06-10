#include "ComOption.h"

//#include "string.h"
//#include <unistd.h>
#include <string>
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <sstream>
#include <fstream>
#include <iostream>

#include "ComFunction.h"


//RunComOption m_runComOption;
//RunOption m_runOption;
ComOption m_comOption;


std::map<std::string, bool(*)()> check_map;
std::map<std::string, bool(*)(char*)> set_map;
std::map<std::string, void(*)()> print_map;
std::map<std::string, void(*)()> print_help_map;
std::map<std::string, std::string> op_labels;
std::map<std::string, std::string> op_labels2;

int ComOption::max_length_of_repeat_unit = 10;
char* ComOption::module_type[MODULE_NUM] = {(char*)"f5index", (char*)"getfeat", (char*)"train", (char*)"predict", (char*)"detect"};

const char* op_label_getfeat[] = {"-l",
                            "-m",
                            "-f",
                            "-q",
                            "-b",
                            "-p",
                            "-i",
                            "-s",
                            "-r",
                            "-u"};
const char* op_label_train[] = {"-g",
                          "-s",
                          "-u",
                          "-e"};
const char* op_label_predict[] = {"-g",
                            "-s",
                            "-u",
                            "-n"};
const char* op_label_detect[] = {"-l",
                           "-m",
                           "-f",
                           "-q",
                           "-b",
                           "-p",
                           "-n",
                           "-i",
                           "-s",
                           "-r",
                           "-u"};
const char* op_label_f5index[] = {"-t",
                            "-m",
                            "-u",
                            "-c",
                            "-j"};
//int op_label_elem_len[MODULE_NUM] = {3, 10, 4, 4, 11};
int op_label_elem_len[MODULE_NUM] = {get_array_size(op_label_f5index),
                                     get_array_size(op_label_getfeat),
                                     get_array_size(op_label_train),
                                     get_array_size(op_label_predict),
                                     get_array_size(op_label_detect) };

const char** op_label5[] = {op_label_f5index, op_label_getfeat, op_label_train, op_label_predict, op_label_detect};
const char* default_uniqid[] = {"f5index", "getfeat", "rep_mod", "rep_pred", "rep_det"};


bool m_create_folder(char * mpath){
   if (isExist(mpath)){
      return true;
   }

   char * mk_pref = (char*)"mkdir -p ";
   char * mcmd = (char*)malloc(1+strlen(mk_pref)+strlen(mpath));
   strcpy(mcmd, mk_pref);
   strcpy(mcmd, mpath);
   const int mkerr = system(mcmd);
   free(mcmd);
   if (mkerr==-1){
      fprintf(stderr, "\tError: Cannot create %s\n", mpath);
      return false;
   }
   return true;
}

void init_ComOption()
{
   int cp_res;
   m_comOption.rep_info.len_repeat_unit = 3;
   m_comOption.rep_info.chrn[0] = '\0';
   m_comOption.rep_info.start_pos = -1;
   m_comOption.rep_info.end_pos = -1 ;
   m_comOption.multifast5 = false ;
   m_comOption.fast5_folder[0] = '\0';
   m_comOption.fastq_folder[0] = '\0';
   m_comOption.bam_file[0] = '\0';
   m_comOption.save_path[0] = '\0';
   snprintf(m_comOption.uniq_id, CHAR_SIZE, "m_rep");
   //m_comOption.uniq_id = (char*)"m_rep";
   m_comOption.submod_id = -1;
   m_comOption.submod_str[0] = '\0';
   m_comOption.feat_path[0] = '\0';
   m_comOption.mod_path[0] = '\0';
   m_comOption.epoch = 10;
   m_comOption.seq_sum[0] = '\0';
   m_comOption.fast5_index_file[0] = '\0';
   //m_comOption.basecalled_path = (char*)"workspace/pass/";
   snprintf(m_comOption.basecalled_path, CHAR_SIZE, "workspace/pass/");
   m_comOption.base_index_path[0] = '\0';
   
   m_comOption.group_size = 4;
}

void print_mod()
{
   fprintf(stdout, "\t f5index: Build index for fast5 files. \n");
   fprintf(stdout, "\t getfeat: Get features from fast5 files. \n");
   fprintf(stdout, "\t train: Train a deep learning model with extracted features. \n");
   fprintf(stdout, "\t predict: Predict genomic positins in repeats with extracted features. \n");
   fprintf(stdout, "\t detect: Detect repeat counts from fast5 files. \n");
}
int get_submod_id(char *sub_mod)
{
   int which_mod = -1;
   int mdi = 0;
   for (; mdi<MODULE_NUM; mdi++){
      if (strcmp(sub_mod,ComOption::module_type[mdi])==0){
         which_mod = mdi;
         break;
      }
   }
   return which_mod;
}

const char** get_op_label(int submod_id)
{
   if (submod_id<0 or submod_id>=MODULE_NUM){
      return NULL;
   }
   return op_label5[submod_id];
   //fprintf(stderr, "module_id=%d must be between 0 and %d: 0-%s 1-%s 2-%s 3-%s", submod_id, MODULE_NUM-1, ComOption::module_type[0], ComOption::module_type[1], ComOption::module_type[2], ComOption::module_type[3]);
}


void set_submod_id(char * sub_mod)
{
   m_comOption.submod_id = get_submod_id(sub_mod);
   //m_comOption.submod_str = sub_mod;
   if (!m_cp_str(m_comOption.submod_str, sub_mod, CHAR_SIZE)){return;}
}

void print_repdict(){
   std::map<std::string, std::map<std::string, RepeatRegion> >::iterator chr_dit;
   std::map<std::string, RepeatRegion>::iterator pos_it;
   for (chr_dit=m_comOption.repdict.begin(); chr_dit!=m_comOption.repdict.end(); chr_dit++){
      for (pos_it=chr_dit->second.begin(); pos_it!=chr_dit->second.end(); pos_it++){
         std::cout << "\t\t" << pos_it->second.chrn << ":"<<pos_it->second.start_pos<<"-"<<pos_it->second.end_pos<<"="<<pos_it->second.len_repeat_unit;
      }
   }
}
bool check_repdict(){
   return check_repeat_pos();
}
bool set_repdict(char* mrepfile){
   if (mrepfile!=NULL){
      std::map<std::string, std::map<std::string, RepeatRegion> >::iterator chr_dit;
      std::map<std::string, RepeatRegion>::iterator pos_it;
      std::string line;
      std::ostringstream oss_chr;
      std::ostringstream oss_pos;

      std::ifstream infile(mrepfile);
      while (std::getline(infile, line)){
          RepeatRegion crr;
          oss_chr.clear();
          oss_pos.clear();

          std::istringstream iss(line);
          if (!(iss >> crr.chrn >> crr.start_pos >> crr.end_pos >> crr.len_repeat_unit)) { break; }
          if (crr.end_pos-crr.start_pos<1 || crr.len_repeat_unit<1){
              std::cout<<"\t Warning!!! Incorrect repeat info="<<line<<std::endl;
              continue;
          }
          oss_chr<<crr.chrn;
          oss_pos<<crr.start_pos<<"-"<<crr.end_pos;
          chr_dit = m_comOption.repdict.find(oss_chr.str());
          if (chr_dit==m_comOption.repdict.end()){
             std::map<std::string, RepeatRegion> new_rr_dict;
             new_rr_dict[oss_pos.str()] = crr;
             m_comOption.repdict[oss_chr.str()] = new_rr_dict;
          }else{
             pos_it = chr_dit->second.find(oss_pos.str());
             if (pos_it==chr_dit->second.end()){
                 (chr_dit->second)[oss_pos.str()] = crr;
             }else{
                 std::cout<<"\t Warning!!! Repeat info already in the dict: "<<oss_chr.str()<<":"<<oss_pos.str()<< " ??? "<<chr_dit->first << ":" << pos_it->first<<std::endl;
             }
          }
      }
      infile.close();
   }
   return check_repdict();
}
void print_help_repdict(){
   fprintf(stdout, "\t%s %-14s\t%s\n", "-r", "--repfile", "A file for a list of repeat positions. Default: NULL");
}

void print_seq_sum(){
   if (m_comOption.seq_sum[0]!='\0'){
      std::cout << " Seq_sum=" << m_comOption.seq_sum;
      //fprintf(stdout, " Seq_sum=%", m_comOption.seq_sum);
   }else{
      std::cout << " Seq_num=NULL";
      //fprintf(stdout, " Seq_num=NULL");
   }
}
bool check_seq_sum(){
   if (m_comOption.seq_sum[0]=='\0'){
      fprintf(stderr, "\tError: Fast5 summary file is not provided!\n");
      return false;
   }
   return true;
}
bool set_seq_sum(char* mseq_sum){
   if (!m_cp_str(m_comOption.seq_sum, mseq_sum, CHAR_SIZE)){return false;}
   return check_seq_sum();
}
void print_help_seq_sum(){
   fprintf(stdout, "\t%s %-14s\t%s\n", "-t", "--seq_sum", "Path for module file. Could be \"data/*/*/sequencing_summary.txt\"(\" is necessary if '*' is in the paramater). Fast5 files for being indexed must have the same parent directory as corresponding 'sequencing_summary.txt'. Default: NULL");
}

void print_basecalled_path(){
   std::cout<< " Basecalled_path=" << m_comOption.basecalled_path;
   //fprintf(stdout, " Basecalled_path=%", m_comOption.basecalled_path);
}
bool check_basecalled_path(){
   if (m_comOption.basecalled_path[0]=='\0'){
      fprintf(stderr, "\tError: Path for fast5 files from sequencing_summary.txt is not provided!\n");
      return false;
   }
   return true;
}
bool set_basecalled_path(char* mbasecalled_path){
   if (!m_cp_str(m_comOption.basecalled_path, mbasecalled_path, CHAR_SIZE)){return false;}
   //m_comOption.basecalled_path = mbasecalled_path;
   return check_basecalled_path();
}
void print_help_basecalled_path(){
   fprintf(stdout, "\t%s %-14s\t%s\n", "-c", "--basecallPath", "Path for fast5 files from sequencing_summary.txt. Default: 'workspace/pass/'");
}

void print_fast5_index_file(){
   if (m_comOption.fast5_index_file[0]!='\0'){
      std::cout<< " Fast5 index file=" << m_comOption.fast5_index_file;
      //fprintf(stdout, " Fast5 index file=%", m_comOption.fast5_index_file);
   }else{
      std::cout<< " Fast5 index file=NULL";
      //fprintf(stdout, " Fast5 index file=NULL");
   }
}
bool check_fast5_index_file(){
      if (m_comOption.fast5_index_file[0]=='\0'){
      fprintf(stderr, "\tError: Fast5 index file is not provided!\n");
      return false;
   }else if (!isExist(m_comOption.fast5_index_file)){
      fprintf(stderr, "\tError: Fast5 index file (%s) does not exist!\n", m_comOption.fast5_index_file);
      return false;
   }
   return true;
}
bool set_fast5_index_file(char* mfast5_index_file){
   if (!m_cp_str(m_comOption.fast5_index_file, mfast5_index_file, CHAR_SIZE)){return false;}
   return check_fast5_index_file();
}
void print_help_fast5_index_file(){
   fprintf(stdout, "\t%s %-14s\t%s\n", "-i", "--fast5index", "Index file for fast5 files. Default: NULL");
}

void print_mod_path()
{
   if (m_comOption.mod_path[0]!='\0'){
      std::cout<< " Mod_path=" << m_comOption.mod_path;
      //fprintf(stdout, " Mod_path=%", m_comOption.mod_path);
   }else{
      std::cout<< " Mod_path=NULL";
      //fprintf(stdout, " Mod_path=NULL");
   }
}
bool check_mod_path()
{
   if (m_comOption.mod_path[0]=='\0'){
       fprintf(stderr, "\tError: Mod path is not provided!\n");
      return false;
   }else if (!isExist(m_comOption.mod_path)){
      fprintf(stderr, "\tError: Mod path (%s) does not exist!\n", m_comOption.mod_path);
      return false;
   }
   return true;
}
bool set_mod_path(char* mmod_path)
{
   if (!m_cp_str(m_comOption.mod_path, mmod_path, CHAR_SIZE)){return false;}
   //m_comOption.mod_path = mmod_path;
   return check_mod_path();
}
void print_help_mod_path()
{
   fprintf(stdout, "\t%s %-14s\t%s\n", "-n", "--mod_path", "Path for module file. Default: NULL");
}

void print_epoch()
{
   std::cout<<" Epoch=" << m_comOption.epoch;
   //fprintf(stdout, " Epoch=%d", m_comOption.epoch);
}
bool check_epoch()
{
   if (m_comOption.epoch < 1) {
      fprintf(stderr, "\tError: Epoch=%d is less than 1!\n", m_comOption.epoch);
      return false;
   }
   return true;
}
bool set_epoch(char* mepoch)
{
   m_comOption.epoch = atoi(mepoch);
   return check_epoch();
}
void print_help_epoch()
{
   fprintf(stdout, "\t%s %-14s\t%s\n", "-e", "--epoch", "Epoch for data used for training. Default: 10");
}

void print_feat_path()
{
   if (m_comOption.feat_path[0]!='\0'){
      std::cout<<" Feat_path="<< m_comOption.feat_path;
      //fprintf(stdout, " Feat_path=%", m_comOption.feat_path);
   }else{
      std::cout<<" Feat_path=NULL";
      //fprintf(stdout, " Feat_path=NULL");
   }
}
bool check_feat_path()
{
   if (m_comOption.feat_path[0]=='\0'){
      fprintf(stderr, "\tError: Feature-path is not provided!\n");
      return false;
   }else if (!isExist(m_comOption.feat_path)){
      fprintf(stderr, "\tError: Feature-path (%s) does not exist!\n", m_comOption.feat_path);
      return false;
   }
   return true;
}
bool set_feat_path(char* mfeat_path)
{
   if (!m_cp_str(m_comOption.feat_path, mfeat_path, CHAR_SIZE)){return false;}
   //m_comOption.feat_path = mfeat_path;
   return check_feat_path();
}
void print_help_feat_path()
{
   fprintf(stdout, "\t%s %-14s\t%s\n", "-g", "--feat_path", "Path for feature file. Default: NULL");
}

void print_len_repeat_unit()
{
   std::cout<< " Repeat-length=" << m_comOption.rep_info.len_repeat_unit;
   //fprintf(stdout, " Repeat-length=%d", m_comOption.rep_info.len_repeat_unit);
}
bool check_len_repeat_unit()
{
   if (m_comOption.rep_info.len_repeat_unit<1){
      fprintf(stderr, "\tError: Length=%d of Repeat Unit cannot be less than 1.\n", m_comOption.rep_info.len_repeat_unit);
      return false;
   }else if (m_comOption.rep_info.len_repeat_unit>ComOption::max_length_of_repeat_unit){
      fprintf(stderr, "\tError: Length=%d(>%d) of Repeat Unitis is not supported now. Will be extended in future.\n", m_comOption.rep_info.len_repeat_unit, ComOption::max_length_of_repeat_unit);
      return false;
   }else if (m_comOption.rep_info.len_repeat_unit<3){
      fprintf(stdout, "\tWarning!!! Length=%d of Repeat Unitis is too small and might cause more errors.\n", m_comOption.rep_info.len_repeat_unit);
      return true;
   }
   return true;
}
bool set_len_repeat_unit(char* mlen)
{
   m_comOption.rep_info.len_repeat_unit = atoi(mlen);
   return check_len_repeat_unit();
}
void print_help_len_repeat_unit()
{
   fprintf(stdout, "\t%s %-14s\tLenth of repeat unit ( from 1 to %d ). Default: 3.\n", "-l", "--replen", ComOption::max_length_of_repeat_unit);
}

void print_multifast5()
{
   std::cout<< " multifast5=" << m_comOption.multifast5?"True":"False";
   //fprintf(stdout, " multifast5=%s", m_comOption.multifast5?"True":"False");
}
bool check_multifast5()
{
   return true;
}
bool set_multifast5(char* mf5)
{
   m_comOption.multifast5 = atoi(mf5);
   return check_multifast5();
}
void print_help_multifast5()
{
   fprintf(stdout, "\t%s %-14s\t%s\n", "-m", "--multifast5", "Is multifast5 format. Default: false.");
}

void print_fast5_folder()
{
   if (m_comOption.fast5_folder[0]!='\0'){
      std::cout<< " Fast5 folder=" << m_comOption.fast5_folder;
      //fprintf(stdout, " Fast5 folder=%s", m_comOption.fast5_folder);
   }else{
      std::cout<< " Fast5 folder=NULL";
      //fprintf(stdout, " Fast5 folder=NULL");
   }
}
bool check_fast5_folder()
{
   if (m_comOption.fast5_folder[0]=='\0'){
      fprintf(stderr, "\tError: Fast5 Folder is not provided.\n");
      return false;
   }else if (!isExist(m_comOption.fast5_folder)){
      fprintf(stderr, "\tError: Fast5 Folder (%s) does not exist.\n", m_comOption.fast5_folder);
      return false;
   }
   return true;
}
bool set_fast5_folder(char* f5_folder)
{
   if (!m_cp_str(m_comOption.fast5_folder, f5_folder, CHAR_SIZE)){return false;}
   //m_comOption.fast5_folder = f5_folder;
   return check_fast5_folder();
}
void print_help_fast5_folder()
{
   fprintf(stdout, "\t%s %-14s\t%s\n", "-f", "--fast5", "Fast5 folder (input). Default: NULL.");;
}

void print_fastq_folder()
{
   if (m_comOption.fastq_folder[0]!='\0'){
      std::cout<< " Fastq folder=" << m_comOption.fastq_folder;
      //fprintf(stdout, " Fastq folder=%s", m_comOption.fastq_folder);
   }else{
      std::cout<< " Fastq folder=NULL";
      //fprintf(stdout, " Fastq folder=NULL");
   }
}
bool check_fastq_folder()
{
   if (m_comOption.fastq_folder[0]=='\0'){
       if (m_comOption.bam_file[0]=='\0'){
          fprintf(stderr, "\tError: Both fastq folder and bam file are not provided.\n");
          return false;
       }
   }else if (!isExist(m_comOption.fastq_folder)){
      fprintf(stderr, "\tError: Fastq Folder (%s) does not exist.\n", m_comOption.fastq_folder);
      return false;
   }
   return true;
}
bool set_fastq_folder(char* fq_folder)
{
   if (!m_cp_str(m_comOption.fastq_folder, fq_folder, CHAR_SIZE)){return false;}
   //m_comOption.fastq_folder = fq_folder;
   return check_fastq_folder();
}
void print_help_fastq_folder()
{
   fprintf(stdout, "\t%s %-14s\t%s\n", "-q", "--fastq", "Fastq folder (input). Default: NULL.");
}

void print_bam_file()
{
   if (m_comOption.bam_file[0]!='\0'){
      std::cout<< " Bam file=" << m_comOption.bam_file;
      //fprintf(stdout, " Bam file=%s", m_comOption.bam_file);
   }else{
      std::cout<< " Bam file=NULL";
      //fprintf(stdout, " Bam file=NULL");
   }
}
bool check_bam_file()
{
   if (m_comOption.bam_file[0]=='\0'){
       if (m_comOption.fastq_folder[0]=='\0'){
          fprintf(stderr, "\tError: Both fastq folder and bam file are not provided.\n");
          return false; 
       }
   }else if (!isExist(m_comOption.bam_file)){
      fprintf(stderr, "\tError: Bam file (%s) does not exist.\n", m_comOption.bam_file);
      return false; 
   }
   return true;
}
bool set_bam_file(char* bam_file)
{
   if (!m_cp_str(m_comOption.bam_file, bam_file, CHAR_SIZE)){return false;}
   //m_comOption.bam_file = bam_file;
   return check_bam_file();
}
void print_help_bam_file()
{
   fprintf(stdout, "\t%s %-14s\t%s\n", "-b", "--bam", "Bam file (input). Default: NULL.");
}

void print_repeat_pos()
{
   std::cout<< " Genomice position="<<m_comOption.rep_info.chrn<<":"<<m_comOption.rep_info.start_pos<<"-"<<m_comOption.rep_info.end_pos;
   ///fprintf(stdout, " Genomice position=%s:%d-%d", m_comOption.rep_info.chrn,m_comOption.rep_info.start_pos,m_comOption.rep_info.end_pos);
}
bool check_repeat_pos()
{
   if (m_comOption.rep_info.chrn[0]=='\0' && m_comOption.repdict.empty()){
      fprintf(stderr, "\tError: Chromosome is not provided.\n");
      return false;
   }else if(m_comOption.rep_info.start_pos<0 && m_comOption.repdict.empty()){
      fprintf(stderr, "\tError: No start position is provided.\n");
      return false;
   }else if(m_comOption.rep_info.end_pos<0 && m_comOption.repdict.empty()){
      fprintf(stderr, "\tError: No end position is provided.\n");
      return false;
   }else if (m_comOption.rep_info.start_pos>=m_comOption.rep_info.end_pos && m_comOption.repdict.empty()){
      fprintf(stderr, "\tError: Start position=%d must be smaller than end position=%d.\n", m_comOption.rep_info.start_pos,m_comOption.rep_info.end_pos);
      return false;
   }
   
   std::map<std::string, std::map<std::string, RepeatRegion> >::iterator chr_dit;
   std::map<std::string, RepeatRegion>::iterator pos_it;
   std::ostringstream oss_chr;
   std::ostringstream oss_pos;
   oss_chr<<m_comOption.rep_info.chrn;
   oss_pos<<":"<<m_comOption.rep_info.start_pos<<"-"<<m_comOption.rep_info.end_pos;
   chr_dit = m_comOption.repdict.find(oss_chr.str());
   if (chr_dit==m_comOption.repdict.end()){
      std::map<std::string, RepeatRegion> new_rr_dict;
      new_rr_dict[oss_pos.str()] = m_comOption.rep_info;
      m_comOption.repdict[oss_chr.str()] = new_rr_dict;
   }else{
      pos_it = chr_dit->second.find(oss_pos.str());
      if (pos_it==chr_dit->second.end()){
          (chr_dit->second)[oss_pos.str()] = m_comOption.rep_info;
      }
   }

   return true;
}
bool get_region_from_str(char* pos_str, RepeatRegion * mrr){
   if (pos_str==NULL){
       fprintf(stderr, "\tError: No genomic position is provided.\n");
       return false;
   }
   char *pch;
   pch = strtok (pos_str,":-");
   //mrr->chrn = pch;
   if (!m_cp_str(mrr->chrn, pch, CHAR_SIZE)){return false;}
   pch = strtok (NULL, ":-");
   if (pch==NULL){
      fprintf(stderr, "\tError: No start position is provided.\n");
      return false;
   }
   mrr->start_pos = atoll(pch);

   pch = strtok (NULL, ":-");
   if (pch==NULL){
      fprintf(stderr, "\tError: No end position is provided.\n");
      return false;
   }
   mrr->end_pos = atoll(pch);
   
   return true;
}
bool set_repeat_pos(char* pos_str)
{
   get_region_from_str(pos_str, &(m_comOption.rep_info)); 
   return check_repeat_pos();
}
void print_help_repeat_pos()
{
   fprintf(stdout, "\t%s %-14s\t%s\n", "-p", "--pos", "Genomic region of repeats. Required. Format=chr:start-end. Default: NULL.");
}

void print_save_path()
{
   std::cout<< " Save path=" << m_comOption.save_path;
   //fprintf(stdout, " Save path=%s", m_comOption.save_path);
}
bool check_save_path()
{
   if (m_comOption.save_path[0]=='\0'){
      fprintf(stderr, "\tError: No save path is provided.\n");
      return false;
   }
   return m_create_folder(m_comOption.save_path);
}
bool set_save_path(char* svpath)
{
   if (!m_cp_str(m_comOption.save_path, svpath, CHAR_SIZE)){return false;}
   //m_comOption.save_path = svpath;
   return check_save_path();
}
void print_help_save_path()
{
   fprintf(stdout, "\t%s %-14s\t%s\n", "-s", "--savepath", "Output directory path(output). Required. Default: NULL.");
}

void print_uniq_id()
{
   std::cout<<" Unique file id="<<m_comOption.uniq_id;
   //fprintf(stdout, " Unique file id=%s", m_comOption.uniq_id);
}
bool check_uniq_id( )
{
   if(m_comOption.uniq_id[0]=='\0'){
      fprintf(stderr, "\tError: Unique name is not provided.\n");
      return false;
   }
   return true;
}
bool set_uniq_id(char* muid)
{
   if (!m_cp_str(m_comOption.uniq_id, muid, CHAR_SIZE)){return false;}
   //m_comOption.uniq_id = muid;
   return check_uniq_id();
}
void print_help_uniq_id()
{
   if (m_comOption.submod_id<0 || m_comOption.submod_id>=MODULE_NUM){
      fprintf(stdout, "\t%s %-14s\t%s\n", "-u", "--uniq_id", "Unique file id for output. Default: 'f5index_' for 'f5index', 'getfeat_' for 'getfeat', 'rep_mod_' for 'train', 'rep_pred_' for 'predict',  'rep_det_' for 'detect'.");
   }else{
      fprintf(stdout, "\t%s %-14s\tUnique file id for output. Default: '%s' for <%s>\n", "-u", "--uniq_id", default_uniqid[m_comOption.submod_id], ComOption::module_type[m_comOption.submod_id]);
   }
}

void print_base_index_path(){
   if (m_comOption.base_index_path[0]!='\0'){
      std::cout<<" Base index path="<<m_comOption.base_index_path;
      //fprintf(stdout, " Base index path=%", m_comOption.base_index_path);
   }else{
      std::cout<<" Base index path=NULL";
      //fprintf(stdout, " Base index path=NULL");
   }
}
bool check_base_index_path(){
   if (m_comOption.base_index_path[0]!='\0' && !isExist(m_comOption.base_index_path)){
      std::cout<< "Error!!! The base index path ("<< m_comOption.base_index_path << ") is not provided or does not exist.\n" << std::endl;
      return false;
   }
   return true;
}
bool set_base_index_path(char* mbase_index_path){
   if (!m_cp_str(m_comOption.base_index_path, mbase_index_path, CHAR_SIZE)){return false;}
   return check_base_index_path();
}
void print_help_base_index_path(){
  fprintf(stdout, "\t%s %-14s\t%s\n", "-j", "--baseIndexPath", "Base path for index fast5 only. Must be the prefix of '-t'. Default: NULL");
}

void print_help(char* m_prog, char * sub_mod)
{
   if (sub_mod==NULL){
      fprintf(stdout, "Usage: %s [f5index|getfeat|train|predict|detect]\n", m_prog);
      print_mod();
      return;
   }

   int which_mod = get_submod_id(sub_mod);
   if (which_mod<0){
      fprintf(stderr, "Error !!! <%s> is not supported.\n", sub_mod);
      fprintf(stdout, "Usage: %s [f5index|getfeat|train|predict|detect]\n", m_prog);
      print_mod();
      return;
   }

   m_comOption.submod_id = which_mod;
   if (!m_cp_str(m_comOption.submod_str, sub_mod, CHAR_SIZE)){return;}
   fprintf(stdout, "Usage: %s %s\n", m_prog, ComOption::module_type[which_mod]);
   const char** cur_op_label = get_op_label(which_mod);
   int opi = 0;
   for (; opi<op_label_elem_len[which_mod]; opi++){
      (*(print_help_map[cur_op_label[opi]]))();
   }
}

void print_options()
{
   if (m_comOption.submod_id<0){
      fprintf(stderr, "Error!!! Module_id=%d must be between 0 and %d: 0-%s 1-%s 2-%s 3-%s 4-%s\n", m_comOption.submod_id, MODULE_NUM-1, ComOption::module_type[0], ComOption::module_type[1], ComOption::module_type[2], ComOption::module_type[3], ComOption::module_type[4]);
      return;
   }

   const char** cur_op_label = get_op_label(m_comOption.submod_id);
   int opi = 0;
   for (; opi<op_label_elem_len[m_comOption.submod_id]; opi++){
      (print_map[cur_op_label[opi]])();
      fprintf(stdout, "\n");
   }
}


bool set_options(int argc, char * argv[])
{
   bool all_correct = true;

   set_submod_id(argv[1]);
   if (m_comOption.submod_id<0){
      return true;
   }

   const char** cur_op_label = get_op_label(m_comOption.submod_id);

   int argi = 2;
   int opi;
   int set_bit = 0;
   std::map<std::string, std::string>::iterator op_it;
   std::string cur_opl;
   while (argi<argc){
      if (argv[argi][0]=='-'){
         if (argi+1>=argc || (argv[argi+1][0]=='-')){
            fprintf(stderr, "\tNo value is provided for <%s>.\n", argv[argi]);
            all_correct = false;
            argi++;
            continue;
         }

         cur_opl.clear();
         if (strlen(argv[argi])<2){
            fprintf(stderr, "\t<%s> is not supported.\n", argv[argi]);
            all_correct = false;
            if (argi+1<argc && argv[argi+1][0]!='-'){
               argi++;
            }
         }else{
            if (argv[argi][1]=='-'){
               op_it = op_labels2.find(argv[argi]);
               if (op_it==op_labels2.end()){
                  fprintf(stderr, "\t<%s> is not supported.\n", argv[argi]);
                  all_correct = false;
                  if (argi<argc && argv[argi+1][0]!='-'){
                     argi++;
                  }
               }else{
                  cur_opl = op_it->second;
               }
            }else{
               cur_opl = argv[argi];
            }
            if (cur_opl.size()>0){
               argi++;
               opi = 0;
               for (; opi<op_label_elem_len[m_comOption.submod_id]; opi++){
                  if (cur_opl.compare(cur_op_label[opi])==0){
                     //std::cout << "Now:" << opi << " " << cur_opl << " " << cur_op_label[opi] << " " << (*(set_map[cur_op_label[opi]])) << " " << all_correct << std::endl;
                     all_correct = ((set_map[cur_op_label[opi]])(argv[argi])) && all_correct;
                     //(print_map[cur_op_label[opi]])(); std::cout<<std::endl;
                     //std::cout<< " for check " << m_comOption.seq_sum << " " << m_comOption.basecalled_path << " " << m_comOption.base_index_path << " -----" << &(m_comOption.seq_sum) << " " << &(m_comOption)  << std::endl;
                     //print_seq_sum(); std::cout<<std::endl;
                     set_bit |= 1<<opi;
                     break;
                  }
               }
               if (opi == op_label_elem_len[m_comOption.submod_id]){
                  fprintf(stderr, "\t<%s> is not supported.\n", cur_opl.c_str());
                  all_correct = false;
               }
            }
         }
      }else{
         fprintf(stderr, "\t<%s> is not supported.\n", argv[argi]);
         all_correct = false;
      }
      argi++;
   }

   if (set_bit==0){
      return true;
   }
   opi = 0;
   for (; opi<op_label_elem_len[m_comOption.submod_id]; opi++){
      if ((set_bit & (1<<opi))==0){
         //std::cout<< " not set " << cur_op_label[opi] << " check now" <<std::endl;
         all_correct = ((check_map[cur_op_label[opi]])()) && all_correct;
      }
   }

   if (strcmp("m_rep", m_comOption.uniq_id)==0){
      if (!(m_comOption.submod_id<0 || m_comOption.submod_id>=MODULE_NUM)){
          snprintf(m_comOption.uniq_id, CHAR_SIZE, default_uniqid[m_comOption.submod_id]);
      }
   }

   std::cout<< "The setting below are provided or set by default. " << std::endl;
   print_options();

   return !all_correct;
}

void set_function()
{
   set_map["--replen"] = set_len_repeat_unit;
   set_map["-l"] = set_len_repeat_unit;                        op_labels["-l"] = "--replen";
   print_map["-l"] = print_len_repeat_unit;
   print_help_map["-l"] = print_help_len_repeat_unit;
   check_map["-l"] = check_len_repeat_unit;

   set_map["--multifast5"] = set_multifast5;
   set_map["-m"] = set_multifast5;                             op_labels["-m"] = "--multifast5";
   print_map["-m"] = print_multifast5;
   print_help_map["-m"] = print_help_multifast5;
   check_map["-m"] = check_multifast5;

   set_map["--fast5"] = set_fast5_folder;
   set_map["-f"] = set_fast5_folder;                           op_labels["-f"] = "--fast5";
   print_map["-f"] = print_fast5_folder;
   print_help_map["-f"] = print_help_fast5_folder;
   check_map["-f"] = check_fast5_folder;

   set_map["--fastq"] = set_fastq_folder;
   set_map["-q"] = set_fastq_folder;                            op_labels["-q"] = "--fastq";
   print_map["-q"] = print_fastq_folder;
   print_help_map["-q"] = print_help_fastq_folder;
   check_map["-q"] = check_fastq_folder;

   set_map["--bam"] = set_bam_file;
   set_map["-b"] = set_bam_file;                                op_labels["-b"] = "--bam";
   print_map["-b"] = print_bam_file;
   print_help_map["-b"] = print_help_bam_file;
   check_map["-b"] = check_bam_file;

   set_map["--pos"] = set_repeat_pos;
   set_map["-p"] = set_repeat_pos;                              op_labels["-p"] = "--pos";
   print_map["-p"] = print_repeat_pos;
   print_help_map["-p"] = print_help_repeat_pos;
   check_map["-p"] = check_repeat_pos;

   set_map["--savepath"] = set_save_path;
   set_map["-s"] = set_save_path;                               op_labels["-s"] = "--savepath";
   print_map["-s"] = print_save_path;
   print_help_map["-s"] = print_help_save_path;
   check_map["-s"] = check_save_path;

   set_map["--uniq_id"] = set_uniq_id;
   set_map["-u"] = set_uniq_id;                                 op_labels["-u"] = "--uniq_id";
   print_map["-u"] = print_uniq_id;
   print_help_map["-u"] = print_help_uniq_id;
   check_map["-u"] = check_uniq_id;

   set_map["--feat_path"] = set_feat_path;
   set_map["-g"] = set_feat_path;                                 op_labels["-g"] = "--feat_path";
   print_map["-g"] = print_feat_path;
   print_help_map["-g"] = print_help_feat_path;
   check_map["-g"] = check_feat_path;

   set_map["--mod_path"] = set_mod_path;
   set_map["-n"] = set_mod_path;                                 op_labels["-n"] = "--mod_path";
   print_map["-n"] = print_mod_path;
   print_help_map["-n"] = print_help_mod_path;
   check_map["-n"] = check_mod_path;

   set_map["--epoch"] = set_epoch;
   set_map["-e"] = set_epoch;                                 op_labels["-e"] = "--epoch";
   print_map["-e"] = print_epoch;
   print_help_map["-e"] = print_help_epoch;
   check_map["-e"] = check_epoch;

   set_map["--seq_sum"] = &set_seq_sum;
   set_map["-t"] = &set_seq_sum;                                 op_labels["-t"] = "--seq_sum";
   print_map["-t"] = &print_seq_sum;
   print_help_map["-t"] = &print_help_seq_sum;
   check_map["-t"] = &check_seq_sum;

   set_map["--baseIndexPath"] = set_base_index_path;
   set_map["-j"] = set_base_index_path;                                 op_labels["-j"] = "--baseIndexPath";
   print_map["-j"] = print_base_index_path;
   print_help_map["-j"] = print_help_base_index_path;
   check_map["-j"] = check_base_index_path;

   set_map["--basecallPath"] = set_basecalled_path;
   set_map["-c"] = set_basecalled_path;                                 op_labels["-c"] = "--basecallPath";
   print_map["-c"] = print_basecalled_path;
   print_help_map["-c"] = print_help_basecalled_path;
   check_map["-c"] = check_basecalled_path;

   set_map["--fast5index"] = set_fast5_index_file;
   set_map["-i"] = set_fast5_index_file;                                 op_labels["-i"] = "--fast5index";
   print_map["-i"] = print_fast5_index_file;
   print_help_map["-i"] = print_help_fast5_index_file;
   check_map["-i"] = check_fast5_index_file;

   set_map["--repfile"] = set_repdict;
   set_map["-r"] = set_repdict;                                 op_labels["-r"] = "--repfile";
   print_map["-r"] = print_repdict;
   print_help_map["-r"] = print_help_repdict;
   check_map["-r"] = check_repdict;

   std::map<std::string, std::string>::iterator op_it;
   for (op_it=op_labels.begin(); op_it!=op_labels.begin(); op_it++){
       op_labels2[op_it->second] = op_it->first;
   }
}


