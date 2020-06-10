#ifndef COMOPTION_H_
#define COMOPTION_H_

#include <map>
#include <stdint.h>
#include <string>
#include <vector>

#include <sstream>
#include <fstream>
#include <iostream>

#include "ComFunction.h"
#include "ComStruct.h"

#define F5INDEX 0
#define GETFEAT 1
#define TRAIN 2
#define PREDICT 3
#define DETECT 4


typedef struct ComOption{
   int group_size;

   static int max_length_of_repeat_unit;
   static char* module_type[MODULE_NUM];

   RepeatRegion rep_info;
 
   bool multifast5;
   char fast5_folder[CHAR_SIZE];
   char fastq_folder[CHAR_SIZE];
   char bam_file[CHAR_SIZE];

   char save_path[CHAR_SIZE];
   char uniq_id[CHAR_SIZE];

   int submod_id;
   char submod_str[CHAR_SIZE];
   char feat_path[CHAR_SIZE];

   char mod_path[CHAR_SIZE];
   int epoch;

   char seq_sum[CHAR_SIZE];
   char fast5_index_file[CHAR_SIZE];
   
   char basecalled_path[CHAR_SIZE];
   char base_index_path[CHAR_SIZE];

   std::map<std::string, std::map<std::string, RepeatRegion> > repdict;
/*
 * 1. add variables
 * 2. add set_*, print_*, print_help_* functions below
 * 3. implement the functions in .c
 * 4. add to set_function;
 * 5. revise option list for different module in .c
 * 6. change number of option
 * */   
} ComOption;


//extern RunComOption m_runComOption;
//extern RunOption m_runOption;
extern ComOption m_comOption;


void set_submod_id(char * sub_mod);

bool check_repdict();
void print_repdict();
bool set_repdict(char* mrepfile);
void print_help_repdict();

bool check_fast5_index_file();
void print_fast5_index_file();
bool set_fast5_index_file(char* mfast5_index_file);
void print_help_fast5_index_file();

bool check_basecalled_path();
void print_basecalled_path();
bool set_basecalled_path(char* mbasecalled_path);
void print_help_basecalled_path();

bool check_seq_sum();
void print_seq_sum();
bool set_seq_sum(char* mseq_sum);
void print_help_seq_sum();

bool check_mod_path();
void print_mod_path();
bool set_mod_path(char* mmod_path);
void print_help_mod_path();

bool check_epoch();
void print_epoch();
bool set_epoch(char* mepoch);
void print_help_epoch();

bool check_feat_path();
void print_feat_path();
bool set_feat_path(char* mfeat_path);
void print_help_feat_path();

bool check_len_repeat_unit();
void print_len_repeat_unit();
bool set_len_repeat_unit(char* mlen);
void print_help_len_repeat_unit();

bool check_multifast5();
void print_multifast5();
bool set_multifast5(char* mf5);
void print_help_multifast5();

bool check_fast5_folder();
void print_fast5_folder();
bool set_fast5_folder(char* f5_folder);
void print_help_fast5_folder();

bool check_fastq_folder();
void print_fastq_folder();
bool set_fastq_folder(char* fq_folder);
void print_help_fastq_folder();

bool check_bam_file();
void print_bam_file();
bool set_bam_file(char* bam_file);
void print_help_bam_file();

bool check_repeat_pos();
void print_repeat_pos();
bool set_repeat_pos(char* pos_str);
void print_help_repeat_pos();

bool check_save_path();
void print_save_path();
bool set_save_path(char* svpath);
void print_help_save_path();

bool check_uniq_id();
void print_uniq_id();
bool set_uniq_id(char* muid);
void print_help_uniq_id();

bool check_base_index_path();
void print_base_index_path();
bool set_base_index_path(char* mbase_index_path);
void print_help_base_index_path();


void init_ComOption();
void print_options();
void set_function();
void print_help(char* m_prog, char* sub_mod);
bool set_options(int argc, char * argv[]);

#endif /* COMOPTION_H_ */
