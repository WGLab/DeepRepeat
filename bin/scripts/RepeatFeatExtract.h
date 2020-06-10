#ifndef REPEAT_FEATURE_EXTRACTION_H_
#define REPEAT_FEATURE_EXTRACTION_H_

#include "ComStruct.h"
#include "ComFunction.h"
#include "BamReader.h"
#include "Fast5Reader.h"
#include "Fast5Index.h"

#define Feat_Len 50
#define Labl_Len 4
#define Signal_Bin_Size 0.2 

#define NOT_IN_REPEAT 0
#define DEL_IN_REPEAT 1
#define INS_IN_REPEAT 2
#define NOM_IN_REPEAT 3

#define REP_BASE_TYPE_F5_BOUNDARY 0
#define REP_BASE_TYPE_NONREP   1
#define REP_BASE_TYPE_REPEAT   2
#define REP_BASE_TYPE_UNKNOWN  3
#define REP_BASE_TYPE_OTHERREP 4

class FeatSpaceOther{
public:
   std::string chrn;
   uint16_t map_strand;
   std::string read_id;
};

class FeatSpace{
public:
   //std::string chrn;
   uint64_t ref_pos;
   //uint16_t map_strand;
   uint64_t signal_num;
   double signal_mean;
   double signal_std;
   int16_t repeat_labels[Labl_Len]; // not-in-repeat, deletion, insertion, in-repeat;
   uint16_t signal_distr[Feat_Len];

   int16_t rep_base_type;   

   FeatSpace(int16_t m_rep_base_type=REP_BASE_TYPE_F5_BOUNDARY);
   ~FeatSpace(){;}
};

class RepeatFeatExtract{
   static std::map<std::string, std::map<std::string, RepeatRegion> > repeat_regions_in_genome;
   static std::map<std::string, std::string> readToF5_whole;

   BamReader * _bam_reader; 
   RunOption * _run_option;
   F5Reader * _f5_reader;

   std::string _bam_file;
   std::string _f5_folder;
   std::string _f5_index_file;
   std::string _f5_conf_file;
   bool _is_multif5; 
   double _nb_size;
   uint16_t _nb_relax_size;

   std::string _save_filename;

   int _status;
   double _gap_perc;

   bool _is_save_to_file;
   

   //int _generate_fs(std::string mchr, int16_t c_rep_label, int16_t c_rep_type, std::vector<FeatSpace>& fs_repeat, const Map1BasePos & p_mapbasepos, const F5Reader * p_f5_reader);
   int _generate_fs(const Bam1Record & b1r, int16_t c_rep_label, int16_t c_rep_type, std::vector<FeatSpace>& fs_repeat, const Map1BasePos & p_mapbasepos, const F5Reader * p_f5_reader);
   int _get_repeat_feature(std::vector<FeatSpace>& fs_repeat, std::vector<FeatSpaceOther>& fs_repeat_other);

   void _output_fs_to_file(FeatSpace & _fs_, std::ofstream & m_output_to_file);

public:
   int get_status(){ return _status; }
   static int set_repeat_regions_in_genome(std::string p_repeat_regions_in_genome_file);
   static int reset_repeat_regions_in_genome(){ repeat_regions_in_genome.clear(); return 0;}
   static int set_readToF5_whole(std::string p_readToF5_whole_file);
 
   RepeatFeatExtract(std::string p_bam_file, std::string p_f5_folder, std::string p_f5_conf_file, bool p_is_multif5, uint64_t p_min_qry_alg_len, std::string p_repeat_region_info, uint64_t p_ovlp_size, int p_rg_file=1, double p_nb_size=-5, double p_gap_perc=0.12, uint16_t p_nb_relax_size=7);
   ~RepeatFeatExtract();

   int set_f5_index_file(std::string p_f5_index_file);
   int set_nb_size(double p_nb_size);
   int reset_bam_file(std::string p_bam_file, std::string p_repeat_region_info, int p_rg_file, uint64_t p_ovlp_size);
   int reset_f5_folder(std::string p_f5_folder);
   int reset_is_multif5(std::string p_f5_conf_file, bool p_is_multif5);
   int reset_save_filename(std::string p_save_filename);

   int save_repeat_feature(std::string p_save_filename);
   std::vector<FeatSpace> get_repeat_feature(std::vector<FeatSpaceOther>& fs_repeat_other);
};

#endif
