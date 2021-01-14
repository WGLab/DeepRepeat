#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>

#include <map>
#include <string>

//#include "ComFunction.h"
//#include "ComOption.h"
//#include "Fast5Reader.h"

#include "RepeatFeatExtract.h"

//#include "H5Cpp.h"

int main (int argc, char * argv[])
{
   if (argc<3){//                             1          2        3         4         5              6             7      8
      std::cout<< "Usage: " << argv[0] << " bam-file save-file f5_folder f5_index rep_region nb_rep_region_file config_f p_nbsize" <<std::endl;
      return 1;
   }

   std::string p_bam_file(argv[1]);
   std::string p_f5_folder("data/alb231/htt12/");
   if (argc>3){ p_f5_folder = std::string(argv[3]); }
   std::string p_f5_conf_file("config/fast5_path.config");
   if (argc>7) {p_f5_conf_file = std::string(argv[7]); std::cout<< "p_f5_conf_file " << argv[7]<<std::endl; }
   bool p_is_multif5 = false;
   std::string p_save_filename(argv[2]);
   //std::string p_repeat_region_info("chr4:3074876-3074933:3");
   std::string p_repeat_region_info("chr4:3074876-3074933:CAG:3");
   if (argc>5) { p_repeat_region_info = std::string(argv[5]); }
   uint64_t p_ovlp_size = 50;
   int p_rg_file = 2;

   // p_bam_file: bam_results/   ---> sys.argv[1]
   // p_f5_folder: data/alb231/htt12/
   // p_f5_conf_file: config/fast5_path.config
   // p_is_multif5 : false
   // p_save_filename: sys.argv[2]
   // p_repeat_region_info : 
   // p_ovlp_size : 50
   // p_rg_file :  2

   std::cout<< "Input = " << p_bam_file<< " " << p_save_filename <<std::endl;
   float p_nbsize = -1.5;
   if (argc>8){
      p_nbsize = std::stof(argv[8]);
   }
   int alig_size = 100;
   if (argc>9){
      alig_size = std::stoi(argv[9]);
   }
   //RepeatFeatExtract(std::string p_bam_file, std::string p_f5_folder, std::string p_f5_conf_file, bool p_is_multif5, std::string p_save_filename, std::string p_repeat_region_info, uint64_t p_ovlp_size, int p_rg_file=1);
   RepeatFeatExtract rfe(p_bam_file, p_f5_folder, p_f5_conf_file, p_is_multif5, alig_size, p_repeat_region_info, p_ovlp_size, p_rg_file, p_nbsize);
   //rfe.set_f5_index_file("data/alb231/htt12/htt.f5index");
   if (argc>4) { rfe.set_f5_index_file(argv[4]); std::cout<< "f5index " << argv[4]<<std::endl;}
   else{ rfe.set_f5_index_file("data/alb231/htt12/htt.f5index"); }
   if (argc<6){
      RepeatFeatExtract::set_repeat_regions_in_genome("trf/wg_trf.cag.bed");
   }else{
      RepeatFeatExtract::set_repeat_regions_in_genome(argv[6]);
      std::cout<< "rep_g_region " << argv[6]<<std::endl;
   }
   
   rfe.save_repeat_feature(p_save_filename);

   return EXIT_SUCCESS;
}
