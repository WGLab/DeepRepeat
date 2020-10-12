
#include <algorithm>    // std::random_shuffle
#include <random>       // std::default_random_engine

#include "RepeatFeatExtract.h"
#include <htslib/sam.h>


std::map<std::string, std::map<std::string, RepeatRegion> > RepeatFeatExtract::repeat_regions_in_genome;
std::map<std::string, std::string> RepeatFeatExtract::readToF5_whole;


int RepeatFeatExtract::set_repeat_regions_in_genome(std::string p_repeat_regions_in_genome_file){
   BamReadOption _bam_op;
   _bam_op.set_repdict(p_repeat_regions_in_genome_file);
   repeat_regions_in_genome =  _bam_op.repdict; // _bam_op.get_repdict();
   /*for (std::map<std::string, std::map<std::string, RepeatRegion> >::iterator chr_it=repeat_regions_in_genome.begin(); chr_it!=repeat_regions_in_genome.end(); chr_it++){
      for (std::map<std::string, RepeatRegion>::iterator pos_it=chr_it->second.begin(); pos_it!=chr_it->second.end(); pos_it++){
          std::cout<<"repeat_regions_in_genome: " << chr_it->first << ":" <<pos_it->first<<":"<< pos_it->second.len_repeat_unit << std::endl;
      }
   }*/ 
   return 0;
}
int RepeatFeatExtract::set_readToF5_whole(std::string p_readToF5_whole_file){
   std::map<std::string, bool> readOfInterest;
   read_f5index(readOfInterest, readToF5_whole, p_readToF5_whole_file);
   return 0;
}

/////////////////////////////////

RepeatFeatExtract::RepeatFeatExtract(std::string p_bam_file, std::string p_f5_folder, std::string p_f5_conf_file, bool p_is_multif5, uint64_t p_min_qry_alg_len, std::string p_repeat_region_info, uint64_t p_ovlp_size, int p_rg_file, double p_nb_size, double p_gap_perc, uint16_t p_nb_relax_size){
   _bam_file = p_bam_file;
   _bam_reader = new BamReader(p_bam_file.c_str());
   if (!_bam_reader->check_bam_status()){
      std::cout<< "Error!!! Cannot open bam file="<<p_bam_file<<std::endl;
      _status |= 2;
   }else{
      if (p_rg_file==1) { 
         _bam_reader->bamReadOp.set_repdict(p_repeat_region_info);
      }else if (p_rg_file==2){
         _bam_reader->bamReadOp.set_repdict1(p_repeat_region_info, ":-");
      }
      if (p_rg_file>0){
         _bam_reader->bamReadOp.set_w_specifiedRegion(true, p_ovlp_size);
      }
      _bam_reader->bamReadOp.set_w_supplementary(true);
      _bam_reader->bamReadOp.set_w_secondary(false);
      _bam_reader->bamReadOp.set_min_read_len(p_min_qry_alg_len);
   }

   _run_option = new RunOption();
   _run_option->group_size = 4;
   _f5_reader = NULL; 
   reset_is_multif5(p_f5_conf_file, p_is_multif5);

   //reset_save_filename(p_save_filename); 
   reset_f5_folder(p_f5_folder);
   set_nb_size(p_nb_size);
   
   _gap_perc = p_gap_perc;
   _nb_relax_size = p_nb_relax_size;
   _status = 0;
}

int RepeatFeatExtract::set_f5_index_file(std::string p_f5_index_file){
   _f5_index_file = p_f5_index_file;
   std::map<std::string, bool> readOfInterest;
   read_f5index(readOfInterest, readToF5_whole, p_f5_index_file);
   return _status;
}
int RepeatFeatExtract::set_nb_size(double p_nb_size){
   if (p_nb_size<0){
      _nb_size = p_nb_size;
   }else {
      if (p_nb_size<10){
         std::cout<< "Warning!!!! neighborhood size is too small!!!" << p_nb_size << std::endl;
      }
      _nb_size = p_nb_size;
   }
   return _status;
}

int RepeatFeatExtract::reset_bam_file(std::string p_bam_file, std::string p_repeat_region_info, int p_rg_file, uint64_t p_ovlp_size){
   _status = 0;
   _bam_file = p_bam_file;
   _bam_reader->resetBam(p_bam_file.c_str());

   if (!_bam_reader->check_bam_status()){
      std::cout<< "Error!!! Cannot open bam file="<<p_bam_file<<std::endl;
      _status |= 2;
   }else{
      _bam_reader->bamReadOp.repdict.clear();
      if (p_rg_file==1) { _bam_reader->bamReadOp.set_repdict(p_repeat_region_info);}
      else if (p_rg_file==2) { _bam_reader->bamReadOp.set_repdict1(p_repeat_region_info, ":-"); }
      if (p_rg_file>0){  _bam_reader->bamReadOp.set_w_specifiedRegion(true, p_ovlp_size); }
      _bam_reader->bamReadOp.set_w_supplementary(true);
      _bam_reader->bamReadOp.set_w_secondary(false);
   }

   return _status;
}

int RepeatFeatExtract::reset_f5_folder(std::string p_f5_folder){
   _status = 0;
   _f5_folder = p_f5_folder;
   return _status;
}
int RepeatFeatExtract::reset_is_multif5(std::string p_f5_conf_file, bool p_is_multif5){
   _status = 0;
   _f5_conf_file = p_f5_conf_file;
   F5Path f5p(p_f5_conf_file.c_str());
   F5Reader::f5path = f5p;
   _is_multif5 = p_is_multif5;

   if (_f5_reader!=NULL){ delete _f5_reader; }
   _f5_reader = NULL;
   if (_is_multif5){ _f5_reader = new F5ReaderMulti(_run_option); }
   else { _f5_reader = new F5Reader1(_run_option); }

   return _status;
}
int RepeatFeatExtract::reset_save_filename(std::string p_save_filename){
   _status = 0;
   _save_filename = p_save_filename;
   return _status;
}

FeatSpace::FeatSpace(int16_t m_rep_base_type){
   ref_pos = 0;
   signal_num = 0;
   signal_mean = 0;
   signal_std = 0;
   int il;
   for (il=0; il<Labl_Len; il++){ repeat_labels[il] = 0; }
   for (il=0; il<Feat_Len; il++){ signal_distr[il] = 0; }
   rep_base_type = m_rep_base_type;
}

//int RepeatFeatExtract::_generate_fs(std::string mchr, int16_t c_rep_label, int16_t c_rep_type, std::vector<FeatSpace>& fs_repeat, const Map1BasePos & p_mapbasepos, const F5Reader * p_f5_reader){
int RepeatFeatExtract::_generate_fs(const Bam1Record & b1r, int16_t c_rep_label, int16_t c_rep_type, std::vector<FeatSpace>& fs_repeat, const Map1BasePos & p_mapbasepos, const F5Reader * p_f5_reader){
   FeatSpace fs(c_rep_type);
   fs.repeat_labels[c_rep_label] = 1;
   //if (c_rep_label<Labl_Len){
   //   fs.repeat_labels[c_rep_label] = 1;
   //}
   //fs.chrn = b1r.map_chr; //mchr;
   //fs.map_strand = b1r.map_strand;
   fs.ref_pos = p_mapbasepos.ref_pos;

   uint64_t m_qry_pos = p_mapbasepos.qry_pos;
   if (_is_multif5){ ; }
   else{
      if (m_qry_pos<2 || m_qry_pos-2>=p_f5_reader->f5event_list.size()){
         fs_repeat.push_back(fs);
         return _status;
      }
      m_qry_pos -= 2;
   }
 
   //std::cout<<" Debug " << p_mapbasepos.qry_pos << "/" <<p_f5_reader->f5event_list.size() <<std::endl; 
   //std::cout<<" Debug " << p_f5_reader->f5event_list[m_qry_pos].start << "-" << p_f5_reader->f5event_list[m_qry_pos].length << "-" << p_f5_reader->f5event_list[m_qry_pos].start+p_f5_reader->f5event_list[m_qry_pos].length <<std::endl;
   fs.signal_num = p_f5_reader->f5event_list[m_qry_pos].length;
   fs.signal_mean = p_f5_reader->f5event_list[m_qry_pos].mean;
   fs.signal_std = p_f5_reader->f5event_list[m_qry_pos].stdv;
   for (uint64_t si=p_f5_reader->f5event_list[m_qry_pos].start; si<p_f5_reader->f5event_list[m_qry_pos].start+p_f5_reader->f5event_list[m_qry_pos].length; si++){
       int16_t dis_bin = (p_f5_reader->f5signal_list[si]+5)/Signal_Bin_Size;
       if (dis_bin<0){ dis_bin = 0; }
       if (dis_bin>Feat_Len-1) { dis_bin = Feat_Len-1; }
       fs.signal_distr[ dis_bin ] += 1;
   }
   //std::cout<<" Debug 2" <<std::endl;
   if (p_f5_reader->f5event_list[m_qry_pos].length>0){
      for (int il=0; il<Feat_Len; il++){ 
         fs.signal_distr[il] = fs.signal_distr[il]*255/p_f5_reader->f5event_list[m_qry_pos].length;
      }
   }
   fs_repeat.push_back(fs);

   /*
   std::cout<<"fs = "<<fs.chrn<<":"<<fs.ref_pos<<"/"<<m_qry_pos<<":" << p_f5_reader->f5event_list[m_qry_pos].model_state[2]<<"<";
   for (int il=0; il<Labl_Len; il++){ std::cout << " "<<fs.repeat_labels[il]; }
   std::cout<< " > tl=" << fs.rep_base_type << " " << c_rep_label;
   std::cout<<" " << fs.signal_mean << "+" << fs.signal_std << "/" << fs.signal_num << std::endl;
   */

   return _status;
}

int RepeatFeatExtract::_get_repeat_feature(std::vector<FeatSpace>& fs_repeat, std::vector<FeatSpaceOther>& fs_repeat_other){
   RepeatFeatExtract::repeat_regions_in_genome.insert(_bam_reader->bamReadOp.repdict.begin(), _bam_reader->bamReadOp.repdict.end());

   int32_t m_nb_size = 0; // _nb_size
   std::ofstream m_output_to_file;
   std::ofstream m_output_other_to_file;
   std::vector<FeatSpace> _this_fs_;

   if (_is_save_to_file){
      m_output_to_file.open(_save_filename.c_str(), std::ofstream::out);
      m_output_other_to_file.open((_save_filename+".more").c_str(), std::ofstream::out);
   }

   /* comment so that read record one-by-one
   std::cout<< "Read bam " << std::endl;
   int m_reader_status = _bam_reader->readBam();
   if (m_reader_status!=0){
      std::cout<< "Error Cannot read bam from " << _bam_file <<std::endl;
      _status |= 4;
      return _status;
   }
   std::cout<< "Read bam done " << _bam_reader->br_list.size() << " _status=" << _status << std::endl;
   
   int try_times = 1;
   int r_i = 0;
   for (; try_times>0; try_times--,r_i++){
       ;//std::cout<< _bam_reader->Basic_Bam1Record_toString(_bam_reader->br_list[r_i]); 
   }
   */

   if (!_is_save_to_file){
      fs_repeat.clear();
      fs_repeat_other.clear();
      // comment for reverse alignment
      //fs_repeat.push_back(FeatSpace());
   }
   std::map<uint64_t, bool> repeat_pos_info;
   std::map<uint64_t, bool> repeat_pos_dup_info;
   std::map<uint64_t, bool>::iterator rep_pos_inf_it;
   std::map<std::string, std::string>::iterator r_to_f5_it;
   /* only when read all records at the same time
   //std::random_shuffle(_bam_reader->br_list.begin(), _bam_reader->br_list.end());
   std::shuffle(_bam_reader->br_list.begin(), _bam_reader->br_list.end(), std::default_random_engine(7));
   */
   bool c_output_fs = false;
   int r_i=0;
   bool _is_first = true;
   FeatSpace _empty_fs_;

   /* for read all records once
   for(r_i=0; r_i<_bam_reader->br_list.size(); r_i++){// for each alignment;
   */
   while(_bam_reader->read1RecordFromBam()==0){
   //for(r_i=0; r_i<3; r_i++){
      //if (_is_save_to_file){
         _this_fs_.clear();
      //}
      _bam_reader->br_list.clear();
      _bam_reader->br_list.push_back(_bam_reader->br);

      std::vector<RepeatRegion> c_rep_regions_in_whole;
      std::cout<< "For: " << _bam_reader->Basic_Bam1Record_toString(_bam_reader->br_list[r_i]) << std::endl;
      _bam_reader->bamReadOp.check_specifiedRegion_multi(_bam_reader->br_list[r_i].map_chr.c_str(), _bam_reader->br_list[r_i].ref_start_pos, _bam_reader->br_list[r_i].ref_end_pos, c_rep_regions_in_whole, RepeatFeatExtract::repeat_regions_in_genome); // get all repeat regions in this alignment;
      if (c_rep_regions_in_whole.size()==0){
         //std::cout<< "Error! cannot find repeat regions for whole genome in --"<< RepeatFeatExtract::repeat_regions_in_genome.size() << "--" <<  _bam_reader->Basic_Bam1Record_toString(_bam_reader->br_list[r_i]);
         continue;
      }
      for (int crri=0; crri<c_rep_regions_in_whole.size(); crri++){
         std::cout<< "c_rep_regions_in_whole_" << crri << " "  << c_rep_regions_in_whole[crri].chrn << ":"<<c_rep_regions_in_whole[crri].start_pos <<"-"<<c_rep_regions_in_whole[crri].end_pos<<":"<<c_rep_regions_in_whole[crri].len_repeat_unit<<std::endl;
      }
      std::vector<RepeatRegion> c_rep_regions;
      //_bam_reader->bamReadOp.check_specifiedRegion_multi(_bam_reader->br_list[r_i].map_chr.c_str(), _bam_reader->br_list[r_i].ref_start_pos, _bam_reader->br_list[r_i].ref_end_pos, c_rep_regions, _bam_reader->bamReadOp.get_repdict()); // get all repeat regions in this alignment;
      _bam_reader->bamReadOp.check_specifiedRegion_multi(_bam_reader->br_list[r_i].map_chr.c_str(), _bam_reader->br_list[r_i].ref_start_pos, _bam_reader->br_list[r_i].ref_end_pos, c_rep_regions, _bam_reader->bamReadOp.repdict);
      if (c_rep_regions.size()==0){
         //std::cout<< "Error! cannot find repeat regions in --"<< _bam_reader->bamReadOp.repdict.size() << "--" <<  _bam_reader->Basic_Bam1Record_toString(_bam_reader->br_list[r_i]);
         continue;
      }
      for (int crri=0; crri<c_rep_regions.size(); crri++){
         std::cout<< "c_rep_regions_" << crri << " " << c_rep_regions[crri].chrn << ":"<<c_rep_regions[crri].start_pos <<"-"<<c_rep_regions[crri].end_pos<<":"<<c_rep_regions[crri].len_repeat_unit<<std::endl;
      }
      repeat_pos_info.clear();
      repeat_pos_dup_info.clear();
      for (std::vector<RepeatRegion>::iterator crr_it=c_rep_regions_in_whole.begin(); crr_it!=c_rep_regions_in_whole.end(); crr_it++){
          int m_nb_size = crr_it->end_pos - crr_it->start_pos;
          if (_nb_size<0) { m_nb_size = m_nb_size*(-_nb_size);
          }else { m_nb_size = _nb_size; }

          uint64_t start_m = (crr_it->start_pos>m_nb_size)? crr_it->start_pos-m_nb_size : 0;
          uint64_t end_m = crr_it->end_pos + m_nb_size;
          for (uint64_t mi=start_m; mi<end_m; mi++){
              bool is_in_repeat = mi>=crr_it->start_pos;
              is_in_repeat = is_in_repeat && (mi<crr_it->end_pos);
              rep_pos_inf_it = repeat_pos_info.find(mi);
              if (rep_pos_inf_it != repeat_pos_info.end()){
                 is_in_repeat = is_in_repeat || (rep_pos_inf_it->second); 
                 repeat_pos_dup_info[mi] = true;
              }
              repeat_pos_info[mi] = is_in_repeat;
          }
      }

      r_to_f5_it = readToF5_whole.find(_bam_reader->br_list[r_i].qry_name); // get fast5;
      if (r_to_f5_it == readToF5_whole.end()){ 
          std::cout<< "Error! Cannot find read_it= " << _bam_reader->br_list[r_i].qry_name << std::endl;
          continue;
      }
     
      _f5_reader->resetFast5File((_f5_folder+"/"+(r_to_f5_it->second)).c_str()); // read fast5: associate signals to fq
     std::cout<< "Read fast5 file = " <<  _f5_folder+"/"+(r_to_f5_it->second)<<std::endl;
     if ( _f5_reader->readF5data() != 0){ // f5event_list;      f5signal_list
        std::cout<< "Error! cannot read fast5 file= " << _f5_folder+"/"+(r_to_f5_it->second) << std::endl;
        continue;
     } 

     uint64_t pos_it;
     double left_total, right_total;
     double left_gap, right_gap;

     for (std::vector<RepeatRegion>::iterator crr_it=c_rep_regions.begin(); crr_it!=c_rep_regions.end(); crr_it++){
        c_output_fs = false;

        int m_nb_size = crr_it->end_pos - crr_it->start_pos;
        if (_nb_size<0) { m_nb_size = m_nb_size*(-_nb_size);
        }else { m_nb_size = _nb_size; }
        int m_nb_relax_size = crr_it->len_repeat_unit * 2 + 1;

        std::cout << "Repeat Info: "<< crr_it->chrn<<":"<<crr_it->start_pos<<"-"<<crr_it->end_pos<<":"<<crr_it->len_repeat_unit<< " nb=" << m_nb_size<<" r_nb=" << m_nb_relax_size << " f5_ev_size=" << _f5_reader->f5event_list.size() << std::endl;
        uint64_t start_m = (crr_it->start_pos>m_nb_size)? crr_it->start_pos-m_nb_size : 0;
        uint64_t end_m = crr_it->end_pos + m_nb_size;

        left_total = 0;
        right_total = 0;
        left_gap = 0;
        right_gap = 0;
        for (pos_it=0; pos_it<_bam_reader->br_list[r_i].map_pos_detail.size(); pos_it++){
           if (_bam_reader->br_list[r_i].map_pos_detail[pos_it].ref_pos < start_m) { continue; }
           if (_bam_reader->br_list[r_i].map_pos_detail[pos_it].ref_pos >= end_m) { break; }
           rep_pos_inf_it = repeat_pos_dup_info.find(_bam_reader->br_list[r_i].map_pos_detail[pos_it].ref_pos);
           if (_bam_reader->br_list[r_i].map_pos_detail[pos_it].ref_pos >= start_m && _bam_reader->br_list[r_i].map_pos_detail[pos_it].ref_pos<crr_it->start_pos-m_nb_relax_size){
              if (rep_pos_inf_it==repeat_pos_dup_info.end() || !(repeat_pos_info[_bam_reader->br_list[r_i].map_pos_detail[pos_it].ref_pos])){
                  left_total += 1;
                  if ((_bam_reader->br_list[r_i].map_pos_detail[pos_it].map_type & 15) != BAM_CMATCH && (_bam_reader->br_list[r_i].map_pos_detail[pos_it].map_type & 15) != BAM_CEQUAL && (_bam_reader->br_list[r_i].map_pos_detail[pos_it].map_type & 15) != BAM_CDIFF ){
                      left_gap += 1;
                  }
              }
           }else if (_bam_reader->br_list[r_i].map_pos_detail[pos_it].ref_pos >= crr_it->end_pos+m_nb_relax_size && _bam_reader->br_list[r_i].map_pos_detail[pos_it].ref_pos < end_m){
              if (rep_pos_inf_it==repeat_pos_dup_info.end() || !(repeat_pos_info[_bam_reader->br_list[r_i].map_pos_detail[pos_it].ref_pos])){
                  right_total += 1;
                  if ((_bam_reader->br_list[r_i].map_pos_detail[pos_it].map_type & 15) != BAM_CMATCH && (_bam_reader->br_list[r_i].map_pos_detail[pos_it].map_type & 15) != BAM_CEQUAL && (_bam_reader->br_list[r_i].map_pos_detail[pos_it].map_type & 15) != BAM_CDIFF ){
                     right_gap += 1;
                  }
              }
           }  
        }
           
        if (left_total>0){ left_gap /= left_total; }
        if (right_total>0){ right_gap /= right_total; }
          
        std::cout << " nb_info L=" << left_total << "/" << left_gap << " R=" << right_total << "/" << right_gap << std::endl;
 
        int16_t c_rep_label;
        int16_t c_rep_type;
        bool save_id = false;
        for (pos_it=0; pos_it<_bam_reader->br_list[r_i].map_pos_detail.size(); pos_it++){
           if (_bam_reader->br_list[r_i].map_pos_detail[pos_it].ref_pos < start_m) { continue; }
           if (_bam_reader->br_list[r_i].map_pos_detail[pos_it].ref_pos >= end_m) { break; }
           if ((_bam_reader->br_list[r_i].map_pos_detail[pos_it].map_type & 15) == BAM_CDEL ) { continue; }
           //if ((_bam_reader->br_list[r_i].map_pos_detail[pos_it].map_type & 15) == BAM_CSOFT_CLIP ) { continue; }
           //if ((_bam_reader->br_list[r_i].map_pos_detail[pos_it].map_type & 15) == BAM_CHARD_CLIP ) { continue; }
           if ((_bam_reader->br_list[r_i].map_pos_detail[pos_it].map_type & 15) == BAM_CPAD ) { continue; }

           c_rep_label = Labl_Len + 1;
           c_rep_type = REP_BASE_TYPE_UNKNOWN;
           rep_pos_inf_it = repeat_pos_dup_info.find(_bam_reader->br_list[r_i].map_pos_detail[pos_it].ref_pos);
           if (_bam_reader->br_list[r_i].map_pos_detail[pos_it].ref_pos >= start_m && _bam_reader->br_list[r_i].map_pos_detail[pos_it].ref_pos <crr_it->start_pos) {
              if (rep_pos_inf_it==repeat_pos_dup_info.end() || !(repeat_pos_info[_bam_reader->br_list[r_i].map_pos_detail[pos_it].ref_pos]) ) {
                  c_rep_label = NOT_IN_REPEAT; 
                  if (left_gap<_gap_perc && left_total>m_nb_size/2){ c_rep_type = REP_BASE_TYPE_NONREP; }
                  else { c_rep_type = REP_BASE_TYPE_UNKNOWN; }
              } else { 
                  c_rep_label = NOM_IN_REPEAT; 
                  if (left_gap<_gap_perc && left_total>m_nb_size/2){ c_rep_type = REP_BASE_TYPE_OTHERREP; }
                  else { c_rep_type = REP_BASE_TYPE_UNKNOWN; }
              }
           }else if (_bam_reader->br_list[r_i].map_pos_detail[pos_it].ref_pos >= crr_it->end_pos && _bam_reader->br_list[r_i].map_pos_detail[pos_it].ref_pos < end_m) {
              if (rep_pos_inf_it==repeat_pos_dup_info.end() || !(repeat_pos_info[_bam_reader->br_list[r_i].map_pos_detail[pos_it].ref_pos]) ){
                 c_rep_label = NOT_IN_REPEAT;
                 if (right_gap<_gap_perc && right_total>m_nb_size/2){ c_rep_type = REP_BASE_TYPE_NONREP; }
                 else { c_rep_type = REP_BASE_TYPE_UNKNOWN;; }
              } else {
                 c_rep_label = NOM_IN_REPEAT;
                 if (right_gap<_gap_perc && right_total>m_nb_size/2){ c_rep_type = REP_BASE_TYPE_OTHERREP; }
                 else { c_rep_type = REP_BASE_TYPE_UNKNOWN; }
              }
           }else if (_bam_reader->br_list[r_i].map_pos_detail[pos_it].ref_pos >= crr_it->start_pos && _bam_reader->br_list[r_i].map_pos_detail[pos_it].ref_pos < crr_it->end_pos){ // - crr_it->len_repeat_unit){
              c_rep_label = NOM_IN_REPEAT;
              c_rep_type = REP_BASE_TYPE_REPEAT;
           }

           uint64_t cur_op_len = ((_bam_reader->br_list[r_i].map_pos_detail[pos_it].map_type)>>4);
           if ((_bam_reader->br_list[r_i].map_pos_detail[pos_it].map_type & 15) == BAM_CMATCH || (_bam_reader->br_list[r_i].map_pos_detail[pos_it].map_type & 15) == BAM_CEQUAL || (_bam_reader->br_list[r_i].map_pos_detail[pos_it].map_type & 15) == BAM_CDIFF ) {
              if (pos_it>1 && c_rep_type==REP_BASE_TYPE_REPEAT){
                 if ((_bam_reader->br_list[r_i].map_pos_detail[pos_it-1].map_type & 15) == BAM_CDEL ) {
                    cur_op_len = ((_bam_reader->br_list[r_i].map_pos_detail[pos_it-1].map_type)>>4);
                    if (cur_op_len>(crr_it->len_repeat_unit>3?crr_it->len_repeat_unit/2:2)){
                       c_rep_type = REP_BASE_TYPE_UNKNOWN;
                    } else{
                       c_rep_label = DEL_IN_REPEAT;
                    }
                 }
              }
           //}else if ((_bam_reader->br_list[r_i].map_pos_detail[pos_it].map_type & 15) == BAM_CDEL ) { ;
           }else if ((_bam_reader->br_list[r_i].map_pos_detail[pos_it].map_type & 15) == BAM_CINS || (_bam_reader->br_list[r_i].map_pos_detail[pos_it].map_type & 15) == BAM_CSOFT_CLIP || (_bam_reader->br_list[r_i].map_pos_detail[pos_it].map_type & 15) == BAM_CHARD_CLIP ) {
              if (cur_op_len>(crr_it->len_repeat_unit>3?crr_it->len_repeat_unit/2:2)){ 
                 if (c_rep_type==REP_BASE_TYPE_REPEAT) { c_rep_type = REP_BASE_TYPE_UNKNOWN; } 
              } else{
                 if (c_rep_type==REP_BASE_TYPE_REPEAT) { c_rep_label = INS_IN_REPEAT; }
              }
           }else{
              std::cout<< "Warning!!! not type of flag " << (_bam_reader->br_list[r_i].map_pos_detail[pos_it].map_type & 15) << std::endl;
           }
           //int m_sts = _generate_fs(_bam_reader->br_list[r_i].map_chr, c_rep_label, c_rep_type, fs_repeat, _bam_reader->br_list[r_i].map_pos_detail[pos_it], _f5_reader); 
           /* revised with read record one by one; comment for revserse alignment
              int m_sts = _generate_fs(_bam_reader->br_list[r_i], c_rep_label, c_rep_type, fs_repeat, _bam_reader->br_list[r_i].map_pos_detail[pos_it], _f5_reader);
           */
           int m_sts = _generate_fs(_bam_reader->br_list[r_i], c_rep_label, c_rep_type, _this_fs_, _bam_reader->br_list[r_i].map_pos_detail[pos_it], _f5_reader);
           if (!save_id){
              if (_is_save_to_file){
                 m_output_other_to_file << _bam_reader->br_list[r_i].map_chr << " " << _bam_reader->br_list[r_i].map_strand << " "<< _bam_reader->br_list[r_i].qry_name << "\n";
              }else{
                 FeatSpaceOther fso;
                 fso.chrn = _bam_reader->br_list[r_i].map_chr;
                 fso.map_strand = _bam_reader->br_list[r_i].map_strand;
                 fso.read_id = _bam_reader->br_list[r_i].qry_name;
                 fs_repeat_other.push_back(fso);
              }
              save_id = true;
              c_output_fs = true;;
           }
           if (m_sts>0){
              std::cout<< "Warning!!! cannot get feature for " <<r_i<<"/" <<pos_it<< " status=" << m_sts << "/" << _status << "\n"; 
              std::cout<<"\t"<<_bam_reader->Basic_Bam1Record_toString(_bam_reader->br_list[r_i]);
              std::cout<<"\t" << c_rep_type << " " <<  c_rep_label << " " << (_bam_reader->br_list[r_i].map_pos_detail[pos_it].map_type & 15) << " " << (_bam_reader->br_list[r_i].map_pos_detail[pos_it].map_type >> 4) << " " << _bam_reader->br_list[r_i].map_pos_detail[pos_it].ref_pos << "/" <<_bam_reader->br_list[r_i].map_pos_detail[pos_it].qry_pos << std::endl;
           }
        }
        if (c_output_fs){
            if (_is_save_to_file){
               if (_is_first){
                  _output_fs_to_file(_empty_fs_, m_output_to_file);
                  _is_first = false;
               }
               if (_bam_reader->br_list[r_i].map_strand==0){
                  for(std::vector<FeatSpace>::iterator _t_it=_this_fs_.begin(); _t_it!=_this_fs_.end(); _t_it++){
                      _output_fs_to_file(*_t_it, m_output_to_file);
                  }
               }else{
                  for(std::vector<FeatSpace>::reverse_iterator _t_rit=_this_fs_.rbegin(); _t_rit!=_this_fs_.rend(); _t_rit++){
                      _output_fs_to_file(*_t_rit, m_output_to_file);
                  }
               }
               _output_fs_to_file(_empty_fs_, m_output_to_file);
            }else{
               if (_is_first){
                   fs_repeat.push_back(FeatSpace());
                   _is_first = false;
               }               
               if (_bam_reader->br_list[r_i].map_strand==0){
                  fs_repeat.insert(fs_repeat.end(), _this_fs_.begin(), _this_fs_.end());
               }else{
                  fs_repeat.insert(fs_repeat.end(), _this_fs_.rbegin(), _this_fs_.rend()); 
               }
               fs_repeat.push_back(FeatSpace());
            }
            //comment for revserse alignment
            //fs_repeat.push_back(FeatSpace());
        }else{
            std::cout<< "Warning!!! No output feature for " <<r_i<<"/" <<crr_it->chrn<<":"<<crr_it->start_pos<<"-"<<crr_it->end_pos<< " " << _f5_folder+"/"+(r_to_f5_it->second) << " " << _bam_reader->Basic_Bam1Record_toString(_bam_reader->br_list[r_i]) << std::endl;
        }
     } 
     //fs_repeat.push_back(FeatSpace()); 
   }

   if (_is_save_to_file){
      m_output_to_file.close();
      m_output_other_to_file.close();
   }

   return _status;
}

void RepeatFeatExtract::_output_fs_to_file(FeatSpace & _fs_, std::ofstream & m_output_to_file){
      int il;
      m_output_to_file << _fs_.ref_pos;
      m_output_to_file << " " << _fs_.rep_base_type;
      for (il=0; il<Labl_Len; il++){
         m_output_to_file << " " << _fs_.repeat_labels[il];
      }
      m_output_to_file << " " << _fs_.signal_num<< " " << _fs_.signal_mean<<" " <<_fs_.signal_std;
      for(il=0; il<Feat_Len; il++){
         m_output_to_file << " " << _fs_.signal_distr[il];
      }
      m_output_to_file << "\n";
}

int RepeatFeatExtract::save_repeat_feature(std::string p_save_filename){
   reset_save_filename(p_save_filename);
   std::vector<FeatSpace> fs_repeat;
   std::vector<FeatSpaceOther> fs_repeat_other;
   _is_save_to_file = true;
 
   if (_get_repeat_feature(fs_repeat, fs_repeat_other)!=0){
      std::cout<<"find repeat error " << fs_repeat.size()<<std::endl;
   }

   return _status;

   std::ofstream m_output_to_file(p_save_filename.c_str(), std::ofstream::out);
   std::ofstream m_output_other_to_file((p_save_filename+".more").c_str(), std::ofstream::out);
   int il;
   bool new_file = false; 
   std::vector<FeatSpaceOther>::iterator fso_it = fs_repeat_other.begin();
   for(std::vector<FeatSpace>::iterator fs_it=fs_repeat.begin(); fs_it!=fs_repeat.end(); fs_it++){
      if (fs_it->rep_base_type==REP_BASE_TYPE_F5_BOUNDARY){
         new_file = true;
      }
      if (new_file && fs_it->rep_base_type!=REP_BASE_TYPE_F5_BOUNDARY){
         new_file = false;
         m_output_other_to_file << fso_it->chrn << " " << fso_it->map_strand << " "<< fso_it->read_id << "\n";
         fso_it++;
      }

      m_output_to_file << fs_it->ref_pos;
      m_output_to_file << " " << fs_it->rep_base_type;
      for (il=0; il<Labl_Len; il++){
         m_output_to_file << " " << fs_it->repeat_labels[il];
      }
      m_output_to_file << " " << fs_it->signal_num<< " " << fs_it->signal_mean<<" " <<fs_it->signal_std;
      for(il=0; il<Feat_Len; il++){
         m_output_to_file << " " << fs_it->signal_distr[il];
      }  
      m_output_to_file << "\n";
   }
   m_output_to_file.close();
   m_output_other_to_file.close();

   return _status;
}

std::vector<FeatSpace> RepeatFeatExtract::get_repeat_feature(std::vector<FeatSpaceOther>& fs_repeat_other){
   std::vector<FeatSpace> fs_repeat;
   //std::vector<FeatSpaceOther> fs_repeat_other;
   _is_save_to_file = false;
   
   if (_get_repeat_feature(fs_repeat, fs_repeat_other)!=0){
      std::cout<<"find repeat error " << fs_repeat.size()<<std::endl;
   }

   return fs_repeat;
}

RepeatFeatExtract::~RepeatFeatExtract(){
  if (_bam_reader!=NULL){ delete _bam_reader; }
  if (_run_option!=NULL){ delete _run_option; }
  if (_f5_reader!=NULL){ delete _f5_reader; }
}


