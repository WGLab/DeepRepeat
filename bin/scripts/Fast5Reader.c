#include "Fast5Reader.h"

#include <iostream>     // std::cout, std::right, std::endl
#include <fstream>
#include <sstream>
#include <math.h>
#include <algorithm>

#include "ComFunction.h"

int F5Path::read_config(const char * path_config_file){
   //std::cout<< " read_config " << path_config_file <<std::endl;
   int error_num = 0;
   std::ifstream m_in(path_config_file);
   if (!(m_in.good())) { 
      std::cout<< " Cannot open file "<< path_config_file << std::endl;
      status = 1;
      return error_num+1;
   }

   std::string line;
   while (std::getline(m_in, line)){
      if (line.size()==0){continue;}
      if (line[0]=='#'){continue;}
      //std::cout<< "\t read " << line<<std::endl;
      size_t first_eq_sign = line.find_first_of('=');
      std::string path_name = line.substr(0,first_eq_sign);
      std::string path_value = line.substr(first_eq_sign+1);
      if (path_value.size()==0){ error_num += 1; }
      if (path_name.compare("Fastq_path")==0){
          fastq_path = path_value;
      }else if (path_name.compare("Signal_path")==0){
          signal_path = path_value;
      }else if (path_name.compare("Event_path")==0){
          event_path = path_value;
      }else if (path_name.compare("Channel_path")==0){
          channel_path = path_value;
      }else {
         std::cout<<"Warning!!! '"<<path_name<<"' is not supported."<<std::endl;
      }
   }
   m_in.close();
   status = 0;
   return error_num;
}

F5Path::F5Path(const char * path_config_file){
   read_config(path_config_file);
}

int F5Path::get_status(){
   if (status==0){
      if (fastq_path.size()==0) { status = 2;}
      else if (event_path.size()==0) { status = 3;}
      else if (signal_path.size()==0) { status = 4;}
      else if (channel_path.size()==0) { status = 5; }
   }
   return status;
}

std::string F5Path::digiti_str = "digitisation";
std::string F5Path::offset_str = "offset";
std::string F5Path::range_str  = "range";
std::string F5Path::samplr_str = "sampling_rate";

std::string F5Path::mean_str = "mean";
std::string F5Path::stdv_str = "stdv";
std::string F5Path::start_str = "start";
std::string F5Path::length_str = "length";
std::string F5Path::model_state_str = "model_state";
std::string F5Path::move_str = "move";

/*F5Path::F5Path(const F5Path & m_f5p){
   fastq_path = m_f5p.fastq_path;
   event_path = m_f5p.event_path;
   signal_path = m_f5p.signal_path;
   channel_path = m_f5p.channel_path;
}*/

void F5Path::toPrint(){
   if (get_status()==0){
      std::cout<< "\t Fastq_path = " << fastq_path << std::endl;
      std::cout<< "\t Event_path = " << event_path << std::endl;
      std::cout<< "\t Signal_path = " << signal_path << std::endl;
      std::cout<< "\t Channel_path = " << channel_path << std::endl;
   }else {
      std::cout<< "Warning!! cannot Print! Error happend in read! " << status << std::endl;
   }
   //std::cout<< "Output done!" << std::endl;
}
std::string F5Path::get_fastq_path(){ return fastq_path; }
std::string F5Path::get_event_path(){ return event_path; }
std::string F5Path::get_channel_path(){ return channel_path; }
std::string F5Path::get_signal_path(){ return signal_path; }

///////////////////////////////////////////////////////////////
//
F5Channel::F5Channel(){
   digitisation = -1;
   offset = -1;
   range = -1;
   sampling_rate = -1;

   status = 0;
}
int F5Channel::reset(){
   status = 0;
   return status;
}
int F5Channel::get_status(){
   if (status==0){
      if (digitisation==-1) {status +=1;}
      if (offset==-1){status +=1;}
      if (range==-1){status +=1;}
      if (sampling_rate==-1){status +=1;}
   }

   return status;
}

void F5Channel::toPrint(){
    if (get_status()==0){
        std::cout << "digitisation=" << digitisation << " offset=" << offset << " range=" << range << " sampling_rate=" <<sampling_rate << std::endl;
    }else {
        std::cout<< "Warning!! cannot Print! Error happend in read! " << status << std::endl;
   }
}

int F5Fastq::reset(){
   name_line.clear();
   fq_name.clear();
   fq.clear();
   qual.clear();

   return 0;
}

///////////////////////////////////////////////////////////////

F5Path F5Reader::f5path;

std::string F5Event_toString(const F5Event & m_f5event){
   std::ostringstream oss_tostr;
   
   oss_tostr << "Signals from " << m_f5event.start << " with length=" << m_f5event.length << ": mean="<<m_f5event.mean << " std=" << m_f5event.stdv << " for "<< m_f5event.model_state << " move="<< m_f5event.move; // << "\n";

   return oss_tostr.str();
}

std::string F5Event_toString_basic(const F5Event & m_f5event){
   std::ostringstream oss_tostr;

   oss_tostr << m_f5event.start << "+" << m_f5event.length << " " << std::setprecision (5) << m_f5event.mean<<"+" << m_f5event.stdv << " " << m_f5event.model_state << "_" << m_f5event.move;

   return oss_tostr.str();
}

std::string F5Event_toString(const F5Event & m_f5event, const std::vector<double>& f5s){
   std::ostringstream oss_tostr;

   oss_tostr << m_f5event.start << "+" << m_f5event.length << " " << std::setprecision (5) << m_f5event.mean<<"+"<<std::setprecision (5) << m_f5event.stdv << " " << m_f5event.model_state << "_" << m_f5event.move;
   for(uint64_t fei=m_f5event.start; fei<m_f5event.start+m_f5event.length; fei++){
      oss_tostr << " " << std::setprecision (5) << f5s[fei];
   }
   
   return oss_tostr.str();
}


int F5Reader::readF5data(){
   int error_num = 0;

   f5event_list.clear();
   f5signal_list.clear();

   if (m_runOption_pointer->group_size<1){ m_runOption_pointer->group_size = 4; }

   error_num += readF5Channel();
   if (error_num == 0){
      error_num += readF5Fastq();
   }

   if (error_num == 0){
      error_num += readF5Signals();
   }
   if (error_num == 0){
      error_num += normalize_signal();
   }

   if (error_num == 0){
      error_num += readF5Events();
   }
   if (error_num == 0){
      error_num += associate_signal_to_fq();   
   }

   return error_num;
}

int F5Reader::_init_(const char* m_f5file, RunOption *op_runOption_pointer){
   status = 0;
   f5file = m_f5file;
   if (op_runOption_pointer!=NULL){
      m_runOption_pointer = op_runOption_pointer;
   }
   if (m_runOption_pointer==NULL){
      std::cout<< "Error!!! RunOption_pointer is NULL" << std::endl;
      status = 41;
   }
   try{
       H5::Exception::dontPrint();
       m_f5_file_ptr = new H5::H5File(f5file, H5F_ACC_RDONLY);
   }catch( H5::FileIException error )
   {
      std::cout<< "Error!!! Canot open file "<< f5file << std::endl;
      error.printError();
      status = 1;
   }
}
int F5Reader::_del_(){
   if (m_f5_file_ptr!=NULL){
      delete m_f5_file_ptr;
   }
}

F5Reader::F5Reader(RunOption *op_runOption_pointer){
   m_f5_file_ptr = NULL;
   m_runOption_pointer = op_runOption_pointer;
   if (m_runOption_pointer==NULL){
      std::cout<< "Error!!! RunOption_pointer is NULL" << std::endl;
      status = 41;
   }
}

F5Reader::F5Reader(const char* m_f5file, RunOption *op_runOption_pointer){
   m_f5_file_ptr = NULL;
   _init_(m_f5file, op_runOption_pointer);
}
int F5Reader::resetFast5File(const char* m_f5file, RunOption *op_runOption_pointer){
   _del_();
   _init_(m_f5file, op_runOption_pointer);
   f5event_list.clear();
   f5signal_list.clear();

   f5fq.reset();
   f5channel.reset();
}
 
F5Reader::~F5Reader(){
   _del_();
}

//////////
int F5ReaderMulti::readF5Channel(){
   _readF5Channel_(m_runOption_pointer->read_id+"/"+f5path.get_channel_path());
   return status;
}

int F5ReaderMulti::readF5Fastq(){
   _readF5Fastq_(m_runOption_pointer->read_id+"/"+f5path.get_fastq_path());
   return status;
}

int F5ReaderMulti::readF5Events(){
   _readF5Events_(m_runOption_pointer->read_id+"/"+f5path.get_event_path());
   return status;
}

int F5ReaderMulti::readF5Signals(){
   _readF5Signals_(m_runOption_pointer->read_id+"/"+f5path.get_signal_path());
   return status;
}

///////////////////////////
int F5Reader1::readF5Channel(){
   _readF5Channel_(f5path.get_channel_path());
   return status;
}

int F5Reader1::readF5Fastq(){
   _readF5Fastq_(f5path.get_fastq_path());
   return status;
}

int F5Reader1::readF5Events(){
   _readF5Events_(f5path.get_event_path());
   return status;
}

int F5Reader1::readF5Signals(){
   try{
       H5::Exception::dontPrint();
       H5::Group signal_g = m_f5_file_ptr->openGroup(f5path.get_signal_path());

       if (signal_g.getNumObjs()==0){
          status = 21;
          std::cout<< " No signal group for "<< f5file << std::endl;
       }else{
          m_runOption_pointer->read_num = signal_g.getObjnameByIdx(0);
          _readF5Signals_(f5path.get_signal_path()+"/"+m_runOption_pointer->read_num+"/Signal");
       }
       //std::cout<< "In readF5Signals "<< signal_g.getNumObjs() << " " <<signal_g.getObjnameByIdx(0) << std::endl;
   }catch( H5::GroupIException error )
   {
      std::cout << "Error!!! Cannot open channel group from " << f5file << std::endl;
      error.printError();
      status = 22;
   }
   return status;
}


///////////////////
int F5Reader::_readF5Channel_(std::string f5c_str){
   try{
       H5::Exception::dontPrint();
       H5::Group chan_g = m_f5_file_ptr->openGroup(f5c_str);

       H5::Attribute m_at = chan_g.openAttribute(f5path.get_digiti_str());
       H5::DataType m_dt = m_at.getDataType();
       m_at.read(m_dt, &(f5channel.digitisation));

       m_at = chan_g.openAttribute(f5path.get_offset_str());
       m_dt = m_at.getDataType();
       m_at.read(m_dt, &(f5channel.offset));

       m_at = chan_g.openAttribute(f5path.get_range_str());
       m_dt = m_at.getDataType();
       m_at.read(m_dt, &(f5channel.range));

       m_at = chan_g.openAttribute(f5path.get_samplr_str());
       m_dt = m_at.getDataType();
       m_at.read(m_dt, &(f5channel.sampling_rate));

       f5channel.get_status();

       chan_g.close();
   }catch( H5::FileIException error )
   {
      std::cout<< "Error!!! Canot open file "<< f5file << std::endl;
      error.printError();
      status = 1;
   }catch( H5::GroupIException error )
   {
      std::cout << "Error!!! Cannot open channel group from " << f5file << std::endl;
      error.printError();
      status = 11;
   }


   return status;
}

int F5Reader::_readF5Fastq_(std::string f5f_str){
   if (status!=0){ return status; }

   m_runOption_pointer->read_num.clear();
   try{
       H5::Exception::dontPrint();
       H5::DataSet fq_ds = m_f5_file_ptr->openDataSet(f5f_str);
       f5fq.reset();

       if (fq_ds.getStrType().getCset()!=H5T_CSET_ASCII){
          std::cout<< "Warning!!! fq in " << f5file << " not encoded by ASCII." << std::endl;
       }
       H5::DataType mdatatype= fq_ds.getDataType();
       fq_ds.read(f5fq.name_line, mdatatype);
      
       std::string subelemt;
       std::istringstream iss(f5fq.name_line);
       if (!std::getline(iss, f5fq.fq_name, '\n')){
          std::cout<< "Error!!! cannot get fq_name line from " << f5file << std::endl;
          status = 3;
       }else if (!std::getline(iss, f5fq.fq, '\n')){
          std::cout<< "Error!!! cannot get fq seq from " << f5file << std::endl;
          status = 4;
       }else if (!std::getline(iss, subelemt, '\n')){
          std::cout<< "Error!!! cannot get fq+ from " << f5file << std::endl;
          status = 5;
       }else if (!std::getline(iss, f5fq.qual, '\n')){
          std::cout<< "Error!!! cannot get fq qual from " << f5file << std::endl;
          status = 6;
       }

       f5fq.name_line = f5fq.fq_name;
       f5fq.fq_name.clear();

       m_runOption_pointer->read_num.clear();
       size_t end_id = f5fq.name_line.find_first_of(" _", 1);
       f5fq.fq_name = f5fq.name_line.substr(1, end_id-1);
       if (end_id != std::string::npos){
          size_t read_pos = f5fq.name_line.find("read", end_id); 
          if (read_pos != std::string::npos ){
              size_t read_num_pos = f5fq.name_line.find_first_of("0123456789", read_pos);
              if (read_num_pos != std::string::npos ){
                  size_t read_num_pos_end = f5fq.name_line.find_first_not_of("0123456789", read_num_pos);
                  if (read_num_pos_end != std::string::npos ){
                     m_runOption_pointer->read_num = "Read_" + f5fq.name_line.substr(read_num_pos, read_num_pos_end-read_num_pos);
                  }else {
                     m_runOption_pointer->read_num = "Read_" + f5fq.name_line.substr(read_num_pos);
                  }
              }
          }
       }

       f5fq.name_line.shrink_to_fit();
       f5fq.fq_name.shrink_to_fit();
       f5fq.fq.shrink_to_fit();
       f5fq.qual.shrink_to_fit();

       fq_ds.close();
   }catch( H5::FileIException error )
   {
      std::cout<< "Error!!! Canot open file "<< f5file << std::endl;
      error.printError();
      status = 1;
   }catch( H5::DataSetIException error )
   {
      std::cout << "Error!!! Cannot open fq dataset from " << f5file << std::endl;
      error.printError();
      status = 10;
   }

   return status;
}

int F5Reader::_group_difference_(const std::vector<double> & source_data, const int data_size){
   m_runOption_pointer->group_sum[0] = 0;
   int f_s_i;
   //std::cout<<m_runOption_pointer->group_sum[0];
   for (f_s_i=0; f_s_i < data_size; f_s_i++){
      m_runOption_pointer->group_sum[f_s_i+1] = source_data[f_s_i] + m_runOption_pointer->group_sum[f_s_i];
      //std::cout<<" " << m_runOption_pointer->group_sum[f_s_i+1];
   }
   //std::cout<< std::endl;
   for (f_s_i=m_runOption_pointer->group_size; f_s_i <= data_size- m_runOption_pointer->group_size; f_s_i++){
       m_runOption_pointer->group_dif[f_s_i-m_runOption_pointer->group_size] = fabs(m_runOption_pointer->group_sum[f_s_i]*2 - m_runOption_pointer->group_sum[f_s_i-m_runOption_pointer->group_size] - m_runOption_pointer->group_sum[f_s_i+m_runOption_pointer->group_size]);
      //std::cout<< f_s_i << " " << m_runOption_pointer->group_size << " " << m_runOption_pointer->group_sum[f_s_i] << " " << m_runOption_pointer->group_sum[f_s_i-m_runOption_pointer->group_size] << " " << m_runOption_pointer->group_sum[f_s_i+m_runOption_pointer->group_size] << " " << m_runOption_pointer->group_dif[f_s_i-m_runOption_pointer->group_size]<< std::endl;
   }

   /*std::cout<< " in _group_difference_ "<< m_runOption_pointer->group_size << " " << data_size <<std::endl;
   for (uint64_t mtesti=m_runOption_pointer->f5events[0].start; mtesti< m_runOption_pointer->f5events[4].start+m_runOption_pointer->f5events[4].length; mtesti++){
       std::cout<< " " << m_runOption_pointer->group_dif[mtesti];
   }
   std::cout<<std::endl;
   */
   return 0;
}

int F5Reader::group_difference(){
   if (m_runOption_pointer->group_size<2){
      std::cout<< " Error!!! group_size="<< m_runOption_pointer->group_size<< " is too small." << std::endl;
      status = 31;
      return status;
   }
   return _group_difference_(f5signal_list, signal_size);
}


int F5Reader::associate_signal_to_fq(){
   std::cout <<"signals-size=" <<signal_size << " events-size=" << events_size << "/fq=" << f5fq.fq.size() << " " << f5signal_list.size() << " " << m_runOption_pointer->group_size <<std::endl;
   double m_median;
   double median_abs_dev;
   _get_median_mad(f5signal_list, 0, f5signal_list.size(), m_median, median_abs_dev);
   for(int si=0; si<f5signal_list.size(); si++){
      double cur_nv = (f5signal_list[si] - m_median)/median_abs_dev;
      if (cur_nv > 5) { cur_nv = 5;}
      else if (cur_nv<-5) { cur_nv = -5; }
      f5signal_list[si] = cur_nv;
   }

   /*F5Event cur_f5etest = m_runOption_pointer->f5events[0];
   cur_f5etest.model_state[2] = 'N';
   std::cout<<F5Event_toString(cur_f5etest)<<std::endl;
   std::cout<<F5Event_toString(m_runOption_pointer->f5events[0])<<std::endl;
   return 0;*/
   
   group_difference();

   //std::vector<RankPos> largest_dif_pos = get_top_N_extreme(m_runOption_pointer->group_dif, m_runOption_pointer->f5events[0].start, m_runOption_pointer->f5events[4].start+m_runOption_pointer->f5events[4].length, 2, m_runOption_pointer->group_size);
   // for(int rpi=0; rpi<largest_dif_pos.size(); rpi++){
   //    std::cout<< largest_dif_pos[rpi].value << " / " << largest_dif_pos[rpi].pos<< std::endl;
   // }
   // return 0;

   if (status>0){
      std::cout<<"Error happened. Cannot further process the association of signals to fq" << std::endl;
      return status;
   }

   uint64_t last_Unprocessed_signal_ind = 0, last_Unprocessed_event_ind=0, f5_e_i=0, f5_fq_i=1, num_splitPoints=0;
   uint64_t dif_pos_ind;
   std::vector<uint64_t> less_split_point;
   for (f5_e_i=0; f5_e_i<events_size; f5_e_i++){
       if (m_runOption_pointer->f5events[f5_e_i].move>0){
          f5_fq_i += m_runOption_pointer->f5events[f5_e_i].move;
          num_splitPoints = m_runOption_pointer->f5events[f5_e_i].move;

          if (m_runOption_pointer->f5events[f5_e_i].model_state[2] == f5fq.fq[f5_fq_i]){;}
          else{
              std::cout<< f5_e_i << " " << m_runOption_pointer->f5events[f5_e_i].model_state <<"|Ev_fq!|" <<f5fq.fq[f5_fq_i] << " " << f5_fq_i<<std::endl;
          }

          if (f5_e_i==0){
             last_Unprocessed_signal_ind = m_runOption_pointer->f5events[f5_e_i].start;
             last_Unprocessed_event_ind = f5_e_i;
             num_splitPoints = 0;
          }else{
             //std::cout<<"Res_" << f5_e_i << ": "<< last_Unprocessed_event_ind << "/" << last_Unprocessed_signal_ind << " > "<<F5Event_toString(m_runOption_pointer->f5events[f5_e_i], f5signal_list) << std::endl;
             std::vector<RankPos> largest_dif_pos = get_top_N_extreme(m_runOption_pointer->group_dif, last_Unprocessed_signal_ind, m_runOption_pointer->f5events[f5_e_i].start+m_runOption_pointer->f5events[f5_e_i].length, num_splitPoints, m_runOption_pointer->group_size);
             /*std::cout<<"signal_sp: < ";
             std::cout<<m_runOption_pointer->f5events[last_Unprocessed_event_ind].start << "==" << last_Unprocessed_signal_ind;
             for (int rpi=0; rpi<largest_dif_pos.size(); rpi++){
                 std::cout<<", "<<largest_dif_pos[rpi].pos;
             }
             std::cout<<" >";
             for (uint64_t espi=last_Unprocessed_signal_ind; espi< m_runOption_pointer->f5events[f5_e_i].start+m_runOption_pointer->f5events[f5_e_i].length; espi++){
                 for (std::vector<RankPos>::iterator rp_it=largest_dif_pos.begin(); rp_it!=largest_dif_pos.end(); rp_it++){
                     if (rp_it->pos == espi){
                        std::cout<< " |" ;
                     }
                 }
                 std::cout<< " " << f5signal_list[espi];
             }
             std::cout<<std::endl;*/
             if (largest_dif_pos.size()<num_splitPoints){
                std::cout << "   " <<f5_e_i << " " << last_Unprocessed_signal_ind <<"/" <<last_Unprocessed_event_ind << " " << m_runOption_pointer->f5events[f5_e_i].start << "+" << m_runOption_pointer->f5events[f5_e_i].length << " sp"<< largest_dif_pos.size() << "/" << num_splitPoints << "="; 
                for (dif_pos_ind=0; dif_pos_ind <largest_dif_pos.size(); dif_pos_ind++){
                    std::cout<< " " << largest_dif_pos[dif_pos_ind].pos << "/" << largest_dif_pos[dif_pos_ind].value;
                }
                std::cout<<std::endl; 
             }

             for (dif_pos_ind=0; dif_pos_ind < num_splitPoints; dif_pos_ind++){
                F5Event cur_f5e;
                if (dif_pos_ind<largest_dif_pos.size()){ cur_f5e = m_runOption_pointer->f5events[last_Unprocessed_event_ind]; }
                else { cur_f5e = m_runOption_pointer->f5events[f5_e_i]; }
                cur_f5e.start = last_Unprocessed_signal_ind;
                if (dif_pos_ind<largest_dif_pos.size()){ cur_f5e.length = largest_dif_pos[dif_pos_ind].pos - cur_f5e.start; }
                else { cur_f5e.length = m_runOption_pointer->f5events[f5_e_i].start+m_runOption_pointer->f5events[f5_e_i].length - cur_f5e.start; }

                cur_f5e.mean = st_mean(cur_f5e.start, cur_f5e.length, f5signal_list);
                cur_f5e.stdv = st_std(cur_f5e.start, cur_f5e.length, cur_f5e.mean, f5signal_list);
                cur_f5e.move = 1;
                if (dif_pos_ind>0){
                   cur_f5e.model_state[0] = '-'; cur_f5e.model_state[1] = '-'; cur_f5e.model_state[3] = '-'; cur_f5e.model_state[4] = '-';
                   cur_f5e.model_state[2] = f5fq.fq[f5_fq_i - (num_splitPoints - dif_pos_ind)];
                }
                f5event_list.push_back(cur_f5e);

                if (dif_pos_ind<largest_dif_pos.size()){
                   last_Unprocessed_signal_ind = largest_dif_pos[dif_pos_ind].pos;
                   last_Unprocessed_event_ind = f5_e_i;
                }else{
                   less_split_point.push_back(f5event_list.size());
                }
                /*if (largest_dif_pos.size()<num_splitPoints){
                   std::cout<<"\t" << last_Unprocessed_signal_ind <<"/" <<last_Unprocessed_event_ind << " " << F5Event_toString_basic(f5event_list[f5event_list.size()-1]) <<std::endl;
                }*/
             }   
             num_splitPoints = 0;
          }
       }
   } 
   if (last_Unprocessed_event_ind <events_size){
      F5Event cur_f5e = m_runOption_pointer->f5events[last_Unprocessed_event_ind];
      cur_f5e.start = last_Unprocessed_signal_ind;
      cur_f5e.length = m_runOption_pointer->f5events[events_size-1].start+m_runOption_pointer->f5events[events_size-1].length - cur_f5e.start;
      cur_f5e.mean = st_mean(cur_f5e.start, cur_f5e.length, f5signal_list);
      cur_f5e.stdv = st_std(cur_f5e.start, cur_f5e.length, cur_f5e.mean, f5signal_list);
      cur_f5e.move = 1;
      f5event_list.push_back(cur_f5e);
   }

   for(int lspi=0; lspi<less_split_point.size(); lspi++){
      if (f5event_list[less_split_point[lspi]-1].start!=f5event_list[less_split_point[lspi]].start){
          std::cout << "Warning!!! Incorrect start points for less_split_point=" << lspi << " " << F5Event_toString_basic(f5event_list[less_split_point[lspi]-1]) << " ||| " << F5Event_toString_basic(f5event_list[less_split_point[lspi]]) << " ||| " << F5Event_toString_basic(f5event_list[less_split_point[lspi]+1]) << std::endl;
      }
 
      f5event_list[less_split_point[lspi]-1].length = f5event_list[less_split_point[lspi]].length/2;
      f5event_list[less_split_point[lspi]].start = f5event_list[less_split_point[lspi]-1].start + f5event_list[less_split_point[lspi]-1].length;
      f5event_list[less_split_point[lspi]].length = f5event_list[less_split_point[lspi]].length - f5event_list[less_split_point[lspi]-1].length;     

      f5event_list[less_split_point[lspi]-1].mean = st_mean(f5event_list[less_split_point[lspi]-1].start, f5event_list[less_split_point[lspi]-1].length, f5signal_list);
      f5event_list[less_split_point[lspi]-1].stdv = st_std(f5event_list[less_split_point[lspi]-1].start, f5event_list[less_split_point[lspi]-1].length, f5event_list[less_split_point[lspi]-1].mean, f5signal_list);

      f5event_list[less_split_point[lspi]].mean = st_mean(f5event_list[less_split_point[lspi]].start, f5event_list[less_split_point[lspi]].length, f5signal_list);
      f5event_list[less_split_point[lspi]].stdv = st_std(f5event_list[less_split_point[lspi]].start, f5event_list[less_split_point[lspi]].length, f5event_list[less_split_point[lspi]].mean, f5signal_list); 
   }
   /*
   f5_fq_i = 1;
   int m_i;
   for (f5_e_i=0; f5_e_i<events_size; f5_e_i++){
       if (m_runOption_pointer->f5events[f5_e_i].move>0){
          f5_fq_i += m_runOption_pointer->f5events[f5_e_i].move;
          m_i = m_runOption_pointer->f5events[f5_e_i].move-1;
          while(m_i>0){
              std::cout<< f5_e_i << " -" << "=" << f5fq.fq[f5_fq_i-m_i] << " " << F5Event_toString_basic(m_runOption_pointer->f5events[f5_e_i-m_i]) << " >>> " << F5Event_toString_basic(f5event_list[f5_fq_i - 2 - m_i]) << std::endl;
              m_i--;
          }
          std::cout<< f5_e_i << " " << m_runOption_pointer->f5events[f5_e_i].model_state[2] << "=" << f5fq.fq[f5_fq_i] << " " << F5Event_toString_basic(m_runOption_pointer->f5events[f5_e_i]) << " >>> " << F5Event_toString_basic(f5event_list[f5_fq_i - 2]) << std::endl;

       }
   }
   uint64_t opsize = 50;
   if (f5event_list.size()<opsize*2){
      opsize = f5event_list.size()/2;
   }
   std::cout<< "size: "<< f5event_list.size() << " " << f5fq.fq.size() << std::endl;
   std::cout<< "Bseq_ev=" ;
   uint64_t op_ind;
   for (op_ind=0; op_ind<opsize; op_ind++){
      std::cout<< f5event_list[op_ind].model_state[2];
   }
   std::cout<<std::endl;
   std::cout<< "Bseq_fq=" << f5fq.fq.substr(2,opsize) <<std::endl;
   std::cout<< "Eseq_ev=";
   for (op_ind=f5event_list.size()-opsize; op_ind<f5event_list.size(); op_ind++){
      std::cout<< f5event_list[op_ind].model_state[2];
   }
   std::cout<<std::endl;
   std::cout<< "Eseq_fq=" << f5fq.fq.substr(f5fq.fq.size()-opsize-2, opsize) <<std::endl;
   */

   //////////////////////////////////////////////////////////
   // 
   //
   //
   //
   //
   //
   //
   ///////////////////////////////////////////////

   return status;

   ///////////////////////////////////////////////////////////////

   //uint64_t f5_e_i = 0, f5_fq_i = 0;
   uint64_t f5_fq_gap_i;
   for (f5_e_i=0; f5_e_i<events_size; f5_e_i++){
      float mean = st_mean(m_runOption_pointer->f5events[f5_e_i].start, m_runOption_pointer->f5events[f5_e_i].length, f5signal_list);
      if (fabs(mean-m_runOption_pointer->f5events[f5_e_i].mean)>0.001){
         std::cout<< "test f5e="<< f5_e_i <<":" << F5Event_toString(m_runOption_pointer->f5events[f5_e_i]) << " cal:" << mean << "\t\t cal stdv=" << st_std(m_runOption_pointer->f5events[f5_e_i].start, m_runOption_pointer->f5events[f5_e_i].length, mean, f5signal_list) << std::endl;      
      }
      m_runOption_pointer->f5events[f5_e_i].mean = mean;
      m_runOption_pointer->f5events[f5_e_i].stdv = st_std(m_runOption_pointer->f5events[f5_e_i].start, m_runOption_pointer->f5events[f5_e_i].length, mean, f5signal_list);
   }

   std::cout<< std::endl;
   for (f5_e_i=0; f5_e_i<100;  f5_e_i++){
      std::cout<< " check signals="<< f5_e_i <<":" << F5Event_toString(m_runOption_pointer->f5events[f5_e_i]);
      for (int opi=m_runOption_pointer->f5events[f5_e_i].start; opi<m_runOption_pointer->f5events[f5_e_i].start+m_runOption_pointer->f5events[f5_e_i].length; opi++){
          std::cout<< " " << f5signal_list[opi];
      }
      std::cout<< std::endl;
   }
   std::cout<< std::endl;

   std::string e_seq;
   std::vector<uint64_t> ev_ind_list;
   std::map<uint64_t, uint64_t > ev_ind_to_fq_ind;
   std::map<uint64_t, uint64_t> ev_pre;
   uint64_t last_move1_ind;
   std::map<uint64_t, uint64_t> ev_suf;
   bool is_small_stdv_move0 = false;
   uint64_t ind_small_stdv_move0;
   for (f5_e_i=0; f5_e_i<events_size; f5_e_i++){
      if (m_runOption_pointer->f5events[f5_e_i].move==0){;}
      else if (m_runOption_pointer->f5events[f5_e_i].move>0){
         last_move1_ind = f5_e_i;
         if( is_small_stdv_move0 ){
            is_small_stdv_move0 = false;
            ev_suf[ind_small_stdv_move0] = f5_e_i;
            std::cout<<" A>>"<<F5Event_toString(m_runOption_pointer->f5events[f5_e_i])<< std::endl ;
         }
         if (m_runOption_pointer->f5events[f5_e_i].move>2){
            std::cout<< "Warning!!! more than 2 move:" << F5Event_toString(m_runOption_pointer->f5events[f5_e_i])<< std::endl;
         }
         
         f5_fq_i += m_runOption_pointer->f5events[f5_e_i].move;

         if (m_runOption_pointer->f5events[f5_e_i].move>1) {
            for (f5_fq_gap_i=0; f5_fq_gap_i < m_runOption_pointer->f5events[f5_e_i].move-1; f5_fq_gap_i++){
               e_seq.push_back('-');
            }
         }
         if (e_seq.size()==0) {
            e_seq.push_back(m_runOption_pointer->f5events[f5_e_i].model_state[0]);
            e_seq.push_back(m_runOption_pointer->f5events[f5_e_i].model_state[1]);
         }
         
         e_seq.push_back(m_runOption_pointer->f5events[f5_e_i].model_state[2]);
      }

      if (m_runOption_pointer->f5events[f5_e_i].stdv<0.5){
          ev_ind_list.push_back(f5_e_i);
          ev_ind_to_fq_ind[f5_e_i] = f5_fq_i-1;
          is_small_stdv_move0 = false;
          if (m_runOption_pointer->f5events[f5_e_i].move==0){
             ev_pre[f5_e_i] = last_move1_ind;
             std::cout<< " B>>"<<F5Event_toString(m_runOption_pointer->f5events[last_move1_ind])<< std::endl;
             is_small_stdv_move0 = true;
             ind_small_stdv_move0 = f5_e_i;
          }
          std::cout<< f5_e_i << " " << f5_fq_i-1 << " " << ev_ind_list.size() <<" " << ev_ind_to_fq_ind.size() << " mean=" << m_runOption_pointer->f5events[f5_e_i].mean << " stdv=" << m_runOption_pointer->f5events[f5_e_i].stdv << " move=" <<  m_runOption_pointer->f5events[f5_e_i].move << std::endl;
      }
   }
   f5_e_i = events_size-1;
   if (m_runOption_pointer->f5events[f5_e_i].move==0){
      e_seq.push_back(m_runOption_pointer->f5events[f5_e_i].model_state[3]);
      e_seq.push_back(m_runOption_pointer->f5events[f5_e_i].model_state[4]);
   }

   std::cout<< "size: "<< e_seq.size() << " " << f5fq.fq.size() << std::endl;
   std::cout<< "Bseq=" << e_seq.substr(0,50) << std::endl;
   std::cout<< "Bseq=" << f5fq.fq.substr(0,50) <<std::endl;
   std::cout<< "Eseq=" << e_seq.substr(e_seq.size()-50) << std::endl;
   std::cout<< "Eseq=" << f5fq.fq.substr(f5fq.fq.size()-50) <<std::endl;

   //double m_median;
   //double median_abs_dev;
   _get_median_mad(f5signal_list, 0, f5signal_list.size(), m_median, median_abs_dev);
   for(int si=0; si<f5signal_list.size(); si++){
      double cur_nv = (f5signal_list[si] - m_median)/median_abs_dev;
      if (cur_nv > 5) { cur_nv = 5;}
      else if (cur_nv<-5) { cur_nv = -5; }
      f5signal_list[si] = cur_nv;
   }
   return status;
}

int F5Reader::_readF5Events_(std::string f5e_str){
   try{
       H5::Exception::dontPrint();
       //f5e_str = "/Analyses/Basecall_1D_001/BaseCalled_template/Events";

       H5::DataSet event_ds = m_f5_file_ptr->openDataSet(f5e_str);
       H5::DataSpace event_sp = event_ds.getSpace(); 
       events_size = event_sp.getSimpleExtentNpoints(); 
       if (events_size >= MAX_F5EVENT_SIZE){
           std::cout << "Error!!! Too much events " << events_size << " more than" << MAX_F5EVENT_SIZE << f5file << std::endl;
          status = 13;
       }else{

          H5::CompType et_type(sizeof(F5Event));
          et_type.insertMember(F5Path::get_mean_str(), HOFFSET(F5Event, mean), H5::PredType::IEEE_F32LE);
          et_type.insertMember(F5Path::get_stdv_str(), HOFFSET(F5Event, stdv), H5::PredType::IEEE_F32LE);
          et_type.insertMember(F5Path::get_start_str(), HOFFSET(F5Event, start), H5::PredType::STD_U64LE);
          et_type.insertMember(F5Path::get_length_str(), HOFFSET(F5Event, length), H5::PredType::STD_U64LE);
          et_type.insertMember(F5Path::get_model_state_str(), HOFFSET(F5Event, model_state), H5::StrType(H5::PredType::C_S1, KMER_SIZE));
          et_type.insertMember(F5Path::get_move_str(), HOFFSET(F5Event, move), H5::PredType::STD_I32LE);

          event_ds.read(m_runOption_pointer->f5events, et_type);
          
          /*for (int si = 0; si < 50; si++){
             if (m_runOption_pointer->f5events[si].move==0) { continue; }
             float mean = st_mean(m_runOption_pointer->f5events[si].start, m_runOption_pointer->f5events[si].length, f5signal_list);
             std::cout<< "test f5e="<< si <<":" << F5Event_toString(m_runOption_pointer->f5events[si]) << " cal:" << mean << "\t\t cal stdv=" << st_std(m_runOption_pointer->f5events[si].start, m_runOption_pointer->f5events[si].length, mean, f5signal_list) << std::endl;
          }*/
          /*float mean = st_mean(m_runOption_pointer->f5events[0].start, m_runOption_pointer->f5events[0].length, f5signal_list);
          std::cout<< "test f5e:" << F5Event_toString(m_runOption_pointer->f5events[0]) << " cal:" << mean << " " << st_std(m_runOption_pointer->f5events[0].start, m_runOption_pointer->f5events[0].length, mean, f5signal_list) << std::endl;  
          mean = st_mean(m_runOption_pointer->f5events[1].start, m_runOption_pointer->f5events[1].length, f5signal_list);
          std::cout<< "test f5e:" << F5Event_toString(m_runOption_pointer->f5events[1]) << " cal:" << mean << " " << st_std(m_runOption_pointer->f5events[1].start, m_runOption_pointer->f5events[1].length, mean, f5signal_list) << std::endl;
          mean = st_mean(m_runOption_pointer->f5events[2].start, m_runOption_pointer->f5events[2].length, f5signal_list);
          std::cout<< "test f5e:" << F5Event_toString(m_runOption_pointer->f5events[2]) << " cal:" << mean << " " << st_std(m_runOption_pointer->f5events[2].start, m_runOption_pointer->f5events[2].length, mean, f5signal_list) << std::endl; */
       }
       event_ds.close();
   }catch( H5::FileIException error )
   {
      std::cout<< "Error!!! Canot open file "<< f5file << std::endl;
      error.printError();
      status = 1;
   }catch( H5::DataSetIException error )
   {
      std::cout << "Error!!! Cannot open event dataset from " << f5file << std::endl;
      error.printError();
      status = 11;
   }catch( H5::DataSpaceIException error )
   {
      std::cout << "Error!!! Cannot open event dataspace from " << f5file << std::endl;
      error.printError();
      status = 14;
   }catch( H5::DataTypeIException error )
   {
      std::cout << "Error!!! Cannot open event datatype from " << f5file << std::endl;
      error.printError();
      status = 15;
   }

   return status;
}

int F5Reader::_readF5Signals_(std::string f5s_str){
   try{
      //std::cout<< " iin f5s_str" << f5s_str << std::endl;
      H5::Exception::dontPrint();
      H5::DataSet signal_ds = m_f5_file_ptr->openDataSet(f5s_str); //path.get_signal_path()+"/"+m_runOption_pointer->read_id+"/Signal");

      H5::DataType mdatatype= signal_ds.getDataType();
      signal_ds.read(f5fq.name_line, mdatatype);

      ////  H5T_STD_I16LE 
      // size would be wrong if not H5T_STD_I16LE 
      //signal_size = signal_ds.getInMemDataSize()/2;
      //signal_size = signal_ds.getStorageSize();
      //
      //This size is correct.
      signal_size = signal_ds.getSpace().getSimpleExtentNpoints();
      if (signal_size >= MAX_F5SIGNAL_SIZE){
          std::cout << "Error!!! Too much signals " << signal_size << " more than" << MAX_F5SIGNAL_SIZE << f5file << std::endl;
          status = 12;
      }else{
          signal_ds.read(m_runOption_pointer->f5signals, mdatatype);
      }

      signal_ds.close();
   }catch( H5::FileIException error )
   {
      std::cout<< "Error!!! Canot open file "<< f5file << std::endl;
      error.printError();
      status = 1;
   }catch( H5::DataSetIException error )
   {
      std::cout << "Error!!! Cannot open signal dataset from " << f5file << std::endl;
      error.printError();
      status = 11;
   }

   return status;
}

int F5Reader::normalize_signal(){
   if (status!=0){ 
      return status;
   }

   if (f5signal_list.size()>0){
       std::cout << "Warning!!! signals have been normalized. Will redo it. " << f5file << std::endl;
       f5signal_list.clear();
   }

   // for test;
   /*double min_signals = 1000, max_signals = 0;
   std::map<int, int64_t> signal_dict;
   for (int sdi=0; sdi<200; sdi++){
       signal_dict[sdi] = 0;
   }*/
   // test edn

   for (int si=0; si<signal_size; si++){
      f5signal_list.push_back((m_runOption_pointer->f5signals[si] + f5channel.offset)*f5channel.range/f5channel.digitisation);

      // fortest
      /*if(f5signal_list[si]<min_signals) { min_signals = f5signal_list[si]; }
      if(f5signal_list[si]>max_signals) { max_signals = f5signal_list[si]; }
      signal_dict[int(f5signal_list[si])] += 1;
      if (f5signal_list[si]<0){
         std::cout<< " negative value = " << si << std::endl;
      }
      if (si>0 && f5signal_list[si-1]<0){
         if(si>1){
            f5signal_list[si-1] = (f5signal_list[si-2] + f5signal_list[si])/2;
         }else{
            f5signal_list[si-1] = f5signal_list[si];
         }
      }*/ 
      // test end;
   } 

   // for test
   /*std::cout<< "min_signals=" << min_signals << " max_signals=" << max_signals << std::endl;
   for (int sdi=0; sdi<200; sdi++){
      if (signal_dict[sdi]!=0){
          std::cout << sdi << " " << signal_dict[sdi] <<std::endl;
      }
   }*/
   // test end

   /*
   double m_median;
   double median_abs_dev;
   //_get_median_mad(f5signal_list, 5000, f5signal_list.size()-5000, m_median, median_abs_dev);
   //std::cout << " median="<<m_median << " mad=" << median_abs_dev <<std::endl;
   _get_median_mad(f5signal_list, 0, f5signal_list.size(), m_median, median_abs_dev);
   //std::cout << " median="<<m_median << " mad=" << median_abs_dev <<std::endl;
   for(int si=0; si<f5signal_list.size(); si++){
      double cur_nv = (f5signal_list[si] - m_median)/median_abs_dev;
      if (cur_nv > 5) { cur_nv = 5;}
      else if (cur_nv<-5) { cur_nv = -5; }
      f5signal_list[si] = cur_nv;
   }*/
}
int F5Reader::_get_median_mad(const std::vector<double> & m_data, uint64_t start_ind, uint64_t end_ind, double & m_median, double & median_abs_dev){
   std::vector<double> find_mean;
   for (int si=start_ind; si<end_ind; si++){
      find_mean.push_back(m_data[si]);
   }
   // for test
   /*for (int testi=0; testi<20; testi++){
      std::cout<< int(f5signal_list[testi])<< " " ;
   }
   std::cout<< std::endl;*/
   // test end
   std::sort (find_mean.begin(), find_mean.end());
   // for test
   /*for (int testi=0; testi<20; testi++){
      std::cout<< int(f5signal_list[testi])<< " " ;
   }
   std::cout<< std::endl;*/
   // test end;

   if (find_mean.size()%2==1){
      int md_ind = find_mean.size()/2;
      m_median = (find_mean[md_ind] + find_mean[md_ind+1])/2;
   }else{
      m_median = find_mean[find_mean.size()/2];
   }   

   for(int si=0; si<find_mean.size(); si++){
      find_mean[si] = fabs(find_mean[si] - m_median);
   }
   std::sort (find_mean.begin(), find_mean.end());
   // for test
   /*for (int testi=0; testi<20; testi++){
      std::cout<< int(f5signal_list[testi])<< " " ;
   }
   std::cout<< std::endl;*/
   // test end;
   if (find_mean.size()%2==1){
      int md_ind = find_mean.size()/2;
      median_abs_dev = (find_mean[md_ind] + find_mean[md_ind+1])/2;
   }else{
      median_abs_dev = find_mean[find_mean.size()/2];
   }

   return status;
}


