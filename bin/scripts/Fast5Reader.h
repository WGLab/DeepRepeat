#ifndef FAST5READER_H_
#define FAST5READER_H_

#include <stdint.h>

#include <string>
#include <vector>

#include "H5Cpp.h"
#include "ComFunction.h"
#include "ComStruct.h"

class F5Path{
protected:
   std::string fastq_path;
   std::string event_path;
   std::string channel_path;
   std::string signal_path;
   int status;
   int read_config(const char * path_config_file);
   
   static std::string digiti_str;
   static std::string offset_str;
   static std::string samplr_str;
   static std::string range_str;

   static std::string mean_str;
   static std::string stdv_str;
   static std::string start_str;
   static std::string length_str;
   static std::string model_state_str;
   static std::string move_str;
public:
   std::string get_fastq_path();
   std::string get_event_path();
   std::string get_channel_path();
   std::string get_signal_path();
   int get_status();

   void toPrint();

   //F5Path(const F5Path & m_f5p); 
   F5Path(const char * path_config_file);
   F5Path(){;}
   ~F5Path(){;}

   static std::string get_digiti_str(){ return digiti_str; }
   static std::string get_offset_str(){ return offset_str; }
   static std::string get_samplr_str(){ return samplr_str; }
   static std::string get_range_str(){ return range_str; }

   static std::string get_mean_str(){ return mean_str; }
   static std::string get_stdv_str(){ return stdv_str; }
   static std::string get_start_str(){ return start_str; }
   static std::string get_length_str(){ return length_str; }
   static std::string get_model_state_str(){ return model_state_str; }
   static std::string get_move_str(){ return move_str; }
};

class F5Fastq {
public:
   std::string name_line;
   std::string fq_name;
   std::string fq;
   std::string qual;
  
   int reset();
};

/*
typedef struct F5Event{
   float mean;
   float stdv;
   uint64_t start;
   uint64_t length;
   char model_state[KMER_SIZE];
   uint32_t move;
//   float weights;
//   float p_model_state;
//   char model_state[5];
//   float p_mp_state;
//   float p_A;
//   float p_C;
//   float p_G
//   float p_T; 
} F5Event;*/
std::string F5Event_toString(const F5Event & m_f5event);
std::string F5Event_toString_basic(const F5Event & m_f5event);
std::string F5Event_toString(const F5Event & m_f5event, const std::vector<double>& f5s);


class F5Channel{
   int status;
public:
   double digitisation; 
   double offset;
   double range;
   double sampling_rate;

   F5Channel();
   int get_status();
   void toPrint();
   int reset();
};

class F5Reader{
public:
   F5Channel f5channel;
   F5Fastq f5fq;
   std::vector<F5Event> f5event_list;
   std::vector<double> f5signal_list;
protected:
   RunOption *m_runOption_pointer;
   
   std::string f5file;

   H5::H5File *m_f5_file_ptr;
   int status;
   hsize_t signal_size;
   hsize_t events_size;

   int _readF5Channel_(std::string f5c_str);
   int _readF5Fastq_(std::string f5f_str);
   int _readF5Events_(std::string f5e_str);
   int _readF5Signals_(std::string f5s_str);

   int _init_(const char* m_f5file, RunOption *op_runOption_pointer);
   int _del_();
   int _get_median_mad(const std::vector<double> & m_data, uint64_t start_ind, uint64_t end_ind, double & m_median, double & median_abs_dev);
   //int _group_difference_(const std::vector<double> & source_data, const int data_size);
public:
   static F5Path f5path;
   virtual int readF5Channel(){;}
   virtual int readF5Fastq(){;}
   virtual int readF5Events(){;}
   virtual int readF5Signals(){;}

   virtual int readF5data();

   F5Reader(const char* m_f5file, RunOption *op_runOption_pointer);
   int resetFast5File(const char* m_f5file, RunOption *op_runOption_pointer=NULL);
   F5Reader(RunOption *op_runOption_pointer);
   virtual ~F5Reader();
   
   int normalize_signal();
   F5Channel get_F5Channel(){ return f5channel; }
   F5Fastq get_F5Fastq() { return f5fq; }
   hsize_t get_signal_size(){ return signal_size; }
   hsize_t get_events_size(){ return events_size; }
   int get_status(){ return status; }
   
   int associate_signal_to_fq();
   int group_difference();
   int _group_difference_(const std::vector<double> & source_data, const int data_size);
};

class F5ReaderMulti:public F5Reader{
public:
   int readF5Channel();
   int readF5Fastq();
   int readF5Events();
   int readF5Signals();

   F5ReaderMulti(const char* m_f5file, RunOption *op_runOption_pointer):F5Reader(m_f5file, op_runOption_pointer){;}
   F5ReaderMulti(RunOption *op_runOption_pointer):F5Reader(op_runOption_pointer){;}
   ~F5ReaderMulti(){;}
};

class F5Reader1:public F5Reader{
public:
   int readF5Channel();
   int readF5Fastq();
   int readF5Events();
   int readF5Signals();
   
   F5Reader1(const char* m_f5file, RunOption *op_runOption_pointer):F5Reader(m_f5file, op_runOption_pointer){;}
   F5Reader1(RunOption *op_runOption_pointer):F5Reader(op_runOption_pointer){;}
   ~F5Reader1(){;}
};


#endif

