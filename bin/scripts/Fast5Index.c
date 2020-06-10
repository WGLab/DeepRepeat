#include "Fast5Index.h"
#include "ComFunction.h"
//#include "ComOption.h"

#include <map>
#include <iostream>

#include <fstream>
#include <sstream>

//#include <gzstream.h>


int read_f5index(const std::map<std::string, bool>& readOfInterest, std::map<std::string, std::string> & readToF5, std::string fast5_index_file){
   int run_status = 0;

   std::map<std::string, std::string>::iterator rTf5_it;
   std::string line;
   std::ifstream infile(fast5_index_file);
   if(!infile){
      fprintf(stderr, "Error!!! Cannot open file %s.\n", fast5_index_file.c_str());
      infile.close();
      return 1;
      //infile.close();
   }   

   std::map<std::string, bool>::iterator roi_it;

   std::cout<< "Read fast5_index_file=" << fast5_index_file << std::endl;
   while (std::getline(infile, line)){
      std::istringstream iss(line);
      std::string fast5file;
      std::string read_id;
      //if (!(iss >> fast5file >> read_id)) { break; }
      if (!(iss >> read_id >> fast5file)) { break; }
      
      if (readOfInterest.size()>0 && readOfInterest.find(read_id)==readOfInterest.end()){
         continue;
      }

      rTf5_it = readToF5.find(read_id);
      if (rTf5_it == readToF5.end()){
          readToF5[read_id] = fast5file;
      }else{
          fprintf(stdout, "Duplicate reads: %s(%s) %s(%s).\n", read_id.c_str(), fast5file.c_str(), rTf5_it->first.c_str(), rTf5_it->second.c_str());
      }
   }
   infile.close();
   std::cout<< "\tThere are " << readToF5.size() << " reads";
   if (readOfInterest.size()>0){
       std::cout << " for " << readOfInterest.size() << " reads of interest."  << std::endl;
   }else{
       std::cout << std::endl;
   }
   return run_status;
}
/*
 * base_index_path : where to save index path 
 * basecalled_path : basecalled workspace/pass path
 * uniq_id: unique file name id;
 * seq_sum : sequencing_summary.txt
 * multifast5 : is multifast5 format
 *
 */
int build_f5index(const std::string base_index_path, const std::string basecalled_path, const std::string uniq_id, const std::string seq_sum, const bool multifast5){
   int run_status = 0;

   // get all sequencing_summary.
   std::vector<std::string> seqsum_list;
   size_t star_pos = seq_sum.find_last_of("*");
   if (star_pos != std::string::npos){
      get_Files_From_Pattern(seq_sum, &seqsum_list);
   }else{
      seqsum_list.push_back(seq_sum);
   }
   std::vector<std::string>::iterator seqs_it;

   // check base_index_path has the same partent folder.
   if (seqsum_list.size()>1 && base_index_path.size()<1){
      fprintf(stderr, "Error!!! No base index path is provided when multiple sequencying_summary.txt are provided.\n");
      return 1;
   }else if (seqsum_list.size()==0){
      fprintf(stderr, "Error!!! Cannot find sequencying_summary.txt using %s\n", seq_sum.c_str());
      return 2;
   }

   // get save path for the index
   size_t pref_pos;
   pref_pos = seqsum_list[0].find_last_of("/");
   // get file to save index;
   std::string m_savePath(base_index_path.size()>0?base_index_path:(pref_pos!=std::string::npos?seqsum_list[0].substr(0,pref_pos):""));
   if (m_savePath.compare("./")==0){
      m_savePath = "";
   }
   size_t pref_path_size = m_savePath.size();
   if (pref_path_size==0 && seqsum_list[0].at(0)=='/'){
      fprintf(stderr, "Error!!! Absolute path of '-t' is used and base index path is not a part of '-t' files (%s).\n", seqsum_list[0].c_str());
      return 3;
   }else if (pref_path_size>0 && seqsum_list[0].substr(0,pref_path_size).compare(m_savePath)!=0){
      fprintf(stderr, "Error!!! Base index path (%s) is not a part of '-t' files (%s).\n", m_savePath.c_str(), seqsum_list[0].c_str());
      return 4;
   }

   if (m_savePath.size()>0 && m_savePath.at(m_savePath.size()-1)!='/'){
      m_savePath.append("/");
   }
   m_savePath.append(uniq_id);
   //m_savePath.append(".f5index.gz");
   m_savePath.append(".f5index");
   std::ofstream m_output_to_file(m_savePath.c_str(), std::ofstream::out);
   //ogzstream m_output_to_file(m_savePath.c_str());
   if (!(m_output_to_file)){
      std::cout<< " Cannot open " << m_savePath << " for write" << std::endl;
      m_output_to_file.close();
      return 5;
   }

   // output the inputs and output files.
   std::cout<<"\nThe input sequencing summary are:"<<std::endl;
   for(int sqs_i=0; sqs_i<seqsum_list.size(); sqs_i++){
       std::cout<< "\t" << seqsum_list[sqs_i]<<std::endl;
   }
   std::cout<<"The basecall path is <"<< basecalled_path << ">" <<std::endl;
   std::cout<<"The output index file is <" << m_savePath << ">" << std::endl;

   for (seqs_it=seqsum_list.begin(); seqs_it!=seqsum_list.end(); seqs_it++){
      // get base path for the sequencying_summary.txt
      std::string fast5_path(*seqs_it);
      fast5_path.erase(fast5_path.find_last_of("/")+1);

      // add further subdirectory path for fast5, and check whether the folder exist.
      fast5_path.append(basecalled_path);
      if (basecalled_path.at(basecalled_path.size()-1)!='/'){
         fast5_path.append("/");
      } 
      //if (!isExist(fast5_path.c_str())){
      //   fprintf(stderr, "Error!!! Folder (%s) does not exist.\n", fast5_path.c_str());
      //   return 6;
      //}

      // get the pattern of fast5 files, and get all fast5 files.
      // By default, check fast5 files directly under the folder if multifast5 is true
      // otherwise, check  fast5 files directly under the subdirectory of the folder
      std::map<std::string, std::string> f5fTopath;
      std::map<std::string, std::string>::iterator f5f_it;
      std::string f5pat(fast5_path);
      if (!multifast5){
         f5pat.append("*/*.fast5");
      }else{
         f5pat.append("*.fast5");
      }
      std::cout<< "   <" << f5pat << "> has ";
      get_Files_From_Pattern(f5pat.c_str(), &f5fTopath);
      std::cout<< f5fTopath.size() << " fast5 files" << std::endl;

      std::string line;
      std::ifstream infile(seqs_it->c_str());

      std::string fast5file;
      std::string read_id;
      std::string run_id;
      int channel;
      float start_time, duration;
      int num_events;
      std::string passes_filtering;
      std::getline(infile, line);
      while (std::getline(infile, line)){
         std::istringstream iss(line);
         //std::cout << ">>>" << line.substr(0, 50) << " --- ";
         if (!(iss >> fast5file >> read_id >> run_id >> channel >> start_time >> duration >> num_events >> passes_filtering )) { break; }
         //std::cout<< passes_filtering << "=" << (passes_filtering.compare("True")!=0) <<std::endl;

         if (passes_filtering.compare("True")!=0){ continue; }
         f5f_it = f5fTopath.find(fast5file);
         if (f5f_it==f5fTopath.end()){
            fprintf(stderr, "Warning!!! Cannot find %s\n", fast5file.c_str());
         }else{
            m_output_to_file << read_id << " " << f5f_it->second.substr(pref_path_size) << std::endl;
         }
      }
      infile.close();
   }
   
   m_output_to_file.close();
   std::cout<<std::endl;
   
   std::cout<< "Base index is saved in <" << m_savePath << ">." << std::endl;
   std::cout<< "Please keep this file together with its indexed fast5 folders." << std::endl;
   std::cout<< "You might need to provide this file for feature extraction and repeat detection." << std::endl;
   
   return run_status;
}



