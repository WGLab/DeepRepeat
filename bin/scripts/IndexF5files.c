#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>

#include <map>
#include <iostream>
#include <string>

#include "ComFunction.h"
#include "Fast5Index.h"

int main (int argc, char * argv[])
{
   if (argc < 5){
      std::cout<< "Usage: " << argv[0] << " Save-index-path basecalled-path unique-file-id seq-sum-file multif5" << std::endl;
      return 1;
   }

   std::string base_index_path(argv[1]);
   std::string basecalled_path(argv[2]);
   std::string uniq_id(argv[3]);
   std::string seq_sum(argv[4]);
   bool multifast5 = atoi(argv[5])>0;

   std::cout<< "Options to be used:"<<std::endl;
   std::cout<< "\tbase_index_path = "<<base_index_path<<std::endl;
   std::cout<< "\tbasecalled_path = "<<basecalled_path<<std::endl;
   std::cout<< "\tuniq_id = "<<uniq_id<<std::endl;
   std::cout<< "\tseq_sum = "<<seq_sum<<std::endl;
   std::cout<< "\tmultifast5 = "<<multifast5<<std::endl;

   build_f5index(base_index_path, basecalled_path, uniq_id, seq_sum, multifast5);

   return EXIT_SUCCESS;
}
