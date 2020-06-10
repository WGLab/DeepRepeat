#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>

#include <map>

#include "ComFunction.h"
#include "ComOption.h"
#include "Fast5Index.h"
#include "BamReader.h"

int main (int argc, char * argv[])
{
   bool has_error = false;
   if (get_cwd()!=0){
      has_error = true;
   }

   /*std::vector<std::string> seqsum_list;
   get_Files_From_Pattern("*.txt", &seqsum_list);
   for (std::vector<std::string>::iterator ss=seqsum_list.begin(); ss!=seqsum_list.end(); ss++){
      std::cout<<*ss<<std::endl;
   }
   std::string mtstr("/a");
   std::cout<< (mtstr.at(0)=='/') <<std::endl;*/

   init_ComOption();
   set_function();
   
   if (argc<2){
      print_help(argv[0], NULL);
      return 1;
   }

   has_error = set_options(argc, argv);
   //std::cout<< " for check in main: " << m_comOption.seq_sum << " " << m_comOption.basecalled_path << " " << m_comOption.base_index_path << ""  << std::endl;

   if (has_error){
      print_help(argv[0], argv[1]);
      return 1;
   }

   /*#define F5INDEX 0
   #define GETFEAT 1
   #define TRAIN 2
   #define PREDICT 3
   #define DETECT 4 */
   std::cout<< "Start to run " << argv[1] << std::endl;
   switch (m_comOption.submod_id){
       case F5INDEX:
           build_f5index();
           break;
       case GETFEAT:
           
           break;
       case TRAIN:

           break;
       case PREDICT:

           break;
       case DETECT:

           break;
       default:
           fprintf(stderr, "The function (%s) is supported.\n", argv[1]);
   }

   return EXIT_SUCCESS;
}
