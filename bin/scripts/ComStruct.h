#ifndef COMSTRUCT_H_
#define COMSTRUCT_H_

#include <string>
#include <map>
#include <vector>


/////////////////////////////////
#define MAX_F5EVENT_SIZE   3000000
#define MAX_F5SIGNAL_SIZE 15000000

#define KMER_SIZE 6

//////////////////////////////////
#define MODULE_NUM 5
#define CHAR_SIZE 1024

#define WhiteSpace "\t\n\v\f\r "

//////////////////////////////////

typedef struct F5Event{
   float mean;
   float stdv;
   uint64_t start;
   uint64_t length;
   char model_state[KMER_SIZE];
   uint32_t move;
} F5Event;


///////////////////////////////////////
typedef struct RankPos{
   double value;
   uint64_t pos;
} RankPos;

///////////////////////////////////
typedef std::basic_string<unsigned char> uc8string;

typedef struct Map1Base{
  char qry_base;
  char ref_base;
} Map1Base;
// 0-based, like bam and bed format
typedef struct Map1BasePos{
  uint64_t qry_pos;
  uint64_t ref_pos;
  uint64_t map_type;
} Map1BasePos;

////////////////////////////

typedef struct RunOption{
  uint64_t group_size; 
  std::string read_num;
  std::string read_id;
  F5Event f5events[MAX_F5EVENT_SIZE];
  int16_t f5signals [MAX_F5SIGNAL_SIZE];
  double group_dif[MAX_F5SIGNAL_SIZE+1];
  double group_sum[MAX_F5SIGNAL_SIZE+1];
  std::map<std::string, std::vector<std::string> > errorInfo;

  std::map<std::string, std::string> readToF5;
  std::map<std::string, bool> readOfInterest;
} RunOption;

// 0-based, like bam and bed format; end-not-included 
typedef struct GenomicRegion{
   int64_t start_pos;
   int64_t end_pos;
   char chrn[CHAR_SIZE];
} GenomicRegion;

// 0-based, like bam and bed format; end-not-included
typedef struct RepeatRegion{
   int64_t start_pos;
   int64_t end_pos;
   char chrn[CHAR_SIZE];
   int len_repeat_unit;
} RepeatRegion;

//////////////////////////////


#endif

