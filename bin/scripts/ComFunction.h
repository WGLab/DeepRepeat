#ifndef COMFUNCTION_H_
#define COMFUNCTION_H_

#include <stdio.h>

#include <map>
#include <vector>
#include <string>

#include <limits>
#include <iomanip>

#include "ComStruct.h"

#define get_array_size(m_a) sizeof(m_a)/sizeof(m_a[0])

#ifdef WINDOWS
    #include <direct.h>
    #define GetCWD _getcwd
#else
    #include <unistd.h>
    #define GetCWD getcwd
#endif


bool isExist(const char *);

int get_Files_From_Pattern(const std::string & file_pat, std::map<std::string, std::string> * file_dict);
int get_Files_From_Pattern(const std::string & file_pat, std::vector<std::string> * pat_files);

extern char m_cwd[];
int get_cwd();

double st_mean(const uint64_t start_pos, const uint64_t length, const std::vector<double>& dvlist);
double st_std(const uint64_t start_pos, const uint64_t length, const double mean, const std::vector<double>& dvlist);

bool compare_RankPos_v (const RankPos& rp1, const RankPos& rp2);
bool compare_RankPos_p (const RankPos& rp1, const RankPos& rp2); 

std::vector<RankPos> get_top_N_extreme(const double * data, const uint64_t m_start, uint64_t m_end, uint16_t topN=1, uint16_t sep_dist=4, int min_max=1);

bool m_cp_str(char * dest, const char * msource, int mlen);

std::vector<std::string> m_split_string(const std::string & m_str, std::string multi_delimiters=WhiteSpace, bool contain_delimiter=false);

#endif
