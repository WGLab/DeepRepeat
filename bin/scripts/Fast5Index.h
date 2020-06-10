#ifndef FAST5INDEX_H_
#define FAST5INDEX_H_

#include <map>
#include <vector>
#include <string>

#include "ComStruct.h"

int read_f5index(const std::map<std::string, bool>& readOfInterest, std::map<std::string, std::string> & readToF5,  std::string fast5_index_file);
int build_f5index(const std::string base_index_path, const std::string basecalled_path, const std::string uniq_id, const std::string seq_sum, const bool multifast5);

#endif
