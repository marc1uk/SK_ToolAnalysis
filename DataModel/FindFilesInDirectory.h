// -*- mode:c++; tab-width:4; -*-
/* vim:set noexpandtab tabstop=4 wrap */
#ifndef FindFilesInDirectory_H
#define FindFilesInDirectory_H

#include "Algorithms.h"
#include "Constants.h"

#include <iostream>
#include <map>
//#include <regex>  // in g++4.8 (in use on sukap/SK container) std::regex is broken
#include <boost/regex.hpp> // since std::regex doesn't work
#include <boost/regex/pattern_except.hpp>

int FindFilesInDirectory(std::string inputdir, std::string pattern, std::vector<std::string> &matches, bool case_sensitive=false, int max_subdir_depth=0, bool use_regex=false, std::vector<std::string>* filenames=nullptr, std::vector<std::vector<std::string>>* output_submatches=nullptr, bool verbose=false);

#endif
