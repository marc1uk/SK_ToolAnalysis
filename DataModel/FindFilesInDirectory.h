// -*- mode:c++; tab-width:4; -*-
/* vim:set noexpandtab tabstop=4 wrap */
#ifndef FindFilesInDirectory_H
#define FindFilesInDirectory_H

#include "Algorithms.h"

#include <iostream>
#include <map>
//#include <regex>  // in g++4.8 (in use on sukap/SK container) std::regex is broken
#include <boost/regex.hpp> // since std::regex doesn't work
#include <boost/regex/pattern_except.hpp>

int FindFilesInDirectory(std::string inputdir, std::string pattern, std::vector<std::string> &matches, bool case_sensitive=false, int max_subdir_depth=0, bool use_regex=false, std::vector<std::string>* filenames=nullptr, std::vector<std::vector<std::string>>* output_submatches=nullptr, bool verbose=false);

namespace algorithms{
//	static const std::map<std::regex_constants::error_type,std::string> regex_err_strings{
//	{std::regex_constants::error_collate, "The expression contained an invalid collating element name."},
//	{std::regex_constants::error_ctype, "The expression contained an invalid character class name."},
//	{std::regex_constants::error_escape, "The expression contained an invalid escaped character, or a trailing escape."},
//	{std::regex_constants::error_backref, "The expression contained an invalid back reference."},
//	{std::regex_constants::error_brack, "The expression contained mismatched brackets ([ and ])."},
//	{std::regex_constants::error_paren, "The expression contained mismatched parentheses (( and ))."},
//	{std::regex_constants::error_brace, "The expression contained mismatched braces ({ and })."},
//	{std::regex_constants::error_badbrace, "The expression contained an invalid range between braces ({ and })."},
//	{std::regex_constants::error_range, "The expression contained an invalid character range."},
//	{std::regex_constants::error_space, "There was insufficient memory to convert the expression into a finite state machine."},
//	{std::regex_constants::error_badrepeat, "The expression contained a repeat specifier (one of *?+{) that was not preceded by a valid regular expression."},
//	{std::regex_constants::error_complexity, "The complexity of an attempted match against a regular expression exceeded a pre-set level."},
//	{std::regex_constants::error_stack, "There was insufficient memory to determine whether the regular expression could match the specified character sequence."}};
	
	static const std::map<boost::regex_constants::error_type,std::string> bregex_err_strings{
	{boost::regex_constants::error_collate, "The expression contained an invalid collating element name."},
	{boost::regex_constants::error_ctype, "The expression contained an invalid character class name."},
	{boost::regex_constants::error_escape, "The expression contained an invalid escaped character, or a trailing escape."},
	{boost::regex_constants::error_backref, "The expression contained an invalid back reference."},
	{boost::regex_constants::error_brack, "The expression contained mismatched brackets ([ and ])."},
	{boost::regex_constants::error_paren, "The expression contained mismatched parentheses (( and ))."},
	{boost::regex_constants::error_brace, "The expression contained mismatched braces ({ and })."},
	{boost::regex_constants::error_badbrace, "The expression contained an invalid range between braces ({ and })."},
	{boost::regex_constants::error_range, "The expression contained an invalid character range."},
	{boost::regex_constants::error_space, "There was insufficient memory to convert the expression into a finite state machine."},
	{boost::regex_constants::error_badrepeat, "The expression contained a repeat specifier (one of *?+{) that was not preceded by a valid regular expression."},
	{boost::regex_constants::error_complexity, "The complexity of an attempted match against a regular expression exceeded a pre-set level."},
	{boost::regex_constants::error_stack, "There was insufficient memory to determine whether the regular expression could match the specified character sequence."},
	{boost::regex_constants::error_bad_pattern, "Invalid regex pattern."}};
	
} // namespace algorithms

#endif
