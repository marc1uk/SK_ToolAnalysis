// -*- mode:c++; tab-width:4; -*-
/* vim:set noexpandtab tabstop=4 wrap */
#include "Constants.h"

std::string G3_process_code_to_string(int process_code){
	if(constants::G3_process_code_to_string.count(process_code)){
		return constants::G3_process_code_to_string.at(process_code);
	} else {
		return "unknown";
	}
}

std::string numnu_code_to_string(int numnu_code){
	if(constants::numnu_code_to_string.count(numnu_code)){
		return constants::numnu_code_to_string.at(numnu_code);
	} else if(numnu_code>5){
		return "fs_other";
	} else {
		return "unknown";
	}
}

std::string neut_mode_to_string(int neut_code){
	if(constants::neut_mode_to_string.count(neut_code)){
		return constants::neut_mode_to_string.at(neut_code);
	} else {
		return "unknown";
	}
}


