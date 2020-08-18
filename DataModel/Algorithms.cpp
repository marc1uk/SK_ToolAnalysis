// -*- mode:c++; tab-width:4; -*-
/* vim:set noexpandtab tabstop=4 wrap */
#include "Algorithms.h"
//#include <libgen.h>  // dirname and basename
//#include <sys/stat.h>  // dirname and basename
//#include <sys/types.h> // for stat() test to see if file or folder
//#include <unistd.h>
//#include <memory>
//#include <exception>
//#include <cstring>  // strncpy
#include <fstream>
#include <sstream>

int ReadListFromFile(std::string filename, std::vector<std::string> &lines, char commentchar, bool trim_whitespace){
	// read each new line into a std::vector<string> and return
	std::ifstream fin (filename.c_str());
	// return if not found or can't be opened
	if(not (fin.is_open() && fin.good())){
		return -1;
	}
	std::string Line;
	// loop over file lines
	while (getline(fin, Line)){
		if (Line[0] == commentchar){
			continue;
		} else{
			// check for trailing comments
			if (Line.find(commentchar) != std::string::npos){
				Line.erase(Line.begin()+Line.find(commentchar), Line.end());
			}
			// trim trailing whitespace
			if(trim_whitespace){
				Line.erase(Line.find_last_not_of(" \t\n\015\014\013")+1);
			}
			// add to vector
			lines.push_back(Line);
		}
	}
	fin.close();
	// return number of lines added
	return lines.size();
}

std::string GetStdoutFromCommand(std::string cmd, int bufsize){
	/*
	  credit: Jeremy Morgan, source:
	  https://www.jeremymorgan.com/tutorials/c-programming/how-to-capture-the-output-of-a-linux-command-in-c/
	*/
	std::string data;
	FILE * stream;
	char buffer[bufsize];
	cmd.append(" 2>&1");
	
	stream = popen(cmd.c_str(), "r");
	if(stream){
		while(!feof(stream)){
			if (fgets(buffer, bufsize, stream) != NULL) data.append(buffer);
		}
		pclose(stream);
	}
	return data;
}

// start capturing output from stderr. call this before any c++ functions that print to std::cerr
// in order to capture those printouts. XXX Don't forget to end_stderr_capture afterwards!!! XXX
std::streambuf* start_stderr_capture(){
	// redirect stderr to a custom stringbuf for capturing return from subsequent function calls
	// =========================================================================================
	// declare a new stringstream to catch cerr and point stderr at it
	// capture the return, which notes where it was previously
	return std::cerr.rdbuf(new std::stringbuf);
}

// stop capturing stderr. call this after your c++ function calls whose output you want to capture.
std::string end_stderr_capture(std::streambuf* previous_buff){
	// stop capturing stderr, retrieve output, and return normal behaviour
	// ===================================================================
	// restore stderr to its previous state, and retreieve our own streambuf
	std::stringbuf* capture_buffer = static_cast<std::stringbuf*>(std::cerr.rdbuf(previous_buff));
	// for backwards compatibility, print whatever was returned via the normal channel
	std::cerr<<capture_buffer->str();
	// get the captured contents
	std::string capture = capture_buffer->str();
	// free the temporary buffer
	delete capture_buffer;
	return capture;
}


