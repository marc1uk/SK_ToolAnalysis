// -*- mode:c++; tab-width:4; -*-
/* vim:set noexpandtab tabstop=4 wrap */
#include "Algorithms.h"
#include "Constants.h"
//#include <libgen.h>  // dirname and basename
#include <sys/stat.h>  // dirname and basename
#include <sys/types.h> // for stat() test to see if file or folder
#include <unistd.h>
//#include <memory>
//#include <exception>
//#include <cstring>  // strncpy
#include <fstream>
#include <sstream>
#include <cctype> // ::tolower

#include "TStyle.h"
#include "TColor.h"
#include "TVector3.h"
#include "TLorentzVector.h"

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

// XXX this doesn't always work, e.g. things like gObjectTable->Print() print to stdout
// but do not get captured. TODO replace with CStdoutRedirector as per getOutputFromFunctionCall
// in Algorithms.h
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

// as above but for stdout rather than stderr
std::streambuf* start_stdout_capture(){
	return std::cout.rdbuf(new std::stringbuf);
}

std::string end_stdout_capture(std::streambuf* previous_buff){
	std::stringbuf* capture_buffer = static_cast<std::stringbuf*>(std::cout.rdbuf(previous_buff));
	std::cout<<capture_buffer->str();
	std::string capture = capture_buffer->str();
	delete capture_buffer;
	return capture;
}

void SetRootColourPlotStyle(){
	const int NRGBs = 5;
	const int NCont = 255;
	
	double stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
	double red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
	double green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
	double blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
	gStyle->SetNumberContours(NCont);
}

double MomentumToEnergy(basic_array<float>& mom, int pdg){
	double mass = PdgToMass(pdg);
	if(mass<0) return -1;
	double momsq = pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2);
	return sqrt(momsq+pow(mass,2));
}

double MomentumToEnergy(TVector3& mom, int pdg){
	double mass = PdgToMass(pdg);
	if(mass<0) return -1;
	double momsq = pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2);
	return sqrt(momsq+pow(mass,2));
}

double Mag2(basic_array<float>& mom){
	return pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2);
}

double Mag(basic_array<float>& mom){
	return sqrt(Mag(mom));
}

bool CheckPath(std::string path, std::string& type){
	struct stat s;
	if(stat(path.c_str(),&s)==0){
		if(s.st_mode & S_IFDIR){        // mask to extract if it's a directory?? how does this work?
			type="d";  //it's a directory
			return true;
		} else if(s.st_mode & S_IFREG){ // mask to check if it's a file??
			type="f"; //it's a file
			return true;
		} else {
			// exists, but neither file nor directory?
			type="???";
			return false;
			//assert(false&&"Check input path: stat says it's neither file nor directory..?");
		}
	} else {
		// does not exist - could be a pattern, e.g. "/path/to/rootfiles_*.root"
		type="none";
		return false;
	}
	return false;
}

std::string ToLower(std::string astring){
	std::transform(astring.begin(), astring.end(), astring.begin(), ::tolower);  // why is the :: needed?
	return astring;
}
