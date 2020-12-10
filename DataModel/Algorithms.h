// -*- mode:c++; tab-width:4; -*-
/* vim:set noexpandtab tabstop=4 wrap */
#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <iostream>
#include <vector>
#include <map>
//#include <regex>  // in g++4.8 (in use on sukap/SK container) std::regex is broken
#include <boost/regex.hpp> // since std::regex doesn't work
#include <boost/regex/pattern_except.hpp>
#include <sstream>
#include <fstream>   // for ofstream

#include "basic_array.h"
#include "OutputRedirector.h"   // for CStdoutRedirector

class TVector3;
class TLorentzVector;

int ReadListFromFile(std::string filename, std::vector<std::string> &lines, char commentchar='#', bool trim_whitespace=true);
std::string GetStdoutFromCommand(std::string cmd, int bufsize=500);
void SetRootColourPlotStyle();
double MomentumToEnergy(basic_array<float[3]>& mom, int pdg);
double MomentumToEnergy(TVector3& mom, int pdg);
double Mag2(basic_array<float>& mom);
double Mag(basic_array<float>& mom);
bool CheckPath(std::string path, std::string& type);
std::string ToLower(std::string astring);

// a header/trailer for trying to capture error messages for calls to functions
// that do not return any error code, but do print to stderr on error.
std::streambuf* start_stderr_capture();
std::string end_stderr_capture(std::streambuf* previous_buff);

namespace algorithms{
	
} // end namespace algorithms

// helper function: to_string with a precision
// particularly useful for printing doubles and floats in the Log function
template <typename T>
std::string toString(const T a_value, const int n = 2){
	std::ostringstream out;
	out.precision(n);
	out << std::fixed << a_value;
	return out.str();
}

template <typename T>
std::string toString(T* a_ptr){
	std::stringstream out;
	out<<a_ptr;
	return out.str();
}

template <typename T>
bool solveQuadratic(const T &a, const T &b, const T &c, T &x0, T &x1){
	T discr = b*b - 4*a*c;
	if(discr < 0){ return false; }
	else if(discr==0){ x0=-0.5*(b/a); x1=x0; }
	else {
		T q = (b>0) ? -0.5*(b+sqrt(discr)) : -0.5*(b-sqrt(discr));
		x0 = q/a;
		x1 = c/q;
	}
	if (x0>x1) std::swap(x0, x1);
	
	return true;
}

template <typename T>
void printVals(T &container, int messagelevel, int verbosity, std::string preamble="", std::string postamble=""){
	// TODO: this is ok for vectors and arrays, but won't work for maps
	// TODO: generalise the Log function to be variadic, and accept containers, which are automatically
	// unwrapped in a suitable way. Could we even handle tuples? Rarely used though
	// if messagelevel==0 (error), print to cerr, otherwise print to cout
	std::ofstream outputbuf;
	if(messagelevel==0){
		outputbuf.copyfmt(std::cerr);
		outputbuf.clear(std::cerr.rdstate());
		outputbuf.basic_ios<char>::rdbuf(std::cerr.rdbuf());
	} else {
		outputbuf.copyfmt(std::cout);
		outputbuf.clear(std::cout.rdstate());
		outputbuf.basic_ios<char>::rdbuf(std::cout.rdbuf());
	}
	outputbuf<<preamble<<" {";
	for(auto it=container.begin(); it!=container.end(); ++it){ outputbuf<<(*it); }
	outputbuf<<"} "<<postamble<<std::endl;
}

//// https://stackoverflow.com/a/16111317
//// version that takes no arguments
//template<typename FunctionType>
//std::string getOutputFromFunctionCall(FunctionType&& FunctionToCall){
//	std::cout<<"version with no args"<<std::endl;
//	// set up capturing of stdout, capturing their current nominal outputs
//	std::stringbuf tempbuf;
//	std::streambuf* nominal_stdout_bufferp = std::cout.rdbuf(&tempbuf);
//	std::streambuf* nominal_stderr_bufferp = std::cerr.rdbuf(&tempbuf);
//	std::streambuf* nominal_stdlog_bufferp = std::clog.rdbuf(&tempbuf);
//	std::cout<<"redirected stream test"<<std::endl;
//	// invoke the requested function
//	std::forward<FunctionType>(FunctionToCall)();
//	// restore the streams, capturing the intermediate buffer
//	std::cout.rdbuf(nominal_stdout_bufferp);
//	std::cerr.rdbuf(nominal_stderr_bufferp);
//	std::clog.rdbuf(nominal_stdlog_bufferp);
//	return tempbuf.str();
//}

//// version that accepts arbitrary arguments?
//template<typename FunctionType, typename... Rest>
//std::string getOutputFromFunctionCall(FunctionType&& FunctionToCall, Rest... rest){
//	std::cout<<"version with args"<<std::endl;
//	// set up capturing of stdout, capturing their current nominal outputs
//	std::stringbuf tempbuf;
//	std::streambuf* nominal_stdout_bufferp = std::cout.rdbuf(&tempbuf);
//	std::streambuf* nominal_stderr_bufferp = std::cerr.rdbuf(&tempbuf);
//	std::streambuf* nominal_stdlog_bufferp = std::clog.rdbuf(&tempbuf);
//	std::cout<<"redirected stream test"<<std::endl;
//	// invoke the requested function
//	std::forward<FunctionType>(FunctionToCall)(rest...);
//	// restore the streams, capturing the intermediate buffer
//	std::cout.rdbuf(nominal_stdout_bufferp);
//	std::cerr.rdbuf(nominal_stderr_bufferp);
//	std::clog.rdbuf(nominal_stdlog_bufferp);
//	return tempbuf.str();
//}

// the above versions do not capture ROOT output.
// that method of redirecting std::cout doesn't capture all output to stdout
// using https://github.com/zpasternack/Redirector instead
// version that takes no arguments
template<typename FunctionType>
std::string getOutputFromFunctionCall(FunctionType&& FunctionToCall){
	// set up capturing of stdout, capturing their current nominal outputs
	CStdoutRedirector theRedirector;
	theRedirector.StartRedirecting();
	// invoke the requested function
	std::forward<FunctionType>(FunctionToCall)();
	// restore the streams, capturing the intermediate buffer
	std::string captured_out = theRedirector.GetOutput();
	theRedirector.ClearOutput();
	theRedirector.StopRedirecting();
	return captured_out;
}

// version that accepts arbitrary arguments?
template<typename FunctionType, typename... Rest>
std::string getOutputFromFunctionCall(FunctionType&& FunctionToCall, Rest... rest){
	// set up capturing of stdout, capturing their current nominal outputs
	CStdoutRedirector theRedirector;
	theRedirector.StartRedirecting();
	// invoke the requested function
	std::forward<FunctionType>(FunctionToCall)(rest...);
	// restore the streams, capturing the intermediate buffer
	std::string captured_out = theRedirector.GetOutput();
	theRedirector.ClearOutput();
	theRedirector.StopRedirecting();
	return captured_out;
}

#endif

