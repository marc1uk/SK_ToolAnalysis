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

#include "basic_array.h"

class TVector3;
class TLorentzVector;

int ReadListFromFile(std::string filename, std::vector<std::string> &lines, char commentchar='#', bool trim_whitespace=true);
std::string GetStdoutFromCommand(std::string cmd, int bufsize=500);
void SetRootColourPlotStyle();
double MomentumToEnergy(basic_array<float[3]>& mom, int pdg);
double MomentumToEnergy(TVector3& mom, int pdg);
double Mag2(basic_array<float>& mom);
double Mag(basic_array<float>& mom);

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

#endif

