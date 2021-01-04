// -*- mode:c++; tab-width:4; -*-
/* vim:set noexpandtab tabstop=4 wrap */
#include "Constants.h"

#include <iostream>

#include "TVector3.h"
#include "TLorentzVector.h"

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

std::string PdgToString(int code){
	if(constants::pdg_to_string.count(code)!=0){
		return constants::pdg_to_string.at(code);
	} else {
		return std::to_string(code);
	}
}

int StringToPdg(std::string name){
	if(constants::string_to_pdg.count(name)!=0){
		return constants::string_to_pdg.at(name);
	} else {
		return -1;
	}
}

std::string G3ParticleCodeToString(int code){
	if(constants::g3_particle_code_to_string.count(code)){
		return constants::g3_particle_code_to_string.at(code);
	} else {
		return std::to_string(code);
	}
}

int StringToG3ParticleCode(std::string name){
	if(constants::string_to_g3_particle_code.count(name)){
		return constants::string_to_g3_particle_code.at(name);
	} else {
		return -1;
	}
}

int G3ParticleCodeToPdg(int code){
	if(constants::g3_particle_code_to_pdg.count(code)){
		return constants::g3_particle_code_to_pdg.at(code);
	} else {
		return -1;
	}
}

int PdgToG3ParticleCode(int code){
	if(constants::pdg_to_g3_particle_code.count(code)){
		return constants::pdg_to_g3_particle_code.at(code);
	} else {
		return -1;
	}
}

void PrintVector(TVector3& avec, bool newline){
	std::string nl = (newline) ? "\n" : "";
	std::cout<<"("<<avec.X()<<", "<<avec.Y()<<", "<<avec.Z()<<")"<<nl;
}

void PrintVector(TLorentzVector& avec, bool newline){
	std::string nl = (newline) ? "\n" : "";
	std::cout<<"("<<avec.X()<<", "<<avec.Y()<<", "<<avec.Z()<<avec.T()<<")"<<nl;
}

double PdgToMass(int code){
	auto particle = constants::particleDb->GetParticle(code);
	if(particle==nullptr) return -1;
	return particle->Mass()*1000.;      // converted to MeV
}




