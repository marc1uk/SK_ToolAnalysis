/* vim:set noexpandtab tabstop=4 wrap */
#include "RootReadTest.h"

#include "TROOT.h"
#include "TInterpreter.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TParameter.h"

#include "type_name_as_string.h"

#include <iostream>
#include <vector>
#include <sstream>

RootReadTest::RootReadTest():Tool(){}


bool RootReadTest::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="") m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	m_variables.Get("inputFile",inputFile);
	m_variables.Get("treeName",treeName);
	m_variables.Get("testFileType",testFileType);
	
	myTreeReader.Load(inputFile, treeName);
	myTreeReader.SetVerbosity(1);
	entrynum=0;
	
	return true;
}

bool RootReadTest::Execute(){
	std::cout<<"ReadRootTest getting entry "<<entrynum<<std::endl;
	
	if(testFileType=="official_ntuple"){
		int get_ok = ReadEntryNtuple(entrynum);
		CheckEntryNtuple();
	}
	
	if(testFileType=="SKROOT"){
		int get_ok = ReadEntrySKROOT(entrynum);
		CheckEntrySKROOT();
	}
	
	entrynum++;
	if(entrynum==3){
		std::cout<<"setting StopLoop"<<std::endl;
		m_data->vars.Set("StopLoop",1);
	}
	
	return true;
}

// Official Ntuple
// ---------------

int RootReadTest::ReadEntryNtuple(long entry_number){
	int bytesread = myTreeReader.GetEntry(entrynum);
	if(bytesread<=0) return false;
	
	int success = 
	(myTreeReader.GetBranchValue("nscndprt", n_secondaries_2)) &&
	(myTreeReader.GetBranchValue("iprtscnd", secondary_PDG_code_2)) &&
	(myTreeReader.GetBranchValue("vtxscnd", secondary_start_vertex_2));
	
	return success;
}

int RootReadTest::CheckEntryNtuple(){
	std::cout<<"we had "<<n_secondaries_2<<" secondaries"<<std::endl;
	
	std::cout<<"we had "<<secondary_PDG_code_2.size()<<" secondary PDG codes: {";
	for(auto&& asecondary : secondary_PDG_code_2){
		std::cout<<asecondary<<", ";
	}
	std::cout<<"\b\b}"<<std::endl;
	
	std::cout<<"we had "<<secondary_start_vertex_2.size()<<" secondary start vertices: {";
	for(int i=0; i<secondary_start_vertex_2.size(); ++i){
		auto&& avertex = secondary_start_vertex_2.at(i);
		std::cout<<"[";
		for(auto&& aval : avertex){
			std::cout<<aval<<", ";
		}
		std::cout<<"\b\b], ";
	}
	std::cout<<"\b\b}"<<std::endl;
	return 1;
}

// SKROOT files
// -------------

int RootReadTest::ReadEntrySKROOT(long entry_number){
	int bytesread = myTreeReader.GetEntry(entrynum);
	if(bytesread<=0) return false;
	
	int success = 
	(myTreeReader.GetBranchValue("MC", mc_info)) &&
	(myTreeReader.GetBranchValue("HEADER", file_header));
	
	return success;
}

int RootReadTest::CheckEntrySKROOT(){
	std::cout<<"MCInfo is at "<<mc_info<<std::endl;
	std::cout<<"we had "<<mc_info->nvc<<" primaries"<<std::endl;
	
	for(int primary_i=0; primary_i<mc_info->nvc; ++primary_i){
		std::cout<<"primary "<<primary_i<<" had pdg code "<<mc_info->ipvc[primary_i]
				 <<" and initial momentum ("<<mc_info->pvc[primary_i][0]
				 <<", "<<mc_info->pvc[primary_i][1]<<", "<<mc_info->pvc[primary_i][2]<<"); ";
	}
	
	std::cout<<"Header is at "<<file_header<<", and indicates nevsk "
			 <<file_header->nevsk<<" and swtrig_id "<<file_header->swtrg_id<<std::endl;
	
	return 1;
}


bool RootReadTest::Finalise(){
	
	return true;
}
