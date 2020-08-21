/* vim:set noexpandtab tabstop=4 wrap */
#include "PlotNeutronCaptures.h"

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TVector3.h"
#include "TLorentzVector.h"

PlotNeutronCaptures::PlotNeutronCaptures():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

bool PlotNeutronCaptures::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="") m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",verbosity);            // how verbose to be
	m_variables.Get("inputFile",inputFile);            // a single specific input file
	m_variables.Get("outputFile",outputFile);          // output file to write
	m_variables.Get("drawPlots",DrawPlots);            // show root plots while working? otherwise just save
	m_variables.Get("maxEvents",maxEvents);            // user limit to number of events to process
	
	// open the input TFile and TTree
	// ------------------------------
	get_ok = myTreeReader.Load(inputFile, "eventtree");
	
	if(DrawPlots){
		// Only one TApplication may exist. Get it, or make it if there isn't one
		int myargc=0;
		intptr_t tapp_ptr=0;
		get_ok = m_data->CStore.Get("RootTApplication",tapp_ptr);
		if(not get_ok){
			if(verbosity>2) std::cout<<toolName<<": making global TApplication"<<std::endl;
			rootTApp = new TApplication("rootTApp",&myargc,0);
			tapp_ptr = reinterpret_cast<intptr_t>(rootTApp);
			m_data->CStore.Set("RootTApplication",tapp_ptr);
		} else {
			if(verbosity>2) std::cout<<toolName<<": Retrieving global TApplication"<<std::endl;
			rootTApp = reinterpret_cast<TApplication*>(tapp_ptr);
		}
		int tapplicationusers;
		get_ok = m_data->CStore.Get("RootTApplicationUsers",tapplicationusers);
		if(not get_ok) tapplicationusers=1;
		else tapplicationusers++;
		m_data->CStore.Set("RootTApplicationUsers",tapplicationusers);
	}
	
	return true;
}


bool PlotNeutronCaptures::Execute(){
	
	Log(toolName+" getting entry "+toString(entrynum),v_debug,verbosity);
	
	// ttree entry is already loaded so just retrieve the desired branches
	get_ok = GetBranches();
	
	// process the data
	get_ok = FillHistos();
	
	// move to next entry
	entrynum++;
	// check if we've hit the user-requested entry limit
	if((maxEvents>0)&&(entrynum==maxEvents)){
		Log(toolName+" hit max events, setting StopLoop",v_message,verbosity);
		m_data->vars.Set("StopLoop",1);
		return 1;
	}
	
	// pre-load the next ttree entry
	get_ok = ReadEntry(entrynum);
	if(get_ok==0){
		return 1; // end of file
	} else if (get_ok<0){
		return 0; // read error
	}
	
	return get_ok;
}

int PlotNeutronCaptures::FillHistos(){
	// TODO
	return 1;
}


bool PlotNeutronCaptures::Finalise(){
	
	// unregister ourselves as a user of the TApplication, and delete it if we're the last one
	if(DrawPlots){
		int tapplicationusers=0;
		get_ok = m_data->CStore.Get("RootTApplicationUsers",tapplicationusers);
		if(not get_ok || tapplicationusers==1){
			if(rootTApp){
				std::cout<<toolName<<": Deleting global TApplication"<<std::endl;
				delete rootTApp;
				rootTApp=nullptr;
			}
		} else if(tapplicationusers>1){
			m_data->CStore.Set("RootTApplicationUsers",tapplicationusers-1);
		}
	}
	
	return true;
}

int PlotNeutronCaptures::ReadEntry(long entry_number){
	// load next entry data from TTree
	int bytesread = myTreeReader.GetEntry(entry_number);
	
	// stop loop if we ran off the end of the tree
	if(bytesread==0){
		Log(toolName+" hit end of input file, stopping loop",v_message,verbosity);
		m_data->vars.Set("StopLoop",1);
	}
	// stop loop if we had an error of some kind
	else if(bytesread<0){
		 if(bytesread==-1) Log(toolName+" IO error loading next input entry!",v_error,verbosity);
		 if(bytesread==-2) Log(toolName+" AutoClear error loading next input entry!",v_error,verbosity);
		 if(bytesread <-2) Log(toolName+" Unknown error "+toString(bytesread)
		                       +" loading next input entry!",v_error,verbosity);
		 m_data->vars.Set("StopLoop",1);
	}
	
	return bytesread;
}

int PlotNeutronCaptures::GetBranches(){
	int success = (
//	(myTreeReader.GetBranchValue("filename",filename))                         &&
//	(myTreeReader.GetBranchValue("water_transparency",water_transparency))     &&
//	(myTreeReader.GetBranchValue("entry_number",entry_number))                 &&
//	(myTreeReader.GetBranchValue("subevent_num",subevent_number))              &&
//	(myTreeReader.GetBranchValue("primary_pdg",primary_pdg))                   &&
//	(myTreeReader.GetBranchValue("primary_energy",primary_energy))             &&
//	(myTreeReader.GetBranchValue("primary_start_pos",primary_start_pos))       &&
//	(myTreeReader.GetBranchValue("primary_end_pos",primary_end_pos))           &&
//	(myTreeReader.GetBranchValue("nuclide_pdg",nuclide_pdg))                   &&
//	(myTreeReader.GetBranchValue("nuclide_creation_pos",nuclide_creation_pos)) &&
//	(myTreeReader.GetBranchValue("nuclide_decay_pos",nuclide_decay_pos))       &&
	(myTreeReader.GetBranchValue("nuclide_daughter_pdg",nuclide_daughter_pdg)) &&
	(myTreeReader.GetBranchValue("neutron_start_pos",neutron_start_pos))       &&
	(myTreeReader.GetBranchValue("neutron_end_pos",neutron_end_pos))           &&
	(myTreeReader.GetBranchValue("neutron_start_energy",neutron_start_energy)) &&
	(myTreeReader.GetBranchValue("neutron_end_energy",neutron_end_energy))     &&
	(myTreeReader.GetBranchValue("neutron_end_process",neutron_end_process))   &&
//	(myTreeReader.GetBranchValue("neutron_n_daughters",neutron_ndaughters))    &&
	(myTreeReader.GetBranchValue("gamma_energy",gamma_energy))                 &&
	(myTreeReader.GetBranchValue("gamma_time",gamma_time))
	);
	
	return success;
}

