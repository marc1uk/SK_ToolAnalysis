/* vim:set noexpandtab tabstop=4 wrap */
#include "MyToolRoot.h"

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"

MyToolRoot::MyToolRoot():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

bool MyToolRoot::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",verbosity);            // how verbose to be
	m_variables.Get("inputFile",inputFile);            // a single specific input file
	m_variables.Get("outputFile",outputFile);          // output file to write
	m_variables.Get("maxEvents",maxEvents);            // user limit to number of events to process
	
	// open the input TFile and TTree
	// ------------------------------
	get_ok = myTreeReader.Load(inputFile, "eventtree");
	
	return true;
}


bool MyToolRoot::Execute(){
	
	Log(toolName+" getting entry "+toString(entrynum),v_debug,verbosity);
	
	// retrieve desired branches
	get_ok = GetBranches();
	
	// process the data
	// TODO
	
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
	
	return true;
}


bool MyToolRoot::Finalise(){
	
	return true;
}

int MyToolRoot::ReadEntry(long entry_number){
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

int MyToolRoot::GetBranches(){
	int success = (
//		(myTreeReader.GetBranchValue("filename",filename))                         &&
//		(myTreeReader.GetBranchValue("gamma_time",gamma_time))
		true
	);
	
	return success;
}
