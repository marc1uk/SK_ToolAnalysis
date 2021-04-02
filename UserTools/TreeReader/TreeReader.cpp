/* vim:set noexpandtab tabstop=4 wrap */
#include "TreeReader.h"
#include "TTree.h"

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"
#include "MTreeSelection.h"

TreeReader::TreeReader():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

bool TreeReader::Initialise(std::string configfile, DataModel &data){
	
	m_data= &data;
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
	// Get the Tool configuration variables
	// ------------------------------------
	LoadConfig(configfile);
	
	// open the input TFile and TTree
	// ------------------------------
	Log(toolName+" creating MTreeReader to read file "+inputFile+", tree "+treeName,v_debug,verbosity);
	get_ok = myTreeReader.Load(inputFile, treeName);
	if(not get_ok){
		Log(toolName+" failed to open reader on tree "+treeName+" in file "+inputFile,v_error,verbosity);
		return false;
	}
	
	// for efficiency of reading, only enable used branches
	Log(toolName+" activating branches",v_debug,verbosity);
	if(ActiveBranches.size() && 
		std::find(ActiveBranches.begin(), ActiveBranches.end(), "*")!=ActiveBranches.end()){
		// only disable unlisted branches if we have a non-empty list of active branches
		// and the key "*" was not specified.
		myTreeReader.OnlyEnableBranches(ActiveBranches);
	}
	
	if(firstEntry>0) entrynum = firstEntry;
	
	// put the reader into the DataModel
	Log(toolName+" registering tree reader "+readerName,v_debug,verbosity);
	m_data->Trees.emplace(readerName,&myTreeReader);
	
	// if we were given a selections file, only read entries that pass the specified cut
	if(selectionsFile!=""){
		// make the MTreeSelection to read the TEntryList file
		myTreeSelections = new MTreeSelection(selectionsFile);
		m_data->Selectors.emplace(readerName,myTreeSelections);
		
		if(cutName=="") cutName = myTreeSelections->GetTopCut();
		Log(toolName+" reading only entries passing cut "+cutName
			+" in selections file "+selectionsFile,v_debug,verbosity);
		
		// scan to the first entry passing our specified cut
		Log(toolName+" scanning to first passing entry",v_debug,verbosity);
		do {
			entrynum = myTreeSelections->GetNextEntry(cutName);
		} while(entrynum<firstEntry);
		Log(toolName+" reading from entry "+toString(entrynum),v_debug,verbosity);
	}
	
	
	return true;
}

bool TreeReader::Execute(){
	
	Log(toolName+" getting entry "+toString(entrynum),v_debug,verbosity);
	
	// load next entry
	// errors are already handled in ReadEntry
	get_ok = ReadEntry(entrynum);
	
	// check if we've hit the user-requested entry limit
	// or if this is the last entry, in which case set StopLoop
	if(myTreeSelections==nullptr){
		entrynum++;
	} else {
		entrynum = myTreeSelections->GetNextEntry(cutName);
	}
	if((maxEntries>0)&&(entrynum>=(maxEntries+firstEntry))){
		Log(toolName+" hit max events, setting StopLoop",v_message,verbosity);
		m_data->vars.Set("StopLoop",1);
		return true;
	}
	else if(myTreeReader.GetTree()->LoadTree(entrynum)<0){
		// use LoadTree to check if the next entry is valid without loading it
		Log(toolName+" reached end of TTree, setting StopLoop",v_message,verbosity);
		m_data->vars.Set("StopLoop",1);
		return true;
	}
	
	return true;
}


bool TreeReader::Finalise(){
	
	if(myTreeSelections) delete myTreeSelections;
	return true;
}

int TreeReader::ReadEntry(long entry_number){
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

int TreeReader::LoadConfig(std::string configfile){
	Log(toolName+" reading configuration",v_debug,verbosity);
	// read the config file
	std::ifstream fin (configfile.c_str());
	std::string Line;
	
	if(not fin.is_open()){
		Log(toolName+" failed to read configuration file "+configfile,v_error,verbosity);
		return -1;
	}
	
	bool settingBranchNames=false;
	std::vector<std::string> ActiveBranches;
	
	// scan over lines in the config file
	while (getline(fin, Line)){
		Log(toolName+" parsing config line \""+Line+"\"",v_debug,verbosity);
		// skip empty lines
		if (Line.empty()) continue;
		std::string LineCopy = Line; // make a copy so we can print it in case of parsing error
		// trim preceding whitespace
		Line.erase(0,Line.find_first_not_of(" \t\015"));
		// skip comment lines
		if(Line[0] == '#') continue;
		// trim line end comments (everything after and including a '#' character)
		if(Line.find('#')!=std::string::npos) Line.erase(Line.find_first_of('#'),std::string::npos);
		// trim trailing whitespace
		if(Line.find_last_not_of(" \t\n\015\014\013")!=std::string::npos)
			Line.erase(Line.find_last_not_of(" \t\n\015\014\013")+1,std::string::npos);
		
		// split apart the key and value
		std::string thekey   = Line.substr(0,Line.find_first_of(" \t\n\015\014\013"));
		std::string thevalue = Line.substr(Line.find_first_of(" \t\n\015\014\013")+1,std::string::npos);
		
		// first check if we're entering or leaving the list of branches to activate
		if (thekey=="StartBranchList"){
			settingBranchNames = true;
		}
		else if(thekey=="EndBranchList"){
			settingBranchNames = false;
		}
		else if(settingBranchNames){
			ActiveBranches.push_back(Line);
		}
		else if(thekey=="verbosity") verbosity = stoi(thevalue);
		else if(thekey=="inputFile") inputFile = thevalue;
		else if(thekey=="treeName") treeName = thevalue;
		else if(thekey=="readerName") readerName = thevalue;
		else if(thekey=="firstEntry") firstEntry = stoi(thevalue);
		else if(thekey=="maxEntries") maxEntries = stoi(thevalue);
		else if(thekey=="selectionsFile") selectionsFile = thevalue;
		else if(thekey=="cutName") cutName = thevalue;
		else {
			Log(toolName+" error parsing config file line: \""+LineCopy
				+"\" - unrecognised variable \""+thekey+"\"",v_error,verbosity);
		}
		
	}
	
	// done parsing, close config file
	fin.close();
	
	return 1;
}


