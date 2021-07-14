/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef TreeReader_H
#define TreeReader_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "MTreeReader.h"
#include "SkrootHeaders.h" // MCInfo, Header etc.
#include "Constants.h"

/**
* \class TreeReader
*
* FIXME
*
* $Author: M.O'Flaherty $
* $Date: 2019/05/28 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class TreeReader: public Tool {
	
	public:
	TreeReader();         ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute();   ///< Execute function used to perform Tool purpose.
	bool Finalise();  ///< Finalise funciton used to clean up resources.
	
	private:
	// functions
	// =========
	int ReadEntry(long entry_number);
	int LoadConfig(std::string configfile);
	int GenerateNewLUN();
	void CloseLUN();
	
	// tool variables
	// ==============
	std::string toolName;
	std::string inputFile="";
	std::string FileListName="InputFileList";
	std::string selectionsFile="";
	std::string cutName="";
	std::string treeName;
	std::string readerName;
	int maxEntries=-1;
	int firstEntry=0;
	int entrynum=0;
	int readEntries=0; // count how many TTree entries we've actually processed
	SKROOTMODE skrootMode = SKROOTMODE::READ; // default to read
	int skreadMode=0;  // 0=skread only, 1=skrawread only, 2=both
	int skreadUser=0;  // 0=auto, 1=skread only, 2=skrawread only, 3=both
	int LUN=10; // This is assumed 10 by some SKROOT algorithms so only override if you know what you're doing!
	std::string skroot_options="31";  // 31 = read HEADER (required).
	int skroot_badopt=23;             // 23 = LOWE default (mask bad chs, dead chs, noisy ID chs and OD chs)
	int skroot_badch_ref_run=0;       // reference run for bad channel list for e.g. MC.
	int sk_geometry=4;                // TODO increment the default to 6.
	std::string outputFile="";        // for when using SKROOT copy mode
	int skip_ped_evts = 1; // automatically skip to next ttree entry if skread returns 1 (pedestal/status entry)
	
	std::vector<std::string> list_of_files;
	
	// verbosity levels: if 'verbosity' < this level, the message type will be logged.
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	std::string logmessage="";
	int get_ok=0;
	
	// variables to read in
	// ====================
	MTreeReader myTreeReader;                  // the TTree reader
	MTreeSelection* myTreeSelections=nullptr;  // a set of entries one or more cuts
	std::vector<std::string> ActiveInputBranches;
	std::vector<std::string> ActiveOutputBranches;
	
};


#endif
