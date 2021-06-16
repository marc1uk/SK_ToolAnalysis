/* vim:set noexpandtab tabstop=4 wrap */
#include "TreeReader.h"
#include "TTree.h"
#include <set>
#include <bitset>

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"
#include "MTreeSelection.h"
#include "fortran_routines.h"

TreeReader::TreeReader():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

const std::vector<std::string> default_branches{
	"HEADER",
	"TQREAL",
	"TQAREAL",
	"LOWE",
	"ATMPD",
	"UPMU",
	"MU",
	"SLE",
	"SWTRGLIST",
	"MC",
	"SECONDARY",
	"TQLIST",
	"ODTQLIST",
	"HWTRGLIST",
	"PEDESTALS",
	"EVENTHEADER",
	"EVENTTRAILER",
	"SOFTWARETRG",
	"QBEESTATUS",
	"DBSTATUS",
	"SPACERS",
	"PREVT0",
	"MISMATCHEDHITS",
	"GPSLIST",
	"T2KGPSLIST",
	"IDODXTLK"
};

bool TreeReader::Initialise(std::string configfile, DataModel &data){
	
	m_data= &data;
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
	// Get the Tool configuration variables
	// ------------------------------------
	LoadConfig(configfile);
	
	// safety check that we were given an input file
	if(inputFile=="" && FileListName==""){
		// unless we are working in SKROOT write mode...
		if(skroot==0 || skrootMode!=SKROOTMODE::WRITE){
			Log(toolName+" error! no InputFile or FileListName given!",v_error,verbosity);
			m_data->vars.Set("StopLoop",1);
			return false;
		}
	} else if(inputFile!=""){
		// single file takes precedence
		list_of_files.emplace_back(inputFile);
	} else {
		get_ok = m_data->CStore.Get(FileListName, list_of_files);
		if(!get_ok){
			Log(toolName+" error! Could not find file list "+FileListName+" in CStore!"
				+" Ensure LoadFileList tool is run before this tool!",v_error,verbosity);
			m_data->vars.Set("StopLoop",1);
			return false;
		}
	}
	
	// safety check that the requested name to associate to this reader is free
	get_ok = m_data->Trees.count(readerName);
	if(get_ok){
		Log(toolName+" error! TreeReader tool used to open file with name "+readerName
			+" but this name is already taken! Each name must be unique!",v_error,verbosity);
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	
	// safety check we have an output file too, if working in SKROOT write or copy mode
	if(skroot==1 && skrootMode!=SKROOTMODE::READ && outputFile==""){
		logmessage = toolName+" error! SKROOT mode is ";
		logmessage += ((skrootMode==SKROOTMODE::WRITE) ? "write" : "copy");
		logmessage += " but no outputFile specified!";
		Log(logmessage,v_error,verbosity);
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	
	// warning check: see if we're given an input when we're in WRITE mode
	if(skroot==1 && skrootMode==SKROOTMODE::WRITE && (inputFile!="" || FileListName!="")){
		Log(toolName+" warning! InputFile or FileListName given, but mode is skroot::write! "
					+"Inputs will be ignored: use another reader instance to read input from another file",
					v_error,verbosity);
	}
	// warning check: see if we've been given an output when we're not in SKROOT::COPY or WRITE mode
	if(outputFile!="" && (skroot!=1 || skrootMode==SKROOTMODE::READ)){
		logmessage  = toolName+" warning! outputFile given, but SKROOT mode is ";
		logmessage += ((skroot==0) ? "not enabled. " : "READ only. ");
		logmessage += "outputFile will be ignored!";
		Log(logmessage,v_warning,verbosity);
	}
	
	// open the input TFile and TTree
	// ------------------------------
	// a lot of SK algorithms retrieve data via skroot_get_* calls behind the scenes,
	// which means if we're processing skroot files we probably need a corresponding TreeManager.
	if(skroot){
		// skroot_open_* invokes the singleton SuperManager class to create a TreeManager
		// to be associated with a given SKROOT file. This TreeManager is just the usual
		// TTree wrapper created via ROOT's MakeClass. 
		// The SuperManager keeps a map of LUNs (unique IDs) for each TreeManager.
		// Unfortunately there is no way to know if a LUN is already in use or not.
		// Calling GetManager will return a nullptr if it doesn't, but then skroot_open_*
		// will fail if you try to subsequently open a new file with this LUN. Doh!
		// TODO fix the SuperManager.
		// For now we'll just keep our own list of LUNs.
		GenerateNewLUN();
		
		// There are 3 modes to the TreeManager:
		// skroot_open_read_ calls the TreeManager constructor with mode = 2;
		// Some branches may be skipped to optimize reading by using ZeroInBranch.
		
		// skroot_open_write_ calls the TreeManager constructor with mode =1;
		// It creates an output TTree with the full set of SKROOT branches.
		// You should call skroot_set_* or TreeManager::Set* methods to populate branches, then
		// call skroot_fill_tree (or TreeManager::fill_tree) to write a new output tree entry.
		// Note: ZeroOutBranch doesn't work in this mode, all output branches will be created.
		// Note: This mode is commented as "zbs 2 root" (or sometimes "root 2 zbs"),
		// but neither skroot functions nor the TreeManager provide any means of handling zbs files.
		// To do that, see e.g. $SKOFL_ROOT/examples/lowe/zbs2skroot.F
		
		// skroot_open_ calls the TreeManager constructor with mode = 0;
		// This allows input file reading (as per mode 2), but also creates an output file
		// with an SKROOT tree set up for copying events from input to output.
		// (i.e. both trees' branch addresses point to the same object).
		// Input read can be optimized via ZeroInBranch - these will not be read in or copied across.
		// Further branches may be omitted from the copy to output via ZeroOutBranch.
		// Output entries are written on each TreeManager::fill_tree() call,
		// so a subset of entries may be copied across.
		
		// create the treemanager, and in write mode, the output file
		switch(skrootMode){
			case SKROOTMODE::READ:  skroot_open_read_(&LUN); break;
			case SKROOTMODE::WRITE: skroot_open_write_(&LUN, outputFile.c_str(), outputFile.size()); break;
			case SKROOTMODE::COPY:  skroot_open_(&LUN, outputFile.c_str(), outputFile.size()); break;
		}
		
		// the following are not relevant for WRITE mode
		if(skrootMode!=SKROOTMODE::WRITE){
			// set the input file(s) to read.
			// These need to be paths to files - glob patterns are not supported.
			// Only the first will be used to retrieve the RareList, but multiple calls can
			// be made to add subsequent files to the TChain.
			// This doesn't yet open the files, it just adds them to a list of strings.
			for(std::string& fname_in : list_of_files){
				skroot_set_input_file_(&LUN, fname_in.c_str(), fname_in.size());
			}
			
			// disable unused input branches.
			int io_dir = 0; // 0 for input branch, 1 for output branch
			if(ActiveInputBranches.size() &&
				std::find(ActiveInputBranches.begin(),ActiveInputBranches.end(),"*")==ActiveInputBranches.end()){
				for(auto&& abranch : default_branches){
					if(std::find(ActiveInputBranches.begin(),ActiveInputBranches.end(),abranch) ==
						ActiveInputBranches.end()){
						skroot_zero_branch_(&LUN, &io_dir, abranch.c_str(), abranch.size());
					}
				}
			}
		}
		// disable unwanted output branches. Only applicable to COPY mode
		if(skrootMode==SKROOTMODE::COPY){
			int io_dir = 1;
			if(ActiveOutputBranches.size() &&
				std::find(ActiveOutputBranches.begin(),ActiveOutputBranches.end(),"*")==ActiveOutputBranches.end()){
				for(auto&& abranch : default_branches){
					if(std::find(ActiveOutputBranches.begin(),ActiveOutputBranches.end(),abranch) ==
						ActiveOutputBranches.end()){
						skroot_zero_branch_(&LUN, &io_dir, abranch.c_str(), abranch.size());
					}
				}
			}
		}
		
		// ok now we perform the actual file opening, TTree cloning, and branch address setting.
		// except in the case of 'write', where all necessary steps are done on construction
		if(skrootMode!=SKROOTMODE::WRITE){
			skroot_init_(&LUN);
		}
		
		// we can access the manager via:
		//TreeManager* mgr = skroot_get_mgr_(&LUN);
		
		// the primary reason to use the TreeManager is to populate the fortran common blocks
		// required by many old SK algorithms. This is not done by the TreeManager itself,
		// but by the fortran routines `skread` and `skrawread`. These internally read data
		// from the ROOT files via the skroot_get_* functions (e.g. via headsk.F),
		// and so depend on there being an underlying TreeManager.
		// There are a few extra steps we need to do to initialize up the fortran common blocks...
		
		// "initialize data structures" it says..? We need this even for ROOT files.
		kzinit_();
		
		// Warning!
		// XXX i'm not sure what of the following (if any) should be called in 'write' mode XXX
		// So... modify if needed. Please explain any changes in comments.
		
		// options for what to read etc.
		skoptn_(const_cast<char*>(skroot_options.c_str()), skroot_options.size());
		
		// options for masking OD / dead / noisy channels
		skbadopt_(&skroot_badopt);
		
		// need to set skgeometry in skheadg common block
		skheadg_.sk_geometry = sk_geometry;
		geoset_();
		
		if(skrootMode!=SKROOTMODE::WRITE){
			// after opening the file with a TreeManager we should still be able to access it
			// with an MTreeReader in parallel, should we want a nicer interface...
			Log(toolName+" creating MTreeReader in parallel to TreeManager",v_debug,verbosity);
			TTree* intree = skroot_get_tree(&LUN);  // n.b. no trailing underscore for this one
			get_ok = myTreeReader.Load(intree);
			if(not get_ok){
				Log(toolName+" failed to open reader on tree "+treeName,v_error,verbosity);
				return false;
			}
			
			bool isMC;
			// skread/tqreal etc detect whether a file is MC or data by checking
			// the `mdrnsk` ("SK run mode" ...) member of the HEADER branch.
			// Unfortunately this is only populated in physics event entries
			// (not pedestal or status entries), which means we need to scan
			// the file for the first suitable entry before we can do the check.
			/*
			// The below works, but it spits out: "error reading the root tree"
			// until it finds a data event. Safe to ignore, but misleading...
			while(true){
				skroot_next_entry_(&LUN,&get_ok);
				if(get_ok){
					std::cerr<<"hit end of tree in checking for MC!"<<std::endl;
					break;
				}
				skcread_(&LUN, &get_ok);
				std::cout<<"entry returned "<<get_ok<<std::endl;
				if(get_ok<1) break;
			}
			const Header* header;
			myTreeReader.Get("HEADER", header);
			if(header==nullptr){
				Log(toolName+" failed to get header of first data event entry!",v_error,verbosity);
			} else {
				isMC = (header->mdrnsk==0 || header->mdrnsk==999999);
			}
			*/
			// We could instead just copy the checks for PDST / RUNINFO entries
			// from the top of headsk.F and avoid using skcread.
			int tmp_entry=0;
			while(true){
				get_ok = myTreeReader.GetEntry(tmp_entry);
				if(get_ok<=0){
					Log(toolName+" error! Hit end of tree while checking if MC!",v_error,verbosity);
					break;
				}
				const Header* header;
				myTreeReader.Get("HEADER", header);
				get_ok &= ((header->nrunsk==0 && header->mdrnsk!=0 && header->mdrnsk!= 999999) ||
				           (std::bitset<8*sizeof(int)>(header->ifevsk).test(19)) ||
				           (header->nrunsk==0 && header->sk_geometry==0 && header->ifevsk!=0) );
				if(get_ok==0){
					// physics entry!
					isMC = (header->mdrnsk==0 || header->mdrnsk==999999);
					break;
				}
				++tmp_entry;
			}
			
			/*
			// A much simpler option is to check for the 'MC' branch, but this could be:
			// 1) dropped during analysis, if the user decides it's not needed
			// 2) created but not filled, if the file is made with skroot_open_write
			// so that's not perfect either.
			isMC = (intree->FindBranch("MC")!=nullptr);
			*/
			
			// put this MC flag into the MTreeReader
			myTreeReader.SetMCFlag(isMC);
			
			// if we're processing MC, we should probably only call `skread`. If we're processing
			// data, we probably need to call `skrawread` as well. (see ReadEntry for more info).
			// We'll use this as a default, but allow ourselves to be overridden by a user.
			if(skreadUser==0){
				skreadMode = (isMC) ? 0 : 2;  // 0=skread only, 1=skrawread only, 2=both
			} else {
				skreadMode = skreadUser-1;
			}
		}
		
	} else {
		
		Log(toolName+" creating MTreeReader to read tree "+treeName,v_debug,verbosity);
		
		get_ok = myTreeReader.Load(list_of_files, treeName);
		if(not get_ok){
			Log(toolName+" failed to open reader on tree "+treeName,v_error,verbosity);
			return false;
		}
		
		// for efficiency of reading, only enable used branches
		Log(toolName+" activating branches",v_debug,verbosity);
		if(ActiveInputBranches.size() && 
			std::find(ActiveInputBranches.begin(), ActiveInputBranches.end(), "*")==ActiveInputBranches.end()){
			// only disable unlisted branches if we have a non-empty list of active branches
			// and the key "*" was not specified.
			myTreeReader.OnlyEnableBranches(ActiveInputBranches);
		}
	}
	
	// put the reader into the DataModel
	Log(toolName+" registering tree reader "+readerName,v_debug,verbosity);
	m_data->Trees.emplace(readerName,&myTreeReader);
	
	// get first entry to process
	if(firstEntry>0) entrynum = firstEntry;
	
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
		} while(entrynum>0 && entrynum<firstEntry);
		if(entrynum<0){
			Log(toolName+" was given both a selections file and a firstEntry,"
				+" but no passing entries were found after the specified starting entry!",v_error,verbosity);
			return false;
		}
		Log(toolName+" reading from entry "+toString(entrynum),v_debug,verbosity);
	}
	
	return true;
}

bool TreeReader::Execute(){
	
	// nothing to do in write mode
	if(skroot==1 && skrootMode==SKROOTMODE::WRITE) return true;
	
	Log(toolName+" getting entry "+toString(entrynum),v_debug,verbosity);
	
	// For SKROOT files most entries are pedestal or status entries that don't actually
	// contain detector data relating to a physics event. We'll usually want to skip these,
	// so if requested (on by default) keep reading until we get an event entry.
	do {
		// load next entry
		get_ok = ReadEntry(entrynum);
		
		// get the index of the next entry to read.
		if(myTreeSelections==nullptr){
			entrynum++;
		} else {
			entrynum = myTreeSelections->GetNextEntry(cutName);
		}
	} while(skroot && skip_ped_evts && get_ok==-99); // skip if this is a pedestal/status entry
	
	++readEntries;  // keep track of the number of entries we've actually returned
	
	// check if we've hit the user-requested limit on number of entries to read
	if((maxEntries>0)&&(readEntries>=maxEntries)){
		Log(toolName+" hit max events, setting StopLoop",v_message,verbosity);
		m_data->vars.Set("StopLoop",1);
	}
	// use LoadTree to check if the next entry is valid without loading it
	// (this checks whether we've hit the end of the TTree/TChain)
	else if(myTreeReader.GetTree()->LoadTree(entrynum)<0){
		Log(toolName+" reached end of TTree, setting StopLoop",v_message,verbosity);
		m_data->vars.Set("StopLoop",1);
	}
	
	return true;
}

bool TreeReader::Finalise(){
	
	if(myTreeSelections) delete myTreeSelections;
	
	if(skroot){
		CloseLUN();                 // deletes tree, file, TreeManager.
		myTreeReader.SetClosed();   // tell the MTreeReader the file is closed.
	}
	// otherwise the MTreeReader destructor will close the file.
	
	// We could check our lunlist to see if there are any remaining TreeReaders,
	// and if this is the last one, clean up the SuperManager.
	// The issue with that is we can't check whether the SuperManager has any
	// outstanding TreeManagers, so if a user has done skroot_open_* themselves,
	// (and didn't create a corresponding entry in the lunlist)
	// then deleting the SuperManager would orphan the corresponding TreeManagers.
	// The best option is probably to just leave this to the OS - i.e. let it leak.
	//skroot_end_();
	
	return true;
}

int TreeReader::ReadEntry(long entry_number){
	int bytesread=1;
	// load next entry data from TTree
	if(skroot){
		
		// Populating fortran common blocks with SKROOT entry data requires using
		// SKRAWREAD and/or SKREAD.
		// These functions call various skroot_get_* functions to retrieve branch data.
		// They then use that data to populate the fortran common blocks.
		// * SKRAWREAD only works on data, not MC
		// * They both only call SOME branch getters
		
		// supposedly skrawread only reads the head, tq and pedestal branches < check this
		// Supposedly skrawread loads some "constant tables" ...
		// These seem to be needed to call e.g. bonsai on data; if only skread
		// is called bonsai complains `neighbsk: error - ISEQ is <0 or > MAXPM! 0`
		// I have no idea what skread does.
		// TODO: Find out exactly what both of these functions do!
		// what do they read, calculate, and populate, along with any side-effects.
		// what does or doesn't get invoked in SKREAD when LUN is negative, and what
		// does this mean for the user. Is there an equivalent behaviour with skrawread?
		
		//   This will result in a segfault in SKROOT copy mode unless ALL remaining
		//   branches are either disabled or loaded by the user (e.g. via skroot_get_entry)
		
		//                         *** IMPORTANT ***
		// Both SKRAWREAD and SKREAD will invoke `skroot_next_entry` before data retrieval
		// *IF* the passed LUN is >0! This increments an entry number member of
		// the TreeManager. Subsequent `skroot_get_*` calls will load data from the
		// TTree entry corresponding to the internal entry number.
		// This allows one to call both SKRAWREAD and SKREAD, but ONLY in the correct
		// order and if SKREAD is given a negative LUN!
		// The sign of LUN also impacts whether various initialization stuff
		// is done, at least in SKREAD (not sure about SKRAWREAD).
		// So you MUST set LUN to be positive for ONE AND ONLY ONE CALL...probably.
		
		//                       *** VERY IMPORTANT ***
		// `skroot_jump_entry_(N) will actually set the internal index to (N-1)!
		// One needs to be very careful with calls to `skroot_jump_entry_`, `skrawread` and `skread`
		// to ensure the entry being read is indeed the entry you want!
		//                      ************************
		
		// =============================================================
		
		// first, set the internal entry number of the TreeManager.
		// remember this will actually set it to (entry_number - 1)
		// but we will then use either SKRAWREAD or SKREAD with a positive LUN
		// to increment it to the one we actually want, while also doing any
		// necessary internal reinitialization of these functions (presumably)
		int entry_temp = static_cast<int>(entry_number);
		skroot_jump_entry_(&LUN, &entry_temp, &get_ok);
		if(get_ok==1){
			// ran off end of TTree!
			bytesread=0;
		}
		
		// skreadMode: 0=skread only, 1=skrawread only, 2=both
		if(bytesread>0 && skreadMode>0){
			skcrawread_(&LUN, &get_ok); // N.B. positive LUN (see above)
			if(get_ok==1){
				Log(toolName+" read error "+toString(get_ok)+" calling skcrawread ",v_error,verbosity);
				// lf_allfit actually continues the read loop if this is encountered,
				// so perhaps this is a recoverable error, or just an error relating to this entry?
				// FIXME if so it may be better to continue to next entry instead of bailing
				bytesread = -1;
			} else if(get_ok==2){
				// this just indicates we've reached the end of the file
				// if we've hit this, it's not good, because downstream tools
				// won't have any valid data!
				bytesread = 0;
			} else if(get_ok!=0) {
				// pedestal or status entry, no detector data, not an actual event
				// 3 = pedestal entry, 4 = runinfo entry.
				//Log(toolName+" skrawread pedestal or status event, skipping",v_debug,verbosity);
				// this happens a lot...
				bytesread = -99;
			}
		}
		if(bytesread>0 && skreadMode!=1){  // skip skread if skrawread had an error
			int LUN2 = LUN;
			if(skreadMode==2) LUN2 = -LUN;  // if we already called skrawread, use a negative LUN
			skcread_(&LUN2, &get_ok);
			if(get_ok==1){
				// error reading entry
				bytesread = -1;
			} else if(get_ok==2){
				// end of file
				bytesread = 0;
			} else if(get_ok!=0) {
				// pedestal or status entry
				bytesread = -99;
			}
		}
		
		// =============================================================
		
		// As mentioned above, neither of these load all TTree branches.
		// To do that we need to call skroot_get_entry.
		skroot_get_entry_(&LUN);
		
	} else {
		// else not using TreeManagers / skread / etc
		bytesread = myTreeReader.GetEntry(entry_number);
	}
	
	// stop loop if we ran off the end of the tree
	if(bytesread==0){
		// not good because downstream tools will not have valid data!
		// we should protect against this in Execute() though.
		Log(toolName+" hit end of input file, stopping loop",v_warning,verbosity);
		m_data->vars.Set("StopLoop",1);
	} else if(bytesread==-99){
		//Log(toolName+" skrawread pedestal or status event, skipping",v_debug,verbosity);
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
	
	bool settingInputBranchNames=false;
	bool settingOutputBranchNames=false;
	
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
		if (thekey=="StartInputBranchList"){
			settingInputBranchNames = true;
		}
		else if(thekey=="EndInputBranchList"){
			settingInputBranchNames = false;
		}
		else if(settingInputBranchNames){
			ActiveInputBranches.push_back(Line);
		}
		else if (thekey=="StartOutputBranchList"){
			settingOutputBranchNames = true;
		}
		else if(thekey=="EndOutputBranchList"){
			settingOutputBranchNames = false;
		}
		else if(settingOutputBranchNames){
			ActiveOutputBranches.push_back(Line);
		}
		else if(thekey=="verbosity") verbosity = stoi(thevalue);
		else if(thekey=="inputFile") inputFile = thevalue;
		else if(thekey=="outputFile") outputFile = thevalue; // when using SKROOT copy mode
		else if(thekey=="FileListName") FileListName = thevalue;
		else if(thekey=="treeName") treeName = thevalue;
		else if(thekey=="readerName") readerName = thevalue;
		else if(thekey=="firstEntry") firstEntry = stoi(thevalue);
		else if(thekey=="maxEntries") maxEntries = stoi(thevalue);
		else if(thekey=="selectionsFile") selectionsFile = thevalue;
		else if(thekey=="cutName") cutName = thevalue;
		else if(thekey=="skroot") skroot = stoi(thevalue);
		else if(thekey=="skrootMode") skrootMode = SKROOTMODE(stoi(thevalue));
		else if(thekey=="skreadMode") skreadUser = stoi(thevalue);
		else if(thekey=="LUN") LUN = stoi(thevalue);
		else if(thekey=="skoptn") skroot_options = thevalue;
		else if(thekey=="skbadopt") skroot_badopt = stoi(thevalue);
		else if(thekey=="SK_GEOMETRY") sk_geometry = stoi(thevalue);
		else if(thekey=="skipPedestals") skip_ped_evts = stoi(thevalue); // is this redundant with skoptn?
		else {
			Log(toolName+" error parsing config file line: \""+LineCopy
				+"\" - unrecognised variable \""+thekey+"\"",v_error,verbosity);
		}
	}
	
	// done parsing, close config file
	fin.close();
	
	return 1;
}

// return a new LUN. We accept a hint, but will only apply it if not already assigned.
int TreeReader::GenerateNewLUN(){
	// each LUN (logic unit number, a fortran file handle (ID) and/or an ID used
	// by the SuperManager to identify the TreeManager associated with a file) must be unique.
	std::map<std::string,int> lunlist;
	m_data->CStore.Get("LUNList",lunlist);
	// check if this LUN is free, otherwise print a warning and assign a different LUN
	// since we map names to LUNs, to check if a LUN is free we need to scan by value, not by key.
	// easiest way to do this while also sorting by LUN is to reverse the map.
	// sorting allows us to immediately know the next free LUN in case this one is in use.
	std::map<int,std::string> revlist;
	for(auto&& apair : lunlist){ revlist.emplace(apair.second,apair.first); }
	if(revlist.count(LUN)){
		int reqLUN=LUN;
		LUN = revlist.rbegin()->first; // get the last assigned LUN
		++LUN;                         // we'll use the next one
		Log(toolName+": Warning! Cannot assign LUN "+toString(LUN)
			+" as it is already taken. Assigning "+toString(reqLUN)+" instead.",v_warning,verbosity);
	}
	//assign the LUN
	lunlist.emplace(readerName,LUN);
	m_data->CStore.Set("LUNList",lunlist);
	return LUN;
}

void TreeReader::CloseLUN(){
	
	// get the map of open TreeManagers
	std::map<std::string,int> lunlist;
	m_data->CStore.Get("LUNList",lunlist);
	
	// sanity check, although all we can do is print warnings if something's gone wrong
	if(lunlist.count(readerName)==0){
		Log(toolName+" error! Could not find LUN list entry for reader "+readerName,v_error,verbosity);
	} else if(lunlist.at(readerName)!=LUN){
		Log(toolName+" error! LUN entry in lunlist is not the same as that in its TreeReader?!",v_error,verbosity);
		// now which do we call close on? Should we erase the lunlist entry? Should not happen!
		lunlist.erase(readerName);
	} else {
		lunlist.erase(readerName);
	}
	
	// close input skroot files, delete the TTreeManager
	// the SuperManager has a nullptr check so shouldn't seg even if this LUN is invalid
	skroot_close_(&LUN);
	// note that if the LUN wasn't valid, that LUN will now be created in the SuperManager's map
	// with a nullptr entry, which completely locks the LUN - it doesn't point to a valid TreeManager,
	// and it's not possible to remove it from the list. TODO fix the SuperManager.
	
}
