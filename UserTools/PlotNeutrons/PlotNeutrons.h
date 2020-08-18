/* vim:set noexpandtab tabstop=4 wrap */
#ifndef PlotNeutrons_H
#define PlotNeutrons_H

#include <string>
#include <iostream>
#include <vector>
#include <array>

#include "Tool.h"

class TApplication;

/**
* \class PlotNeutrons
*
* A simple tool to plot basic characteristics of neutron captures, reading outputs from skdetsim.
* 1. run skdetsim or skdetsim-gd with output option 'SKCNTL-OUTPUTTYPE 2' to generate zbs output file
* 2. convert the zbs file to an hbk file with `fillnt_simple.sh -o (output hbook file) (input zbs file)`
* 3. convert the hbk file to a root file with `h2root (input hbook file) (output root file)`
* This somewhat convoluted process is required to retain information about neutrons and gammas.
* This file processes the output root files from step 3.
*
* $Author: M.O'Flahery $
* $Date: 2020/08/12 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class PlotNeutrons: public Tool {
	
	public:
	
	PlotNeutrons(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Executre function used to perform Tool purpose.
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	
	private:
	std::string toolName;
	bool DrawPlots;
	TApplication* rootTApp;
	
	// verbosity levels: if 'verbosity' < this level, the message type will be logged.
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	std::string logmessage="";
	int get_ok=0;
	
	// file stuff
	std::string inputFile;   // if just passing a single filename directly to this tool
	std::vector<std::string> input_file_names;  // if using upstream LoadFileList tool
	
	// variables to read in
	// ====================
	int MAX_N_PRIMARIES=255;   // starting value, check tot_n_primaries and/or n_outgoing_primaries
	int MAX_N_SECONDARIES=255; // starting value, check n_secondaries_1 and/or n_secondaries_2
	
	double water_transparency;                                  // [cm]
	
	// event meta info
	int run_number;
	int subrun_number;
	int event_number;
//	std::array<int,3> event_date;                               // [year,month,day] array
//	std::array<int,4> time;                                     // [hour,minute,second,???] array
	
	// event level detector info
	int N_hit_ID_PMTs;                                          // "nqisk"
	int total_ID_pes;                                           // "qismsk"
	int max_ID_PMT_pes;                                         // "qimxsk", max # PEs from a single ID PMT?
	
	// neutrino interaction info
//	int nu_intx_mode;                                           // use constants::neut_mode_to_string to decode
//	int tot_n_primaries;                                        // ingoing and outgoing
//	// the following ar aarrays of size tot_n_primaries
//	std::vector<int> primary_pdg;                               // see constants::numnu_code_to_string for mapping
//	std::vector<std::array<double,3>> primary_momentum;         // [GeV/c]
	
	// primary event
	std::array<double,3> primary_event_vertex;                  // [cm]
	double primary_event_dist_from_wall;                        // [cm]
	int n_outgoing_primaries;                                   // should be (tot_n_primaries - 2)...
	// following are arrays of size n_outgoing_primaries
	std::vector<int> primary_G3_code;                           // use constants::g3_to_pdg to map to pdg code
	std::vector<std::array<double,3>> primary_start_mom;        // [units?] this and ipv are arrays of size npar
	
	// secondaries - first secondaries arrays...
	int n_secondaries_1;
	// the following are arrays of size npar2
	std::vector<int> secondary_G3_code_1;                       // 
	std::vector<std::array<double,3>> secondary_start_vertex_1; // array of 3, [cm?] what about time?
	std::vector<double> secondary_start_dist_from_wall_1;       // [cm?]
	std::vector<std::array<double,3>> secondary_start_mom_1;    // [units?]
	std::vector<int> secondary_origin_1;                        // what is "origin"?
	
	// secondaries - second secondaries array...
	int n_secondaries_2;
	// the following are arrays of size nscndprt
	std::vector<int> secondary_PDG_code_2;                      //
	std::vector<std::array<double,3>> secondary_start_vertex_2; // [units?]
	std::vector<double> secondary_start_time_2;                 // [ns]? relative to event start?
	std::vector<std::array<double,3>> secondary_start_mom_2;    // [units?]
	std::vector<int> secondary_gen_process;                     // use constants::G3_process_code_to_string
	std::vector<int> secondary_n_daughters;                     // 
	std::vector<int> secondary_first_daugher_index;             // if >0, 1-based index in this array
	std::vector<int> parent_index;                              // if >0, 1-based index in this array
	
	// further parentage information - Useful?
//	std::vector<int> parent_G3_code;                                 // or is it a PDG code?
//	std::vector<std::array<double,3>> parent_mom_at_sec_creation;    // use w/daughter γ to see n energy @ capture
//	std::vector<std::array<double,3>> parent_init_pos;               // [cm?]
//	std::vector<std::array<double,3>> parent_init_mom;               // [units?]
//	std::vector<int> parent_G3_trackid;                              // how do we use this?
//	std::vector<int> parent_G3_stack_trackid;                        // how do we use this?
//	std::vector<int> parent_trackid;                                 // how do we use this?
//	std::vector<int> parent_track_pid_code;                          // i'm so confused
	
	

};


#endif