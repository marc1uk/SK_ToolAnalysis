#ifndef PurewaterLi9Rate_H
#define PurewaterLi9Rate_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "Algorithms.h"
#include "MTreeReader.h"
#include "MTreeSelection.h"
#include "SkrootHeaders.h"    // MCInfo, Header etc.
#include "thirdredvars.h"     // ThirdRed class

class TSpline3;

//#include "cut_third.h"

/**
* \class PurewaterLi9Rate
*
* Measure rate of Li9 production by extracting the number of Li9 events following muons
* Methods based on 2015 paper by Yang Zhang
*
* $Author: M.O'Flaherty $
* $Date: 2020/12/11 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/

class PurewaterLi9Rate: public Tool {
	
	public:
	
	PurewaterLi9Rate(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Executre function used to perform Tool perpose.
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	
	private:
	
	// file IO stuff
	// =============
	std::string inputFile="";                                   // input file to read
	MTreeReader myTreeReader;                                   // the TTree reader
	MTreeSelection myTreeSelections;                            // record what passes what cuts
	int entry_number=0;                                         // input TTree entry
	int first_entry=0;                                          // first TTree entry to start from
	int MAX_ENTRIES=-1;                                         // max num input TTree entries to process
	std::string outputFile="li9_cuts.root";                     // output file to write
	
	// these are unused
	TFile* outfile=nullptr;                                     // the output TFile
	TTree* outtree=nullptr;                                     // the output TTree
	int WRITE_FREQUENCY=10;                                     // update output tree every N fills
	
	// cut configurations
	// ==================
	int run_min=0;
	int run_max=0;
	float max_closest_muon_dt=0.001;  // 1ms
	float max_closest_lowe_dx=490;    // 490cm
	float ntag_FOM_threshold=0.95;    // BDT FOM threshold for Li9 ntagging
	float li9_lifetime_dtmin;         // range of dt_mu_lowe values to accept
	float li9_lifetime_dtmax;         // for Li9 candidates, [us]
	// ROC curve reading for ntag
//	TSpline3 *sig=nullptr;       // TODO for BDT efficiency
//	TSpline3 *sigsys=nullptr;    // TODO for BDT efficiency
//	TSpline3 *bg=nullptr;        // TODO for BDT efficiency
//	TSpline3 *bgsys=nullptr;     // TODO for BDT efficiency
	
	// variables to read in
	// ====================
	const Header  *HEADER = new Header;
	const LoweInfo *LOWE  = new LoweInfo;
	const ThirdRed *thirdredvars = new ThirdRed;
	int num_neutron_candidates;                      // num pulses: i.e. num neutrons in AFT window
	int max_hits_200ns_AFT;                          // max num hits in 200ns sliding window within AFT trigger
	basic_array<float*> ntag_FOM;                    // ntag BDT output figure of merit
	basic_array<float*> dt_lowe_n;                   // dt positron->neutron
	int num_pre_muons;                               // num muons in 30s preceding lowe event
	int num_post_muons;                              // num muons in 30s following lowe event
//	basic_array<float*> spaloglike;                  // log likelihood of being spallation based on...?
//	basic_array<float*> spaloglike_shfld;            // version based on...?
//	basic_array<float*> spaloglike_kirk;             // version based on...?
//	basic_array<float*> spaloglike_kirk_shfld;       // version based on...?
//	basic_array<float*> spaloglike_li9;              // version based on...?
//	basic_array<float*> spaloglike_kirk_li9;         // version based on...?
	basic_array<int*>   mu_class;                    // muboy muon classification (see enum class at top)
//	basic_array<int*>   mubntrack;                   // num muons found by muboy
	basic_array<int*>   mu_index;                    // index of this muon, of those found by muboy
	basic_array<float*> mu_fit_goodness;             // muboy goodness of fit: >0.4 is ok for single thru muons
//	basic_array<float*> mubffgood;                   // brute force fitter goodness of fit
	basic_array<float*> dt_mu_lowe;                  // time between muon and lowe event [seconds]
//	basic_array<float*> spadt_li9;                   // version based on...?
	basic_array<float*> dlt_mu_lowe;                 // transverse distance between muon and lowe event [cm]
//	basic_array<float*> spadll;                      // longitudinal distance between muon and lowe event [cm]
//	basic_array<float*> spadll_kirk;                 // version based on...?
//	basic_array<float*> sparesq;                     // surplus charge over that of a MIP with reco track length
//	basic_array<float*> spaqpeak;                    // ? 
//	basic_array<float*> spaqpeak_kirk;               // version based on...?
//	basic_array<float*> spamuqismsk;                 // ?
//	basic_array<float*> spamuqismsk_pertrack;        // ?
//	basic_array<float*> spadts;                      // ?
//	basic_array<int*>   candidates;                  // ?
//	basic_array<int*>   muindex;                     // ?
//	basic_array<int*>   neut_flag;                   // ?
//	basic_array<int*>   mult_flag;                   // ?
//	basic_array<float*[3]> neut_shift;               // ?
//	basic_array<float*> neutdiff;                    // ?
	basic_array<float*> closest_lowe_60s;            // closest distance to another lowe event within 60s??
	
	// variables to write out
	// ======================
	// output is handled by the MTreeSelection
	
	// calculation variables
	// =====================
	// varying BDT FOM thresholds
	std::string bdt_outfile;
	int num_bdt_cut_bins;
	float bdt_cut_min, bdt_cut_max;
	std::vector<float> bdt_cut_thresholds;
	int nspatot=0; // sonia's count of li9 candidates before ntag
	std::vector<float> dtimes_li9_ncap;  // sonia's li9 ncapture times
	
	std::vector<float> spall_lifetimes;
	
	// functions
	// =========
	// input
	int DisableUnusedBranches();                                // disable unused branches to speed up reading input
	int ReadEntryNtuple(long entry_number);                     // get next entry from TreeReader
	
	// output - these are currently redundant, output is handled by the MTreeSelection
	int CreateOutputFile(std::string outputFile);               // create output TFile and TTree
	int FillTree();                                             // fill output TTree entries
	void ClearOutputTreeBranches();                             // reset containers to prevent carry-over
	int WriteTree();                                            // update output TFile
	void CloseFile();                                           // close output TFile
	
	// other
	void fill_ntag_roc(TSpline3 **sig, TSpline3 **bg, TSpline3 **sigsys, TSpline3 **bgsys);
	void make_BDT_bins();
	bool apply_third_reduction(const ThirdRed *th, const LoweInfo *LOWE);
	bool Analyse();        // main body
	bool SoniasAnalyse();  // translation of sonia's set of cuts
	
	//void PrintBranches();                                     // print an output event, for debug
	
	// standard tool stuff
	// ===================
	std::string toolName;
	// verbosity levels: if 'verbosity' < this level, the message type will be logged.
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	std::string logmessage="";
	int get_ok=0;
	
	// for testing performance
	bool stopping=false;
	uint64_t num_processed_events=0;
	std::chrono::high_resolution_clock::time_point toolchain_start;
	std::chrono::high_resolution_clock::time_point toolchain_end;
	double loop_times[1000];
	double analyse_times[1000];
	int loop_i=0;
};


#endif
