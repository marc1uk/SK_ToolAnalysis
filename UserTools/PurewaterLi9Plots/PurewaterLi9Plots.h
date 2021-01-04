#ifndef PurewaterLi9Plots_H
#define PurewaterLi9Plots_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "Algorithms.h"
#include "MTreeReader.h"
#include "MTreeSelection.h"
#include "SkrootHeaders.h"    // MCInfo, Header etc.
#include "thirdredvars.h"     // ThirdRed class

class TH1F;
class TF1;
class THStack;

/**
* \class PurewaterLi9Plots
*
* This tool makes plots of the rate of spallation products, tries to select a li9 enhanced sample,
* and fit the rate of li9 decay and associated neutron capture times.
* Ultimately the number of spallation and separately of li9 events are fit to obtain
* the number of events, from which the rate of isotope generation is extracted.
* Event selection is done by the PurewaterLi9Rate tool. This one makes and fits the plots.
*
* $Author: M.O'Flaherty $
* $Date: 2019/05/28 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class PurewaterLi9Plots: public Tool {
	
	public:
	
	PurewaterLi9Plots(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Executre function used to perform Tool purpose.
	bool Finalise(); ///< Finalise funciton used to clean up resources.
	
	private:
	
	// file IO stuff
	// =============
	std::string inputFile="";                                   // input file to read
	std::string cutFile="";                                     // file of entries to use
	MTreeReader* myTreeReader=nullptr;                          // the TTree reader
	MTreeSelection* myTreeSelections=nullptr;                   // record what passes what cuts
	int entry_number=0;                                         // input TTree entry
	int first_entry=0;                                          // first TTree entry to start from
	int MAX_ENTRIES=-1;                                         // max num input TTree entries to process
	uint64_t entries_processed=0;                               // since processed entry numbers aren't sequential
	std::string outputFile="li9_plots.root";                    // output file to write
	bool driver=false;                                          // whether we're controlling the ToolChain
	std::string topCutName;                                     // name of loosest cut we need every entry from
	std::string valuesFile="li9_values.bs";                     // save or read values from BoostStore
	std::string valuesFileMode="";                              // "read", "write" or anything else = neither
	int show_plots=0;                                           // if we want to show plots, not just save them
	
	// output file of plots and fits
	TFile* outfile=nullptr;                                     // the output TFile
	// output file of saved values
	BoostStore* valueStore=nullptr;
	
	int DisableUnusedBranches();
	int ReadEntryNtuple(long entry_number);
	
	bool CreateOutputFile(std::string outputFile);
	void CloseFile();
	
	bool Analyse();        // main body
	bool PlotMuonDt();
	bool PlotMuonDlt();
	bool PlotPaperDt(THStack& ourplots);
	bool PlotPaperDlt(THStack& ourplots);
	bool PlotSpallationDt();
	bool FitSpallationDt(TH1F& dt_mu_lowe_hist_short, int rangenum);
	bool PlotLi9TripletsDt();
	bool BinnedLi9DtChi2Fit(TH1F* li9_muon_dt_hist);
	bool UnbinnedLi9DtLogLikeFit();
	double li9_lifetime_loglike(double* x, double* par);
	bool PlotLi9NtagDt();
	bool PlotLi9BetaEnergy();
//	void MeasureDltSystematic(); // TODO
//	void MeasureIsotopeRates();  // TODO
	
	// helper functions used in FitSpallationDt
	void SetParameterNames(std::string isotope);
	void SetParameterNames(std::string isotope_1, std::string isotope_2);
	void FixLifetime(TF1& func, std::string isotope);
	void PushFitAmp(TF1& func, std::string isotope);
	void PullFitAmp(TF1& func, std::string isotope);
	TF1 BuildFunction(std::vector<std::string> isotopes, double func_min=0, double func_max=30);
	
	// variables to read in
	// ====================
	// TODO reduce to just what's needed
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
	
	// histogram ranges
	// ==================
	float li9_lifetime_dtmin;         // range of dt_mu_lowe values to accept
	float li9_lifetime_dtmax;         // for Li9 candidates, [us]
	
	// variables to save/plot
	// ======================
	// maps of muboy class vs a histogram of mu->lowe time and transverse distance
	std::vector<std::vector<float>> dlt_vals_pre{6};
	std::vector<std::vector<float>> dt_vals_pre{6};
	std::vector<std::vector<float>> dlt_vals_post{6};
	std::vector<std::vector<float>> dt_vals_post{6};
	
	// varying muon->lowe dt thresholds for assessing systematic of dlt cut
	int num_dt_cuts=5; // this is the size of the spall_lifetimes vector in PurewaterLi9Rate tool
	// moreover it's the number of cuts named "pre/post_mu_dt_cut_%d" we have.
	// TODO retrieve the list of cuts from the MTreeSelection, count how many we have of this type?
	std::vector<std::vector<float>> dlt_systematic_dt_cuts_pre{5};
	std::vector<std::vector<float>> dlt_systematic_dt_cuts_post{5};
	// each entry is a different dt cut, inner vector is the dlts of passing events
	// for a given dt cut, the difference between values gives the distribution of *spallation* dlt.
	// across the various dt cuts, the difference between spallation dlt distributions gives
	// the systematic error in dlt cut efficiency.
	// We're particularly interested in the variation in the bin corresponding to cut value of dlt=200cm,
	// so make sure we have a bin edge at that value
	
	// distribution of all muon->lowe times, after dlt cut, to extract rates of isotopes
	std::vector<float> dt_mu_lowe_vals;
	// distribution of muon->lowe times after further Li9 cuts
	std::vector<float> li9_muon_dt_vals;
	// distribution of li9->ntag times for mu+li9+ntag triplets
	std::vector<float> li9_ntag_dt_vals;
	// distribution of lowe energies for mu+li9+ntag triplets
	std::vector<float> li9_e_vals;
	
	// functions used during fitting of the number of isotope events
	std::map<std::string,TF1*> fitfuncs;
	
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
	
};


#endif
