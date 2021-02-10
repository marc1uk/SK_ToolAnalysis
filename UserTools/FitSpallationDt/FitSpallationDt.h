/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef FitSpallationDt_H
#define FitSpallationDt_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "MTreeReader.h"
#include "SkrootHeaders.h" // MCInfo, Header etc.

#include "ColourWheel.h"

class TH1F;
class TF1;

/**
* \class FitSpallationDt
*
* FIXME
*
* $Author: M.O'Flaherty $
* $Date: 2019/05/28 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class FitSpallationDt: public Tool {
	
	public:
	FitSpallationDt();         ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute();   ///< Execute function used to perform Tool purpose.
	bool Finalise();  ///< Finalise funciton used to clean up resources.
	
	private:
	// functions
	// =========
	int ReadEntry(long entry_number);
	int GetBranches();
	bool Analyse();
	int DisableUnusedBranches();
	bool GetEnergyCutEfficiencies();
	bool PlotSpallationDt();
	bool FitDtDistribution(TH1F& dt_mu_lowe_hist_short, TH1F& dt_mu_lowe_hist_log, int rangenum);
	// helper functions used in FitSpallationDt
	std::vector<double> MakeLogBins(double xmin, double xmax, int nbins);
	void FixLifetime(TF1& func, std::string isotope);
	void PushFitAmp(TF1& func, std::string isotope);
	void PushFitAmp(double amp, std::string isotope);
	void PullFitAmp(TF1& func, std::string isotope, bool fix=true);
	void PullPaperAmp(TF1& func, std::string isotope, bool threshold_scaling=true, double fixed_scaling=1.);
	double GetPaperAmp(std::string isotope, bool threshold_scaling, double fixed_scaling);
	void BuildPaperPlot();
	TF1 BuildFunction(std::vector<std::string> isotopes, double func_min=0, double func_max=30);
	TF1 BuildFunction2(std::vector<std::string> isotopes, double func_min=0, double func_max=30);
	
	// tool variables
	// ==============
	std::string toolName;
	std::string outputFile="";
	std::string inputFile;
	std::string treeName;
	int maxEvents;
	double paper_scaling = 1.;
	std::string efficienciesFile="FlukaBetaEfficiencies.bs";    // BoostStore of efficiencies of energy thresholds
	
	MTreeReader myTreeReader; // the TTree reader
	int entrynum=0;
	float dt;                // this is the only variable we need to read
	std::vector<float> my_dt_mu_lowe_vals;
	double livetime=1;
	
	// energy threshold comparison
	// ===========================
	std::map<std::string, double> true_effs_6mev;
	std::map<std::string, double> true_effs_8mev;
	std::map<std::string, double> true_effs_scaling;
	std::map<std::string, double> reco_effs_6mev;
	std::map<std::string, double> reco_effs_8mev;
	std::map<std::string, double> reco_effs_scaling;
	
	// results used in fitting of the number of isotope events
	std::map<std::string,double> fit_amps;
	
	ColourWheel colourwheel;
	
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
	
	// variables to write out
	// ======================
	
};


#endif
