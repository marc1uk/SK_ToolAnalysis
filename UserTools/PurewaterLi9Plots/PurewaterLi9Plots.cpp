/* vim:set noexpandtab tabstop=4 wrap */
#include "PurewaterLi9Plots.h"

#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include "TH1.h"
#include "TGraph.h"
#include "THStack.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMath.h"
#include "TString.h"
#include "TSpline.h"
#include "TEntryList.h"

// For Fitting
#include "Fit/Fitter.h"
//#include "Fit/BinData.h"
#include "Fit/UnBinData.h"
//#include "Fit/Chi2FCN.h"
//#include "Fit/DataOptions.h"
#include "Fit/FitConfig.h"

// For defining the functions
//#include "TList.h"
#include "Math/WrappedMultiTF1.h"
//#include "HFitInterface.h"

#include "Constants.h"

const double li9_endpoint = 14.5; // MeV
const std::map<std::string,double> lifetimes{  // from 2015 paper, seconds
	{"11Be",19.9},
	{"16N",10.3},
	{"15C",3.53},
	{"8Li",1.21},
	{"8B",1.11},
	{"16C",1.08},
	{"9Li",0.26},
	{"9C",0.18},
	{"8He",0.17},
	{"12Be",0.034},
	{"12B",0.029},
	{"13B",0.025},
	{"14B",0.02},
	{"12N",0.016},
	{"13O",0.013},
	{"11Li",0.012},
	{"ncapture",204.8E-6}
};

const std::map<std::string,double> papervals{  // from 2015 paper plot at t->0, represents Ni/τi, with τi in units of 0.006s
	{"11Be",28},
	{"16N",450},
	{"15C",80},
	{"8Li",350},
	{"8B",580},
	{"16C",0},
	{"9Li",350},
	{"9C",75},  // paired, half 150
	{"8He",75}, // paired with above
	{"12Be",0},
	{"12B",75000},
	{"13B",0},
	{"14B",0},
	{"12N",25000},
	{"13O",0},
	{"11Li",0},
	{"const",330},
	// for the comparison function
	{"8He_9C",150},
	{"8Li_8B",930}
};

const std::map<std::string,double> papervals2{  // from 2015 paper table II, converted w/ eqn 2, divided by isotope lifetime in units 0.006s to be consistent with above... TODO; remove the scaling by lifetime and use papervals3
	{"11Be",0},
	{"16N",443},
	{"15C",0},
	{"8Li",375},
	{"8B",489},
	{"16C",0},
	{"9Li",346},
	{"9C",0},
	{"8He",0},
	{"12Be",0},
	{"12B",79264},
	{"13B",0},
	{"14B",0},
	{"12N",25094},
	{"13O",0},
	{"11Li",0},
	{"const",330},   // not in table, use value from plot
	// for the comparison function
	{"8He_9C",0},
	{"8Li_8B",432}
};

const std::map<std::string,double> papervals3{  // from 2015 paper table II, converted w/ eqn 2
	{"11Be",0},
	{"16N",759709},
	{"15C",0},
	{"8Li",75532},
	{"8B",90533},
	{"16C",0},
	{"9Li",15002},
	{"9C",0},
	{"8He",0},
	{"12Be",0},
	{"12B",383107},
	{"13B",0},
	{"14B",0},
	{"12N",66917},
	{"13O",0},
	{"11Li",0},
	{"const",330},   // not in table, use value from plot
	// for the comparison function
	{"8He_9C",0},
	{"8Li_8B",465}
};

const std::map<std::string,double> efficiencies{  // from FLUKA
	{"11Be",38.1},
	{"16N",45.0},
	{"15C",31.8},
	{"8Li",42.8},
	{"8B",51.3},
	{"16C",0},
	{"9Li",39.2},
	{"9C",50.2},  // paired, half 150
	{"8He",22.2}, // paired with above
	{"12Be",0},
	{"12B",45.5},
	{"13B",0},
	{"14B",0},
	{"12N",56.2},
	{"13O",0},
	{"11Li",0}
};

PurewaterLi9Plots::PurewaterLi9Plots():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

bool PurewaterLi9Plots::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",verbosity);            // how verbose to be
	m_variables.Get("inputFile",inputFile);            // the data
	m_variables.Get("cutFile",cutFile);                // the passing events
	m_variables.Get("outputFile",outputFile);          // output file to write
	m_variables.Get("maxEntries",MAX_ENTRIES);         // terminate after processing at most this many events*
	m_variables.Get("topCutName",topCutName);          // first cut we need to see every passing event from.
	m_variables.Get("valuesFileMode",valuesFileMode);  // "write", "read", <anything else>="neither"
	m_variables.Get("valuesFile",valuesFile);          // BoostStore file
	m_variables.Get("bdt_outfile",bdt_outfile);        // ROC curves for the ntag BDT
	//* unused when an upstream tool is driving the toolchain
	//* when using values from a BoostStore the ROOT files do not need to be read
	
	// cut thresholds
	m_variables.Get("li9_lifetime_dtmin",li9_lifetime_dtmin);
	m_variables.Get("li9_lifetime_dtmax",li9_lifetime_dtmax);
	m_variables.Get("li9_ncapture_dtmin",ncap_dtmin);
	m_variables.Get("li9_ncapture_dtmax",ncap_dtmax);
	
	// hacky way to scale the total number of events when comparing to the paper plots
	m_variables.Get("paper_scaling",paper_scaling);
	
	// energy threshold efficiencies, from FLUKA
	m_variables.Get("efficienciesFile",efficienciesFile);
	
	// get the TApplication if we want to show plots on the fly
	m_variables.Get("show_plots",show_plots);
	if(show_plots) data.GetTApp(); // we don't really need it, just trigger its creation
	
	// if reading values from a BoostStore, we don't need to do any loops
	if(valuesFileMode=="read"){
		m_data->vars.Set("StopLoop",1); // jump straight to finalise
	} else {
		// search for a TreeReader from an upstream tool if this toolchain is actively performing the cut,
		// (such as PurewaterLi9Rate)
		intptr_t myTreeSelectionsPtr = 0;
		get_ok = m_data->CStore.Get("MTreeSelection",myTreeSelectionsPtr);
		if(not get_ok || (myTreeSelectionsPtr==0)){
			driver=true;
			// no upstream reader, make our own
			// first we'll need a MTreeReader
			myTreeReader = new MTreeReader(inputFile, "data");
			DisableUnusedBranches();
			
			// make the MTreeSelection to read the TEntryList file
			myTreeSelections = new MTreeSelection(cutFile);
			entry_number = myTreeSelections->GetNextEntry(topCutName);
			ReadEntryNtuple(entry_number);
		} else {
			driver=false;
			myTreeSelections = reinterpret_cast<MTreeSelection*>(myTreeSelectionsPtr);
			// else we're filling the cut while processing
			// get the treereader from the treeselections
			myTreeReader = myTreeSelections->GetTreeReader();
			entry_number = myTreeReader->GetEntryNumber();
			MAX_ENTRIES=-1;   // disable early termination: let the upstream reader finish.
		}
		first_entry=entry_number; // debug print only
	}
	
	// read efficiencies of the different thresholds of old vs new data
	GetEnergyCutEfficiencies();
	
	return true;
}

bool PurewaterLi9Plots::Execute(){
	
	if(valuesFileMode=="read") return true;
	
	++entries_processed;
	if((entries_processed%1000)==0) Log(toolName+" processing entry "+toString(entry_number), v_warning,verbosity);
	
	// Do Li9 analysis
	Log(toolName+" doing analysis",v_debug,verbosity);
	// we must catch any exceptions to ensure our entry number advances
	// otherwise we can end up getting stuck re-processing the same entry endlessly
	// if the tool aborts before 
	try{
		Analyse();
	}
	catch (...){
		Log(toolName+" encountered error doing Analyse!",v_error,verbosity);
	}
	
	// Get next entry
	if(driver){
		// Get the number of the next passing event from the MTreeSelection
		bool failed;
		do {
			try{ entry_number = myTreeSelections->GetNextEntry(topCutName); failed=false; }
			catch (...){ Log(toolName+" failed getting next entry!",v_error,verbosity); failed=true; }
		} while (failed);
	} else {
		// if we're not the event driver (i.e. there is an upstream tool reading the input file)
		// then just retrieve the event number that will be processed next
		entry_number = myTreeReader->GetEntryNumber();
	}
	
	// stop at user-defined limit to the number of events to process
	if((MAX_ENTRIES>0)&&(entries_processed>MAX_ENTRIES)){
		Log(toolName+" reached MAX_ENTRIES, setting StopLoop",v_error,verbosity);
		m_data->vars.Set("StopLoop",1);
	} else {
		// Pre-Load next input entry so we can stop the toolchain if we're
		// about to run off the end of the tree, or if we encounter a read error
		get_ok = ReadEntryNtuple(entry_number);
		if(get_ok==0){
			m_data->vars.Set("StopLoop",1);
			Log(toolName+" Hit end of input file, stopping loop",v_warning,verbosity);
		}
		else if(get_ok==-2){
			Log(toolName+" Error during AutoClear while loading next input ntuple entry!",v_error,verbosity);
			return false;
		}
		else if(get_ok<0){
			Log(toolName+" IO error loading next input ntuple entry!",v_error,verbosity);
			return false;
		}
	}
	
	return true;
}

int PurewaterLi9Plots::ReadEntryNtuple(long entry_number){
	// get next entry from TreeReader
	int bytesread = myTreeReader->GetEntry(entry_number);
	if(bytesread<=0) return bytesread;
	
	// retrieve variables from branches
	// for efficiency, add all used branches to DisableUnusedBranches
	// TODO reduce this to just what is needed for plotting
	int success = 
	(myTreeReader->GetBranchValue("HEADER", HEADER)) &&
	(myTreeReader->GetBranchValue("LOWE", LOWE)) &&
	(myTreeReader->GetBranchValue("ThirdRed", thirdredvars)) &&
	(myTreeReader->GetBranchValue("np", num_neutron_candidates)) &&
	(myTreeReader->GetBranchValue("N200M", max_hits_200ns_AFT)) &&
	(myTreeReader->GetBranchValue("neutron5", ntag_FOM)) &&
	(myTreeReader->GetBranchValue("dt", dt_lowe_n)) &&
	(myTreeReader->GetBranchValue("nmusave_pre", num_pre_muons)) &&
	(myTreeReader->GetBranchValue("nmusave_post", num_post_muons)) &&
	(myTreeReader->GetBranchValue("mubstatus", mu_class)) &&
	(myTreeReader->GetBranchValue("mubitrack", mu_index)) &&
	(myTreeReader->GetBranchValue("mubgood", mu_fit_goodness)) &&
	(myTreeReader->GetBranchValue("spadt", dt_mu_lowe)) &&
	(myTreeReader->GetBranchValue("spadlt", dlt_mu_lowe)) &&
	(myTreeReader->GetBranchValue("multispa_dist", closest_lowe_60s));
	
	return success;
}

int PurewaterLi9Plots::DisableUnusedBranches(){
	// TODO reduce this to just what is needed for plotting
	std::vector<std::string> used_branches{
		// list actively used branches here
		"HEADER",
		"LOWE",
		"ThirdRed",
		"np",
		"N200M",
		"neutron5",
		"dt",
		"nmusave_pre",
		"nmusave_post",
//		"spaloglike",
//		"spaloglike_shfld",
//		"spaloglike_kirk",
//		"spaloglike_kirk_shfld",
//		"spaloglike_li9",
//		"spaloglike_kirk_li9",
		"mubstatus",
//		"mubntrack",
		"mubitrack",
		"mubgood",
//		"mubffgood",
		"spadt",
		"spadlt",
//		"spadll",
//		"spadll_kirk",
//		"sparesq",
//		"spaqpeak",
//		"spaqpeak_kirk",
//		"spamuqismsk",
//		"spamuqismsk_pertrack",
//		"spadts",
//		"spadt_li9",
//		"candidates",
//		"muindex",
//		"neut_flag",
//		"mult_flag",
//		"neut_shift",
//		"neutdiff",
		"multispa_dist"
	};
	
	return myTreeReader->OnlyEnableBranches(used_branches);
}

bool PurewaterLi9Plots::CreateOutputFile(std::string outputFile){
	Log(toolName+": Creating output file "+outputFile,v_message,verbosity);
	outfile = new TFile(outputFile.c_str(), "RECREATE");
	outfile->cd();
	return true;
}

void PurewaterLi9Plots::CloseFile(){
	outfile->Write("*",TObject::kOverwrite);
	Log(toolName+" closing output file",v_debug,verbosity);
	outfile->Close();
	delete outfile;
	outfile = nullptr;
}

// =====================================================================
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bool PurewaterLi9Plots::Analyse(){
	// This gets called for each Execute iteration, to process one lowe event
	
	// we need to check the last cut before any histogram fill calls
	if(not myTreeSelections->GetPassesCut("lowe_energy>6MeV")) return true;
	
	// plot lt, dt distributions as a function of muon class, for both pre- and post-muons
	Log(toolName+" looping over "+toString(num_pre_muons)+" preceding muons and "+toString(num_post_muons)
				+" following muons to fill spallation and control dl, dt histograms",v_debug,verbosity);
	
	// pre muons
	// only consider first muboy muon (only for multi-mu events?)
	std::set<size_t> pre_muboy_first_muons = myTreeSelections->GetPassingIndexes("pre_muon_muboy_i==0");
	for(int mu_i : pre_muboy_first_muons){
		Log(toolName+" filling spallation dt and dlt distributions",v_debug+2,verbosity);
		dlt_vals_pre.at(mu_class[mu_i]).push_back(dlt_mu_lowe[mu_i]);   // FIXME weight by num_pre_muons
		dt_vals_pre.at(mu_class[mu_i]).push_back(dt_mu_lowe[mu_i]);     // FIXME weight by num_pre_muons
		
		// to evaluate systematic on lt cut, apply various dt cuts and see how the lt cut efficiency varies
		// since we're interested in the effect on the spallation sample, which is given by
		// the total - post-muon sample, record both pre- and post- muon samples with various dt cuts
		for(int dt_cut_i=0; dt_cut_i<num_dt_cuts; ++dt_cut_i){
			Log(toolName+" checking nominal dlt cut systematic",v_debug+2,verbosity);
			if(myTreeSelections->GetPassesCut("pre_mu_dt_cut_"+toString(dt_cut_i),mu_i)){
				Log(toolName+" filling spallation dlt distribution for dt cut "
				            +toString(dt_cut_i),v_debug+2,verbosity);
				dlt_systematic_dt_cuts_pre.at(dt_cut_i).push_back(dt_mu_lowe[mu_i]);
			}
		}
	}
	// post muons
	std::set<size_t> post_muboy_first_muons = myTreeSelections->GetPassingIndexes("post_muon_muboy_i==0");
	for(int mu_i : post_muboy_first_muons){
		Log(toolName+" filling spallation dt and dlt distributions",v_debug+2,verbosity);
		dlt_vals_post.at(mu_class[mu_i]).push_back(dlt_mu_lowe[mu_i]);   // FIXME weight by num_post_muons
		dt_vals_post.at(mu_class[mu_i]).push_back(dt_mu_lowe[mu_i]);     // FIXME weight by num_post_muons
		
		for(int dt_cut_i=0; dt_cut_i<num_dt_cuts; ++dt_cut_i){
			Log(toolName+" checking nominal dlt cut systematic",v_debug+2,verbosity);
			if(myTreeSelections->GetPassesCut("post_mu_dt_cut_"+toString(dt_cut_i),mu_i)){
				Log(toolName+" filling spallation dlt distribution for dt cut "
				            +toString(dt_cut_i),v_debug+2,verbosity);
				dlt_systematic_dt_cuts_post.at(dt_cut_i).push_back(dt_mu_lowe[mu_i]);
			}
		}
	}
	// in Finalise we'll substract the two to get dt and dlt distributions for spallation only.
	// we'll also compare across various dt cuts to get the systematic error on the spallation dlt cut.
	
	// the following cuts are based on muon-lowe pair variables, so loop over muon-lowe pairs
	std::set<size_t> spall_mu_indices = myTreeSelections->GetPassingIndexes("dlt_mu_lowe>200cm");
	Log(toolName+" Looping over "+toString(spall_mu_indices.size())
				+" preceding muons to look for spallation events",v_debug,verbosity);
	for(size_t mu_i : spall_mu_indices){
		// record the distribution of dt_mu_lowe
		Log(toolName+" filling mu->lowe dt distribution for spallation entries",v_debug+2,verbosity);
		// this *should* only contain pre muons with dt < 0:
		if(dt_mu_lowe[mu_i]>0){
			Log(toolName+" error! Spallation muon with time "+toString(dt_mu_lowe[mu_i])+" after lowe event!",
				v_warning,verbosity);
			continue;
		}
		dt_mu_lowe_vals.push_back(dt_mu_lowe[mu_i]);  // FIXME weight by num_pre_muons
		// in finalise we'll fit this distribution to estimate the number of events from each isotope
		
		// now check whether this passed the additional Li9 cuts
		if(not myTreeSelections->GetPassesCut("ntag_FOM>0.995")) continue;
		
		// plot distribution of beta energies from passing triplets, compare to fig 4
		Log(toolName+" filling li9 candidate distributions",v_debug+2,verbosity);
		li9_e_vals.push_back(LOWE->bsenergy); // FIXME weight by num_post_muons
		
		// plot distirbution of mu->beta   dt from passing triplets, compare to fig 6
		li9_muon_dt_vals.push_back(fabs(dt_mu_lowe[mu_i])); // FIXME weight by num_post_muons
		
		// plot distribution of beta->ntag dt from passing triplets, compare to fig 5
		// Zhang had no events with >1 ntag candidate: should we only take the first? XXX
		for(size_t neutron_i=0; neutron_i<num_neutron_candidates; ++neutron_i){
			if(myTreeSelections->GetPassesCut("mu_lowe_ntag_triplets",{mu_i, neutron_i})){
				double ncap_time = dt_lowe_n[neutron_i];  // these times are in nanoseconds
				// according to sonias code we need to account for some offset of the AFT trigger timestamps?
				// "Remove 10mus time shift for AFT events and convert time to microseconds"
				// doesn't seem to tie up with what this is actually doing, though
				// adjusted too instead convert presumably ms, to seconds for consistency
				double ncap_time_adjusted = ncap_time < 50000 ? ncap_time : ncap_time - 65000;
				li9_ntag_dt_vals.push_back(ncap_time_adjusted/1E9);
			}
		}
		
	} // end loop over muons
	
	return true;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// =====================================================================

bool PurewaterLi9Plots::Finalise(){
	// if we processed the data, print some stats
	if(myTreeSelections){
		Log(toolName+" processed "+toString(entries_processed)+" from entry "+toString(first_entry)+" to "
			+toString(entry_number),v_warning,verbosity);
		Log(toolName+" event counts trace: ",v_warning,verbosity);
		myTreeSelections->PrintCuts();
	}
	
	// if saving or retrieving values from a BoostStore, do that now
	if(valuesFileMode=="write"){
		// set all the values into the BoostStore
		valueStore = new BoostStore(true,constants::BOOST_STORE_BINARY_FORMAT);
		valueStore->Set("livetime",livetime);
		valueStore->Set("dlt_vals_pre",dlt_vals_pre);
		valueStore->Set("dt_vals_pre",dt_vals_pre);
		valueStore->Set("dlt_vals_post",dlt_vals_post);
		valueStore->Set("dt_vals_post",dt_vals_post);
		valueStore->Set("dlt_systematic_dt_cuts_pre",dlt_systematic_dt_cuts_pre);
		valueStore->Set("dlt_systematic_dt_cuts_post",dlt_systematic_dt_cuts_post);
		valueStore->Set("dt_mu_lowe_vals",dt_mu_lowe_vals);
		valueStore->Set("li9_muon_dt_vals",li9_muon_dt_vals);
		valueStore->Set("li9_ntag_dt_vals",li9_ntag_dt_vals);
		valueStore->Set("li9_e_vals",li9_e_vals);
		std::cout<<"out containers have sizes: "<<std::endl;
		std::cout<<"dlt_vals_pre: {";
		for(auto&& av : dlt_vals_pre){ std::cout<<av.size()<<", "; } std::cout<<"\b\b}\n";
		std::cout<<"dt_vals_pre: {";
		for(auto&& av : dt_vals_pre){ std::cout<<av.size()<<", "; } std::cout<<"\b\b}\n";
		std::cout<<"dlt_vals_post: {";
		for(auto&& av : dlt_vals_post){ std::cout<<av.size()<<", "; } std::cout<<"\b\b}\n";
		std::cout<<"dt_vals_post: {";
		for(auto&& av : dt_vals_post){ std::cout<<av.size()<<", "; } std::cout<<"\b\b}\n";
		std::cout<<"dlt_systematic_dt_cuts_pre: {";
		for(auto&& av : dlt_systematic_dt_cuts_pre){ std::cout<<av.size()<<", "; } std::cout<<"\b\b}\n";
		std::cout<<"dlt_systematic_dt_cuts_post: {";
		for(auto&& av : dlt_systematic_dt_cuts_post){ std::cout<<av.size()<<", "; } std::cout<<"\b\b}\n";
		std::cout<<"dt_mu_lowe_vals = "<<dt_mu_lowe_vals.size()<<std::endl;
		std::cout<<"li9_muon_dt_vals = "<<li9_muon_dt_vals.size()<<std::endl;
		std::cout<<"li9_ntag_dt_vals = "<<li9_ntag_dt_vals.size()<<std::endl;
		std::cout<<"li9_e_vals = "<<li9_e_vals.size()<<std::endl;
		
		// save BoostStore
		valueStore->Save(valuesFile.c_str());
		valueStore->Close(); // necessary to complete the file write!
	} else if(valuesFileMode=="read"){
		valueStore = new BoostStore(true,constants::BOOST_STORE_BINARY_FORMAT);
		valueStore->Initialise(valuesFile.c_str());
		valueStore->Get("livetime",livetime);
		valueStore->Get("dlt_vals_pre",dlt_vals_pre);
		valueStore->Get("dt_vals_pre",dt_vals_pre);
		valueStore->Get("dlt_vals_post",dlt_vals_post);
		valueStore->Get("dt_vals_post",dt_vals_post);
		valueStore->Get("dlt_systematic_dt_cuts_pre",dlt_systematic_dt_cuts_pre);
		valueStore->Get("dlt_systematic_dt_cuts_post",dlt_systematic_dt_cuts_post);
		valueStore->Get("dt_mu_lowe_vals",dt_mu_lowe_vals);
		valueStore->Get("li9_muon_dt_vals",li9_muon_dt_vals);
		valueStore->Get("li9_ntag_dt_vals",li9_ntag_dt_vals);
		valueStore->Get("li9_e_vals",li9_e_vals);
		std::cout<<"in containers have sizes: "<<std::endl;
		std::cout<<"dlt_vals_pre: {";
		for(auto&& av : dlt_vals_pre){ std::cout<<av.size()<<", "; } std::cout<<"\b\b}\n";
		std::cout<<"dt_vals_pre: {";
		for(auto&& av : dt_vals_pre){ std::cout<<av.size()<<", "; } std::cout<<"\b\b}\n";
		std::cout<<"dlt_vals_post: {";
		for(auto&& av : dlt_vals_post){ std::cout<<av.size()<<", "; } std::cout<<"\b\b}\n";
		std::cout<<"dt_vals_post: {";
		for(auto&& av : dt_vals_post){ std::cout<<av.size()<<", "; } std::cout<<"\b\b}\n";
		std::cout<<"dlt_systematic_dt_cuts_pre: {";
		for(auto&& av : dlt_systematic_dt_cuts_pre){ std::cout<<av.size()<<", "; } std::cout<<"\b\b}\n";
		std::cout<<"dlt_systematic_dt_cuts_post: {";
		for(auto&& av : dlt_systematic_dt_cuts_post){ std::cout<<av.size()<<", "; } std::cout<<"\b\b}\n";
		std::cout<<"dt_mu_lowe_vals = "<<dt_mu_lowe_vals.size()<<std::endl;
		std::cout<<"li9_muon_dt_vals = "<<li9_muon_dt_vals.size()<<std::endl;
		std::cout<<"li9_ntag_dt_vals = "<<li9_ntag_dt_vals.size()<<std::endl;
		std::cout<<"li9_e_vals = "<<li9_e_vals.size()<<std::endl;
	}
	
	// create the output file for plots and fits
	// -----------------------
	CreateOutputFile(outputFile);
	
	// prettify plots
	gStyle->SetOptTitle(0); gStyle->SetOptStat(0);
	gStyle->SetOptFit(1111); gStyle->SetStatBorderSize(0);
	gStyle->SetStatX(.89); gStyle->SetStatY(.89);
	
	// subtract distribution of post-muons from pre-muons to obtain lt and dt distributions of mu-lowe pairs
	Log(toolName+" making plots of muon-lowe Dt and Dlt distributions",v_debug,verbosity);
	PlotMuonDt();
	PlotMuonDlt();
	
	// this is for candidates that pass the dlt cut: compare to paper fig 3
	Log(toolName+" making plots of spallation Dt distributions",v_debug,verbosity);
	PlotSpallationDt();
	
	// mu-lowe time diff for li9 triplet candidates: compare to paper Fig 6
	Log(toolName+" making plot of Li9 candidates Dt distribution",v_debug,verbosity);
	PlotLi9LifetimeDt();
	
	Log(toolName+" making plots of Li9 candidate beta spectrum",v_debug,verbosity);
	PlotLi9BetaEnergy();
	
	Log(toolName+" making plot of Li9 ncapture Dt distribution",v_debug,verbosity);
	PlotNcaptureDt();
	
	// measure dlt cut systematic TODO
//	MeasureDltSystematic(); // using variation in dlt_systematic_dt_cuts in bin corresponding to dlt=200cm
	// for each entry in dlt_systematic_dt_cuts_pre, dlt_systematic_dt_cuts_post, subtract post from pre.
	// this gives a set of efficiencies of spallation for varying dts.
	// compare efficiency across various dts in each dlt bin to get the dlt systematic error.
	
	// fit time distribution of muon-lowe pairs with the nominal dlt cut to obtain isotope counts
//	MeasureIsotopeRates();  // TODO
	// split into sum of contributions from each isotope by doing piecewise fit to time distribution:
	// F(t) = sum_i (N_i/τ_i) exp(-dt/τ_i) + const
	verbosity=0;
	
	Log(toolName+" writing outputs and closing file",v_debug,verbosity);
	CloseFile();
	
	Log(toolName+" closing input file",v_debug,verbosity);
	if(myTreeReader) delete myTreeReader;
	myTreeReader=nullptr;
	
	Log(toolName+" closing cut file",v_debug,verbosity);
	if(myTreeSelections) delete myTreeSelections;
	myTreeSelections=nullptr;
	
	Log(toolName+" finished",v_debug,verbosity);
	
	return true;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// =====================================================================

bool PurewaterLi9Plots::PlotMuonDlt(){
	
	THStack my_spall_dlts;
	// make histograms of transverse distance to muon for all pre- and post-muons
	// and their difference to extract the spallation distributions
	for(auto&& aclass : constants::muboy_class_to_name){  // we have 5 muboy classifications
		int mu_class_i = aclass.first;
		const char* mu_class_name = aclass.second.c_str();
		TH1F ahist_pre(TString::Format("dlt_pre_%d",mu_class_i),
					     "All Pre-Muon to Low-E Transverse Distances",8,0,400);
		for(auto&& aval : dlt_vals_pre.at(mu_class_i)) ahist_pre.Fill(aval);
		ahist_pre.Write();
		
		TH1F ahist_post(TString::Format("dlt_post_%d",mu_class_i),
					     "All Post-Muon to Low-E Transverse Distances",8,0,400);
		for(auto&& aval : dlt_vals_post.at(mu_class_i)) ahist_post.Fill(aval);
		ahist_post.Write();
		
		TH1F* dlt_hist_spall = (TH1F*)ahist_pre.Clone(TString::Format("dlt_spall_%s",mu_class_name));
		dlt_hist_spall->Reset(); // clear entries
		// subtract post from pre
		dlt_hist_spall->Add(&ahist_pre,&ahist_post,1,-1);
		std::cout<<"subtracting "<<ahist_post.GetEntries()<<" post entries from "<<ahist_pre.GetEntries()
				 <<" pre entries results in "<<dlt_hist_spall->GetEntries()<<" ("
				 <<(ahist_pre.GetEntries()-ahist_post.GetEntries())<<") entries"<<std::endl;
		dlt_hist_spall->Scale(1./dlt_hist_spall->Integral());
		// XXX we seem to lose half our entries here....??
		dlt_hist_spall->Write();
		dlt_hist_spall->SetDirectory(0);
		my_spall_dlts.Add(dlt_hist_spall);
	}
	
	// create and add the paper versions
	PlotPaperDlt(my_spall_dlts);
	
	my_spall_dlts.Write("spall_dls_stack");
	
	// cleanup
	my_spall_dlts.GetHists()->Delete();
	
	return true;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// =====================================================================

bool PurewaterLi9Plots::PlotMuonDt(){
	// make histograms of time difference to muon for all pre- and post-muons
	// and their difference to extract spallation distributions.
	
	// colour match to the existing paper
	// muboy_classes{ misfit=0, single_thru_going=1, single_stopping=2, multiple_mu=3, also_multiple_mu=4, corner_clipper=5};
	std::map<std::string, EColor> class_colours{
		{"misfit",kBlue},
		{"single_thru_going",kBlack},
		{"single_stopping",kMagenta},
		{"multiple_mu",kRed},
		{"also_multiple_mu",kRed},
		{"corner_clipper",kWhite}
	}; // corner clippers not shown...
	
	// plot over two ranges: full range (30s) and paper plotted range (0.25s)
	std::vector<float> dt_range_full{30,0.25};
	THStack my_spall_dts;
	for(auto&& dt_max : dt_range_full){
		for(auto&& aclass : constants::muboy_class_to_name){  // we have 5 muboy classifications
			int mu_class_i = aclass.first;
			const char* mu_class_name = aclass.second.c_str();
			// need to take the fabs of the time so time 0 is in bin 0 for both pre- and post-
			// in order to be able to subtract the bin counts.
			TH1F ahist_pre(TString::Format("dt_pre_%d_%0.2f",mu_class_i,dt_max),
							 "All Pre-Muon to Low-E Time Differences",10,0,dt_max);
			for(auto&& aval : dt_vals_pre.at(mu_class_i)) ahist_pre.Fill(fabs(aval));
			ahist_pre.Write();
			
			TH1F ahist_post(TString::Format("dt_post_%d_%0.2f",mu_class_i,dt_max),
							 "All Post Muon to Low-E Time Differences",10,0,dt_max);
			for(auto&& aval : dt_vals_post.at(mu_class_i)) ahist_post.Fill(aval);
			ahist_post.Write();
			
			TH1F* dt_hist_spall = 
				(TH1F*)ahist_pre.Clone(TString::Format("dt_spall_%s_%0.2f",mu_class_name,dt_max));
			dt_hist_spall->Reset(); // clear entries
			// subtract post from pre
			dt_hist_spall->Add(&ahist_pre,&ahist_post,1,-1);
			dt_hist_spall->SetLineColor(class_colours.at(mu_class_name));
			dt_hist_spall->Scale(1./dt_hist_spall->Integral());
			dt_hist_spall->Write();
			dt_hist_spall->SetDirectory(0);
			if(dt_max<20.f) my_spall_dts.Add(dt_hist_spall);  // add shorter axis ones for comparison
		}
	}
	
	// create and add the paper ones for overlay
	PlotPaperDt(my_spall_dts);
	
	// write the stack to file so we can easily see all histos
	my_spall_dts.Write("spall_dts_stack");
	
	// TH1::Clone creates a copy that we own, so are responsible for cleanup.
	// THStacks do not take ownership of their histos, and we can't make them.
	// The way to cleanup is therefore:
	my_spall_dts.GetHists()->Delete();
	
	return true;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// =====================================================================

bool PurewaterLi9Plots::PlotSpallationDt(){
	// dt distribution of muons passing dlt<200cm cut
	TH1F dt_mu_lowe_hist("dt_mu_lowe_hist","Spallation Muon to Low-E Time Differences",5000,0,30);
	TH1F dt_mu_lowe_hist_short("dt_mu_lowe_hist_short","Spallation Muon to Low-E Time Differences",500,0,0.25);
	for(auto&& aval : dt_mu_lowe_vals){
		if(aval>0) continue;  // only fit preceding muons (true spallation)
		dt_mu_lowe_hist.Fill(fabs(aval));
		dt_mu_lowe_hist_short.Fill(fabs(aval));
	}
	dt_mu_lowe_hist.Write();
	dt_mu_lowe_hist_short.Write();
	
	// make a histo with equally spaced bins in log scale
	// we need to fit the same histogram (with the same binning) for the intermediate fits,
	// otherwise the bin widths change, the contents change, and the fit parameters aren't the same.
	// 5000 bins gives chi2/NDOF that matches the paper.
	int nbins=5000;
	std::vector<double> binedges = MakeLogBins(0.001, 30, nbins+1);
	TH1F dt_mu_lowe_hist_log("dt_mu_lowe_hist_log","Data;dt(s);Events/0.006 s",
							 nbins, binedges.data());
	for(auto&& aval : dt_mu_lowe_vals){ if(aval>0) continue; dt_mu_lowe_hist_log.Fill(fabs(aval)); }
	// scale each bin by its width to correct back to number of events per equal time interval
	for(int bini=1; bini<dt_mu_lowe_hist_log.GetNbinsX()+1; ++bini){
		// numbers here are from binning of dt_mu_lowe_hist
		dt_mu_lowe_hist_log.SetBinContent(bini,
			dt_mu_lowe_hist_log.GetBinContent(bini)/(dt_mu_lowe_hist_log.GetBinWidth(bini)/(30./5000.)));
	}
//	std::cout<<"this plot will be used for all the time distribution fits:"<<std::endl;
//	dt_mu_lowe_hist_log.Draw();
//	gPad->WaitPrimitive();
	
	// fit this with all the rates....
	// we do the fitting in 5 stages, initially fitting sub-ranges of the distribution
	for(int i=0; i<5; ++i) FitSpallationDt(dt_mu_lowe_hist, dt_mu_lowe_hist_log, i);
	
	// the yield across the whole energy range is obtained from the fit amplitude via:
	// Yi = Ni / (Rmu * T * rho * Lmu)
	// where Rmu is the muon rate, T the live time, rho the density of water and Lmu the
	// measured path length of the muon track...? is this event-wise? Yield... maybe.
	
	// the production rate integrated over the whole energy range is given by:
	// Ri = Ni / (FV * T * eff_i)
	// where FV is the fid vol and T again live time.
	// note this is using the paper definiton, not zhangs definition. just changes defn of Ni<->Ni*eff_i
	std::map<std::string,double> rates;
	for(auto&& anisotope : fit_amps){
		std::string isotope = anisotope.first;
		if(isotope.substr(0,5)=="const") continue;
		// we need to know the efficiency of the selection, which is obtained
		// from the fraction of beta spectrum that's below 6MeV
		std::cout<<"getting efficiency of "<<isotope<<std::endl;
		if(isotope.find("_")==std::string::npos){
			double efficiency = efficiencies.at(anisotope.first);
			rates[anisotope.first] = anisotope.second/efficiency * (fiducial_vol / livetime);
			std::cout<<"efficiency of "<<anisotope.first<<", Rate of production is "<<rates[anisotope.first]
					 <<"/kton/day"<<std::endl;
		} else {
			std::string first_isotope = isotope.substr(0,isotope.find_first_of("_"));
			std::string second_isotope = isotope.substr(isotope.find_first_of("_")+1,std::string::npos);
			double first_efficiency = efficiencies.at(first_isotope);
			double second_efficiency = efficiencies.at(second_isotope);
			rates[first_isotope] = anisotope.second/first_efficiency * (fiducial_vol / livetime);
			rates[second_isotope] = anisotope.second/second_efficiency * (fiducial_vol / livetime);
			std::cout<<"efficiency of "<<first_isotope<<" is "<<first_efficiency
					 <<", calculated rate of production is "<<rates[first_isotope]
					 <<"/kton/day"<<std::endl
					 <<", efficiency of "<<second_isotope<<" is "<<second_efficiency
					 <<", calculated rate of production is "<<rates[second_isotope]
					 <<"/kton/day"<<std::endl;
		}
	}
	// TODO save this map to an ouput file
	// not worth doing till we re-run and determine the livetime
	return true;
}

std::vector<double> PurewaterLi9Plots::MakeLogBins(double xmin, double xmax, int nbins){
	std::vector<double> binedges(nbins+1);
	double xxmin=log10(xmin);
	double xxmax = log10(xmax);
	for(int i=0; i<nbins; ++i){
		binedges[i] = pow(10,xxmin + (double(i)/(nbins-1.))*(xxmax-xxmin));
	}
	binedges[nbins+1] = xmax; // required
	return binedges;
}

bool PurewaterLi9Plots::FitSpallationDt(TH1F& dt_mu_lowe_hist, TH1F& dt_mu_lowe_hist_log, int rangenum){
	// do sub-range fitting based on section B of the 2015 paper
	// formula 1 from the paper, with a sub-range as described in section B2
	// we have 4 sub-ranges to fit
	switch (rangenum){
	case 0:{
		// fit range 50us -> 0.1s with 12B + 12N only
//		TH1F hist_to_fit("hist_to_fit","Spallation Muon to Low-E Time Differences",100,50e-6,0.1);
//		for(auto&& aval : dt_mu_lowe_vals){ if(aval>0) continue; hist_to_fit.Fill(fabs(aval)); }
		// make the function which we'll fit
		std::cout<<"calling BuildFunction for case "<<rangenum<<", 12B+12N"<<std::endl;
		TF1 func_sum = BuildFunction({"12B","12N"},50e-6,0.1);
		// do the fit
		std::cout<<"fitting"<<std::endl;
//		hist_to_fit.Fit(&func_sum,"","",50e-6,0.1);
		dt_mu_lowe_hist_log.Fit(&func_sum,"R","",50e-6,0.1);
		// record the results for the next step
		std::cout<<"recording fit results"<<std::endl;
		PushFitAmp(func_sum,"12B");
		PushFitAmp(func_sum,"12N");
		fit_amps["const_0"] = func_sum.GetParameter("const");
		std::cout<<"recorded results were: "<<std::endl
				 <<"12B: "<<fit_amps["12B"]<<std::endl
				 <<"12N: "<<fit_amps["12N"]<<std::endl
				 <<"const_0: "<<fit_amps["const_0"]<<std::endl;
//		hist_to_fit.Draw();
		dt_mu_lowe_hist_log.Draw();
		//gPad->WaitPrimitive();
		dt_mu_lowe_hist_log.GetListOfFunctions()->Clear();
		break;
		}
	case 1:{
		// fit range 6-30s with 16N + 11B only
//		TH1F hist_to_fit("hist_to_fit","Spallation Muon to Low-E Time Differences",100,6,30);
//		for(auto&& aval : dt_mu_lowe_vals){ if(aval>0) continue; hist_to_fit.Fill(fabs(aval)); }
		std::cout<<"calling BuildFunction for case "<<rangenum<<", 16N+11Be"<<std::endl;
		TF1 func_sum = BuildFunction({"16N","11Be"},6,30);
		std::cout<<"fitting"<<std::endl;
//		hist_to_fit.Fit(&func_sum,"","",6,30);
		dt_mu_lowe_hist_log.Fit(&func_sum,"R","",6,30);
		// record the results for the next step
		std::cout<<"recording the fit results"<<std::endl;
		PushFitAmp(func_sum,"16N");
		PushFitAmp(func_sum,"11Be");
		fit_amps["const_1"] = func_sum.GetParameter("const");
		std::cout<<"recorded results were: "<<std::endl
				 <<"16N: "<<fit_amps["16N"]<<std::endl
				 <<"11Be: "<<fit_amps["11Be"]<<std::endl
				 <<"const_1: "<<fit_amps["const_1"]<<std::endl;
//		hist_to_fit.Draw();
		dt_mu_lowe_hist_log.Draw();
		//gPad->WaitPrimitive();
		dt_mu_lowe_hist_log.GetListOfFunctions()->Clear();
		break;
		}
	case 2:{
		// fit the range 0.1-0.8s with the components previously fit now fixed,
//		TH1F hist_to_fit("hist_to_fit","Spallation Muon to Low-E Time Differences",100,0.1,0.8);
//		for(auto&& aval : dt_mu_lowe_vals){ if(aval>0) continue; hist_to_fit.Fill(fabs(aval)); }
		// allowing additional components Li9 + a combination of 8He+9C
		std::cout<<"calling BuildFunction for case "<<rangenum<<", 9Li+8He_9C+8Li_8B+12B+12N+16N+11Be"<<std::endl;
		TF1 func_sum = BuildFunction({"9Li","8He_9C","8Li_8B","12B","12N","16N","11Be"},0.1,0.8);
		std::cout<<"retrieving past results"<<std::endl;
		// pull fit results from the last two stages
		PullFitAmp(func_sum,"12B");
		PullFitAmp(func_sum,"12N");
		PullFitAmp(func_sum,"16N");
		PullFitAmp(func_sum,"11Be");
		// fit the new components
		std::cout<<"fitting"<<std::endl;
//		hist_to_fit.Fit(&func_sum,"","",0.1,0.8);
		dt_mu_lowe_hist_log.Fit(&func_sum,"R","",0.1,0.8);
		// record the results for the next step
		std::cout<<"recording the fit results"<<std::endl;
		PushFitAmp(func_sum,"9Li");
		PushFitAmp(func_sum,"8He_9C");
		PushFitAmp(func_sum,"8Li_8B");
		fit_amps["const_2"] = func_sum.GetParameter("const");
		std::cout<<"recorded results were: "<<std::endl
				 <<"9Li: "<<fit_amps["9Li"]<<std::endl
				 <<"8He_9C: "<<fit_amps["8He_9C"]<<std::endl
				 <<"8Li_8B: "<<fit_amps["8Li_8B"]<<std::endl
				 <<"const_2: "<<fit_amps["const_2"]<<std::endl;
//		hist_to_fit.Draw();
		dt_mu_lowe_hist_log.Draw();
		//gPad->WaitPrimitive();
		dt_mu_lowe_hist_log.GetListOfFunctions()->Clear();
		
		// debug check, let's see what's going on here
		std::cout<<"Drawing past fits and this one, to see how they look"<<std::endl;
		TF1 func_early = BuildFunction({"12B","12N"},50E-6,30);
		// pull fit results from the last two stages
		PullFitAmp(func_early,"12B");
		PullFitAmp(func_early,"12N");
		PullFitAmp(func_early,"const_0");
		func_early.SetLineColor(kBlue);
		std::cout<<"early"<<std::endl;
		func_early.Draw();
		//gPad->WaitPrimitive();
		TF1 func_late = BuildFunction({"16N","11Be"},50E-6,30);
		// pull fit results from the last two stages
		PullFitAmp(func_late,"16N");
		PullFitAmp(func_late,"11Be");
		PullFitAmp(func_late,"const_1");
		func_late.SetLineColor(kGreen-1);
		std::cout<<"late"<<std::endl;
		func_late.Draw();
		//gPad->WaitPrimitive();
		
		func_sum.SetRange(50E-6,30);
		std::cout<<"int"<<std::endl;
		func_sum.Draw();
		//gPad->WaitPrimitive();
		
		std::cout<<"everything"<<std::endl;
		dt_mu_lowe_hist_log.Draw();
		func_sum.Draw("same");
		func_early.Draw("same");
		func_late.Draw("same");
		//gPad->WaitPrimitive();
		
		break;
		}
	case 3:{
		// fit the range 0.8-6s with the components previously fit now fixed,
//		TH1F hist_to_fit("hist_to_fit","Spallation Muon to Low-E Time Differences",100,0.8,6);
//		for(auto&& aval : dt_mu_lowe_vals){ if(aval>0) continue; hist_to_fit.Fill(fabs(aval)); }
		// allowing additional components 15C + 16N
		std::cout<<"calling BuildFunction for case "<<rangenum<<", 15C+16N+8Li_8B"<<std::endl;
		TF1 func_sum = BuildFunction({"15C","16N","8Li_8B"},0.8,6);
		std::cout<<"retrieving past results"<<std::endl;
		PullFitAmp(func_sum,"8Li_8B");
		// fit the new components
		std::cout<<"fitting"<<std::endl;
//		hist_to_fit.Fit(&func_sum,"","",0.8,6);
		dt_mu_lowe_hist_log.Fit(&func_sum,"R","",0.8,6);
		// record the results for the next step
		std::cout<<"recording fit results"<<std::endl;
		PushFitAmp(func_sum,"15C");
		PushFitAmp(func_sum,"16N");
		fit_amps["const_3"] = func_sum.GetParameter("const");
//		hist_to_fit.Draw();
		dt_mu_lowe_hist_log.Draw();
		//gPad->WaitPrimitive();
		dt_mu_lowe_hist_log.GetListOfFunctions()->Clear();
		break;
		}
	case 4:{
		// final case: release all the fixes but keep the previously fit values as starting points.
		std::cout<<"calling BuildFunction for case "<<rangenum<<", everything"<<std::endl;
		TF1 func_sum = BuildFunction({"12B","12N","16N","11Be","9Li","8He_9C","8Li_8B","15C"},0,30);
		func_sum.SetLineColor(kRed);
		std::cout<<"retrieving past results"<<std::endl;
		PullFitAmp(func_sum,"12B",false);
		PullFitAmp(func_sum,"12N",false);
		PullFitAmp(func_sum,"16N",false);
		PullFitAmp(func_sum,"11Be",false);
		PullFitAmp(func_sum,"9Li",false);
		PullFitAmp(func_sum,"8He_9C",false);
		PullFitAmp(func_sum,"8Li_8B",false);
		PullFitAmp(func_sum,"15C",false);
		
		/*
		// make a histo with equally spaced bins in log scale
		std::vector<double> binedges = MakeLogBins(0.001, 30, 100+1);
		TH1F hist_to_fit("hist_to_fit","Spallation Muon to Low-E Time Differences",100, binedges.data());
		for(auto&& aval : dt_mu_lowe_vals){ if(aval>0) continue; hist_to_fit.Fill(fabs(aval)); }
		// scale each bin by its width to correct back to number of events per equal time interval
		for(int bini=1; bini<hist_to_fit.GetNbinsX()+1; ++bini){
			// numbers here are from binning of dt_mu_lowe_hist
			hist_to_fit.SetBinContent(bini,hist_to_fit.GetBinContent(bini)/(hist_to_fit.GetBinWidth(bini)/(30./5000.)));
		}
		*/
		std::cout<<"fitting"<<std::endl;
//		dt_mu_lowe_hist.Fit(&func_sum,"","",0,30);   // do not fit to the logarithmic binning right away
		TFitResultPtr fitresult =  dt_mu_lowe_hist_log.Fit(&func_sum,"RS","",0,30);
//		std::cout<<"drawing initial fit result"<<std::endl;
//		hist_to_fit.SetLineColor(kRed);
//		hist_to_fit.Draw();
//		func_sum.Draw("same");
//		dt_mu_lowe_hist.Draw("same");
//		gPad->SetLogx();
//		gPad->SetLogy();
//		gPad->WaitPrimitive();
		
		std::cout<<"fixing pars"<<std::endl;
		// we want to constrain the number of each isotope to be >0, which would
		// mean putting a lower limit of 0 on the amplitude in the fit.
		// unfortunately ROOT only lets us put both upper and lower limits (not just one)
		// and Minuit doesn't like limits of 0 and 10E10 (all sorts of errors).
		// So, we put no constraint on the parameter, but our function (from BuildFunction)
		// only uses its magnitude. This means we end up with fit values that are negative,
		// but resulting function is the same even if we fix the sign to positive. So.
		// let's try to fix up this hack by setting all abundances positive as they should be
		for(int pari=0; pari<func_sum.GetNpar(); ++pari){
			std::string parname = func_sum.GetParName(pari);
			if(parname.substr(0,9)=="lifetime_") continue; // skip lifetimes, they're good
			
//			// when adjusting fits interactively they're given a limit,
//			// but it seems like they don't have one in code, so skip this
//			double parmin, parmax;
//			func_sum.GetParLimits(pari,parmin,parmax);
//			double newupper = std::max(abs(parmin),abs(parmax));
//			double parval = func_sum.GetParameter(pari);
//			func_sum.SetParLimits(pari,parmin,newupper);
//			func_sum.SetParameter(pari,abs(parval));
//			func_sum.SetParLimits(pari,0,newupper);
			
			// much simpler in that case
			func_sum.SetParameter(pari,abs(func_sum.GetParameter(pari)));
		}
		
		/*
		// redo the fit but this time to the logarithmic binned histogram
		TFitResultPtr fitresult =  hist_to_fit.Fit(&func_sum,"S","",0,30); // option S, give us a fitresultptr
//		// seems like this is necessary again - unnecessary as we're using the TFitResultPtr
//		for(int pari=0; pari<func_sum.GetNpar(); ++pari){
//			std::string parname = func_sum.GetParName(pari);
//			if(parname.substr(0,9)=="lifetime_") continue; // skip lifetimes, they're good
//			func_sum.SetParameter(pari,abs(func_sum.GetParameter(pari)));
//		}
		dt_mu_lowe_hist.GetListOfFunctions()->Clear(); // clear the old bad fit
		hist_to_fit.Draw();
		func_sum.Draw("same");
		dt_mu_lowe_hist.Draw("same");
		*/
		dt_mu_lowe_hist_log.Draw();
		gPad->SetLogx();
		gPad->SetLogy();
		
		// since we're having issues, let's see how well the old results fit our data
		TF1 func_paper = BuildFunction({"12B","12N","16N","11Be","9Li","8He_9C","8Li_8B","15C"},0,30);
		// scale the paper down to be around the same num events as us, to ease comparison
		std::cout<<"retrieving paper results"<<std::endl;
		PullPaperAmp(func_paper,"12B",true,paper_scaling);
		PullPaperAmp(func_paper,"12N",true,paper_scaling);
		PullPaperAmp(func_paper,"16N",true,paper_scaling);
		PullPaperAmp(func_paper,"11Be",true,paper_scaling);
		PullPaperAmp(func_paper,"9Li",true,paper_scaling);
		PullPaperAmp(func_paper,"8He_9C",true,paper_scaling);
		//PullPaperAmp(func_paper,"8Li_8B",true,paper_scaling);
		PullPaperAmp(func_paper,"8Li",true,paper_scaling);
		PullPaperAmp(func_paper,"8B",true,paper_scaling);
		PullPaperAmp(func_paper,"15C",true,paper_scaling);
		PullPaperAmp(func_paper,"const",true,paper_scaling);
		func_paper.SetLineColor(kViolet);
		func_paper.Draw("same");
		
		gPad->Modified();
		gPad->Update();
		gPad->WaitPrimitive();
		
		std::cout<<"recording fit results"<<std::endl;
		// for some reason TF1::GetParameter() is returning double values truncated to integer
		// is this something that happens when we fix the sign?
		PushFitAmp(func_sum,"12B");
		PushFitAmp(func_sum,"12N");
		PushFitAmp(func_sum,"16N");
		PushFitAmp(func_sum,"11Be");
		PushFitAmp(func_sum,"9Li");
		PushFitAmp(func_sum,"8He_9C");
		PushFitAmp(func_sum,"8Li_8B");
		PushFitAmp(func_sum,"15C");
		fit_amps["const_4"] = func_sum.GetParameter("const");
		
		// so i guess we'll use the TFitResultPtr as that seems to work...
		for(auto&& anisotope : fit_amps){
			if(anisotope.first.substr(0,5)=="const") continue; // not a real isotope
			int par_number = func_sum.GetParNumber(("amp_"+anisotope.first).c_str());
			double amp = fitresult->Parameter(par_number);
			PushFitAmp(abs(amp),anisotope.first);
			std::cout<<"We had "<<amp<<" "<<anisotope.first<<" decays"<<std::endl;
		}
		
		// add the individual isotopic contributions to reproduce the old plot
		std::vector<TF1> indiv_funcs;
		indiv_funcs.reserve(fit_amps.size());
		for(auto&& theisotope : fit_amps){
			if(theisotope.first.substr(0,5)=="const") continue; // not a real isotope
			std::cout<<"amplitude of isotope "<<theisotope.first<<" is "<<theisotope.second<<std::endl;
			//if(theisotope.second<=0.) continue;
			std::string anisotope = theisotope.first;
			TF1 next_func = BuildFunction({anisotope},0,30);
			PullFitAmp(next_func,anisotope);
			next_func.SetLineColor(colourwheel.GetNextColour());
			
			// this lot might not be necessary ...? SetRangeUser wasn't working, but now it is...?
			double ltime;
			// setrangeuser doesn't seem to be working, so try constraining the x range
			if(anisotope.find("_")==std::string::npos){
				// not a pair
				ltime = lifetimes.at(anisotope);
			} else {
				// pair of isotopes
				// it's not possible to generically solve the necessary equation,
				// but it should be good enough just take the average lifetime
				std::string first_isotope = anisotope.substr(0,anisotope.find_first_of("_"));
				std::string second_isotope = anisotope.substr(anisotope.find_first_of("_")+1,std::string::npos);
				ltime = (lifetimes.at(first_isotope) + lifetimes.at(second_isotope)) / 2.;
			}
			double miny = 1E-2;  // must be less than the plot min or it gets cut off early
			double maxx = -ltime * log(miny*(ltime/theisotope.second));
			if(maxx<0) maxx = 30;
			next_func.SetRange(0.0,maxx);
			
//			std::cout<<anisotope<<" has amplitude "<<theisotope.second<<std::endl;
//			std::cout<<"drawing"<<std::endl;
//			next_func.Draw();
////			next_func.GetYaxis()->SetRangeUser(1E-4,1E3);
//			gPad->Modified();
//			gPad->Update();
//			gPad->WaitPrimitive();
			// this works, but it doesn't add them to the legend!!
			//hist_to_fit.GetListOfFunctions()->Add(&indiv_funcs.back());
			
			indiv_funcs.push_back(next_func);
			indiv_funcs.back().SetName(anisotope.c_str());
			indiv_funcs.back().SetTitle(anisotope.c_str());
		}
		// and last but not least the constant background
		TF1 constfunc("constfunc","[0]",0,30);
		constfunc.SetParameter(0,fit_amps["const_4"]);
		indiv_funcs.push_back(constfunc);
		//hist_to_fit.GetListOfFunctions()->Add(&indiv_funcs.back());
		indiv_funcs.back().SetName("const");
		indiv_funcs.back().SetTitle("const");
		
		// draw it all
		//hist_to_fit.Draw();
		//hist_to_fit.GetYaxis()->SetRangeUser(1E-3,1E5);
		dt_mu_lowe_hist_log.Draw();
		dt_mu_lowe_hist_log.GetYaxis()->SetRangeUser(1.,1E5);
		for(auto&& afunc : indiv_funcs){ std::cout<<"drawing "<<afunc.GetName()<<std::endl; afunc.Draw("same"); }
		gStyle->SetOptStat(0);
		gStyle->SetOptFit(0);
		gPad->GetCanvas()->BuildLegend();
		gPad->Modified();
		gPad->Update();
		//gPad->WaitPrimitive();
		
		// reproduce the paper plot including breakdowns, to check our value extraction
		BuildPaperPlot();
		
		gPad->SetLogx(false);
		gPad->SetLogy(false);
		break;
		}
	default:
		Log(toolName+" FitSubRangeDt invoked with invalid range "+toString(rangenum),v_error,verbosity);
	}
	
	return true;
}

void PurewaterLi9Plots::BuildPaperPlot(){
	// add the individual isotopic contributions to reproduce the old plot
	std::vector<TF1> indiv_funcs;
	indiv_funcs.reserve(papervals.size());
	for(auto&& theisotope : papervals){
		if(theisotope.first.substr(0,5)=="const") continue; // not a real isotope
		if(theisotope.second==0) continue; // do not add to plot isotopes with no abundance
		// we also skip some isotopes since not all are plotted, and our map has some repeats
		if((theisotope.first=="9C") || (theisotope.first=="8He") ||
			//(theisotope.first=="8Li_8B")) continue;                         // if using papervals from plot
			(theisotope.first=="8Li") || (theisotope.first=="8B")) continue;  // if using papervals from table
		
		std::cout<<"amplitude of isotope "<<theisotope.first<<" is "<<theisotope.second<<std::endl;
		std::string anisotope = theisotope.first;
		TF1 next_func = BuildFunction({anisotope},0,30);
		
		PullPaperAmp(next_func,anisotope,false);
		next_func.SetLineColor(colourwheel.GetNextColour());
		
		indiv_funcs.push_back(next_func);
		indiv_funcs.back().SetName(anisotope.c_str());
		indiv_funcs.back().SetTitle(anisotope.c_str());
	}
	// and last but not least the constant background
	TF1 constfunc("constfunc","[0]",0,30);
	constfunc.SetParameter(0,papervals.at("const"));
	indiv_funcs.push_back(constfunc);
	indiv_funcs.back().SetName("const");
	indiv_funcs.back().SetTitle("const");
	
	// now the sum of everything
	TF1 func_paper = BuildFunction({"12B","12N","16N","11Be","9Li","8He_9C","8Li_8B","15C"},0.001,30);
	func_paper.SetName("data");
	func_paper.SetTitle("data");
	std::cout<<"retrieving paper results"<<std::endl;
	PullPaperAmp(func_paper,"12B",false);
	PullPaperAmp(func_paper,"12N",false);
	PullPaperAmp(func_paper,"16N",false);
	PullPaperAmp(func_paper,"11Be",false);
	PullPaperAmp(func_paper,"9Li",false);
	PullPaperAmp(func_paper,"8He_9C",false);
	PullPaperAmp(func_paper,"8Li_8B",false);
	//PullPaperAmp(func_paper,"8Li",false);
	//PullPaperAmp(func_paper,"8B",false);
	PullPaperAmp(func_paper,"15C",false);
	PullPaperAmp(func_paper,"const",false);
	func_paper.SetLineColor(kViolet);
	std::cout<<"Drawing paper total"<<std::endl;
	func_paper.Draw();
	func_paper.GetYaxis()->SetRangeUser(1.,3E7);
//	gPad->Modified();
//	gPad->Update();
//	gPad->WaitPrimitive();
	
	// add all the components
	for(auto&& afunc : indiv_funcs){
		std::cout<<"drawing "<<afunc.GetName()<<std::endl;
		afunc.Draw("same");
//		gPad->Modified();
//		gPad->Update();
//		gPad->WaitPrimitive();
	}
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gPad->GetCanvas()->BuildLegend();
	std::cout<<"Done, drawing paper Fig 3"<<std::endl;
	gPad->Modified();
	gPad->Update();
	gPad->WaitPrimitive();
	
}

// this is the true BuildFunction
TF1 PurewaterLi9Plots::BuildFunction(std::vector<std::string> isotopes, double func_min, double func_max){
	std::string total_func="";
	std::string func_name="";
	std::map<std::string,int> parameter_posns;
	int next_par_index=0;
	for(std::string& anisotope : isotopes){
		func_name += anisotope+"_";
		std::cout<<"adding isotope "<<anisotope<<std::endl;
		// first, this 'isotope' may be a degenerate pair, so try to split it apart
		if(anisotope.find("_")==std::string::npos){
			// not a pair
			std::cout<<"not a pair"<<std::endl;
			int first_index=next_par_index;
			// "([0]/[1])*exp(-x/[1])"
			std::string this_func =  "(abs(["+toString(next_par_index)
									+"])/["+toString(next_par_index+1)
									+"])*exp(-x/["+toString(next_par_index+1)+"])";
			next_par_index +=2;
			// add this function to the total function string
			total_func = (total_func=="") ? this_func : total_func+" + "+this_func;
			// add the parameter names to our map
			parameter_posns.emplace("amp_"+anisotope,first_index);
			parameter_posns.emplace("lifetime_"+anisotope,first_index+1);
		} else {
			// it's a pair. For now, only support two isotopes at a time.
			std::cout<<"a pair, splitting into ";
			std::string first_isotope = anisotope.substr(0,anisotope.find_first_of("_"));
			std::string second_isotope = anisotope.substr(anisotope.find_first_of("_")+1,std::string::npos);
			std::cout<<first_isotope<<" and "<<second_isotope<<std::endl;
			// the fit function isn't just the sum of two single isotope functions
			// as they share an amplitude and constant
			//"[0]*0.5*(exp(-x/[1])/[1]+exp(-x/[2])/[2])
			int first_index=next_par_index;
			std::string this_func =  "abs(["+toString(next_par_index)+"])*0.5*"
									+"(exp(-x/["+toString(next_par_index+1)+"])/" // don't increment index
									+"["+toString(next_par_index+1)+"]+"
									+"exp(-x/["+toString(next_par_index+2)+"])/"  // don't increment index
									+"["+toString(next_par_index+2)+"])";
			next_par_index += 3;
			// add this function to the total function string
			total_func = (total_func=="") ? this_func : total_func+" + "+this_func;
			// add the parameter names to our map
			parameter_posns.emplace("amp_"+anisotope,first_index++);
			parameter_posns.emplace("lifetime_"+first_isotope,first_index++);
			parameter_posns.emplace("lifetime_"+second_isotope,first_index);
		}
	}
	// add the constant term
	total_func += " + abs(["+toString(next_par_index)+"])";
	parameter_posns.emplace("const",next_par_index);
	
	// build the total function from the sum of all strings
	func_name.pop_back(); // remove trailing '_'
	TF1 afunc(func_name.c_str(),total_func.c_str(),func_min,func_max);
	
	// OK, propagate parameter names to the function
	for(auto&& next_par : parameter_posns){
		std::cout<<"settin parname "<<next_par.second<<" to "<<next_par.first<<std::endl;
		afunc.SetParName(next_par.second,next_par.first.c_str());
		if(next_par.first.substr(0,9)=="lifetime_"){
			std::string anisotope = next_par.first.substr(9,std::string::npos);
			std::cout<<"fixing lifetime of "<<anisotope<<std::endl;
			FixLifetime(afunc,anisotope);
		} else {
			// for amplitudes / background constant we'll set the lower limit to be 0.
			// This seems to be necessary for some fits otherwise
			// ROOT gives negative amounts of some isotopes.
			//afunc.SetParLimits(next_par.second, 0, 1E9); // FIXME what to use as upper limit????
			// strictly i think we ought to give initial fit values,
			// but thankfully MINUIT manages without them, but it doesn't like having
			// (default) initial values at the limit, so set some small amount
			//afunc.SetParameter(next_par.second,10);
		}
	}
	
	// return the built function
	return afunc;
}

// this is the hack
TF1 PurewaterLi9Plots::BuildFunction2(std::vector<std::string> isotopes, double func_min, double func_max){
	std::string total_func="";
	std::string func_name="";
	std::map<std::string,int> parameter_posns;
	int next_par_index=0;
	for(std::string& anisotope : isotopes){
		func_name += anisotope+"_";
		std::cout<<"adding isotope "<<anisotope<<std::endl;
		// first, this 'isotope' may be a degenerate pair, so try to split it apart
		if(anisotope.find("_")==std::string::npos){
			// not a pair
			std::cout<<"not a pair"<<std::endl;
			int first_index=next_par_index;
			
			double paperval = papervals.at(anisotope) * (lifetimes.at(anisotope)/0.006);
			// scale to match our max
			paperval *= reco_effs_scaling.at(anisotope); // isotope-dependent scaling from change in E thresh
			paperval *= paper_scaling;                   // a fixed scaling... just in case we want to do this
			std::string paperstring = std::to_string(paperval/2.); // fix half the paper val, fit the rest
			// "((x.xx + abs([0]))/[1])*exp(-x/[1])"
			std::string this_func =  "(("+paperstring+"+abs(["+toString(next_par_index)
									+"]))/["+toString(next_par_index+1)
									+"])*exp(-x/["+toString(next_par_index+1)+"])";
			next_par_index +=2;
			// add this function to the total function string
			total_func = (total_func=="") ? this_func : total_func+" + "+this_func;
			// add the parameter names to our map
			parameter_posns.emplace("amp_"+anisotope,first_index);
			parameter_posns.emplace("lifetime_"+anisotope,first_index+1);
		} else {
			// it's a pair. For now, only support two isotopes at a time.
			std::cout<<"a pair, splitting into ";
			std::string first_isotope = anisotope.substr(0,anisotope.find_first_of("_"));
			std::string second_isotope = anisotope.substr(anisotope.find_first_of("_")+1,std::string::npos);
			std::cout<<first_isotope<<" and "<<second_isotope<<std::endl;
			// the fit function isn't just the sum of two single isotope functions
			// as they share an amplitude and constant
			//"(x.xx + abs([0]))*0.5*(exp(-x/[1])/[1]+exp(-x/[2])/[2])
			int first_index=next_par_index;
			double paperval = papervals.at(first_isotope) * (lifetimes.at(first_isotope)/0.006);
			// scale to match our max
			// isotope-dependent scaling from change in E thresh
			paperval *= 0.5*(reco_effs_scaling.at(first_isotope)+reco_effs_scaling.at(second_isotope));
			// // a fixed scaling... just in case we want to do this
			paperval *= paper_scaling;
			// fix the abundance to at least half the paper val, fit the rest
			std::string paperstring = std::to_string(paperval/2.);
			std::string this_func =  "("+paperstring+"+abs(["+toString(next_par_index)+"]))*0.5*"
									+"(exp(-x/["+toString(next_par_index+1)+"])/" // don't increment index
									+"["+toString(next_par_index+1)+"]+"
									+"exp(-x/["+toString(next_par_index+2)+"])/"  // don't increment index
									+"["+toString(next_par_index+2)+"])";
			next_par_index += 3;
			// add this function to the total function string
			total_func = (total_func=="") ? this_func : total_func+" + "+this_func;
			// add the parameter names to our map
			parameter_posns.emplace("amp_"+anisotope,first_index++);
			parameter_posns.emplace("lifetime_"+first_isotope,first_index++);
			parameter_posns.emplace("lifetime_"+second_isotope,first_index);
		}
	}
	// add the constant term
	total_func += " + abs(["+toString(next_par_index)+"])";
	parameter_posns.emplace("const",next_par_index);
	
	// build the total function from the sum of all strings
	func_name.pop_back(); // remove trailing '_'
	TF1 afunc(func_name.c_str(),total_func.c_str(),func_min,func_max);
	
	// OK, propagate parameter names to the function
	for(auto&& next_par : parameter_posns){
		std::cout<<"settin parname "<<next_par.second<<" to "<<next_par.first<<std::endl;
		afunc.SetParName(next_par.second,next_par.first.c_str());
		if(next_par.first.substr(0,9)=="lifetime_"){
			std::string anisotope = next_par.first.substr(9,std::string::npos);
			std::cout<<"fixing lifetime of "<<anisotope<<std::endl;
			FixLifetime(afunc,anisotope);
		} else {
			// for amplitudes / background constant we'll set the lower limit to be 0.
			// This seems to be necessary for some fits otherwise
			// ROOT gives negative amounts of some isotopes.
			//afunc.SetParLimits(next_par.second, 0, 1E9); // FIXME what to use as upper limit????
			// strictly i think we ought to give initial fit values,
			// but thankfully MINUIT manages without them, but it doesn't like having
			// (default) initial values at the limit, so set some small amount
			//afunc.SetParameter(next_par.second,10);
		}
	}
	
	// return the built function
	return afunc;
}

// fix lifetime of TF1 based on isotope name
void PurewaterLi9Plots::FixLifetime(TF1& func, std::string isotope){
	int par_number = func.GetParNumber(("lifetime_"+isotope).c_str());
	func.FixParameter(par_number,lifetimes.at(isotope));
}

// copy fit result amplitude value back into a component TF1
void PurewaterLi9Plots::PushFitAmp(TF1& func, std::string isotope){
//	std::cout<<"fixing fit result for func "<<func.GetName()<<", isotope "<<isotope<<std::endl;
//	std::cout<<"available pars are: "<<std::endl;
//	for(int i=0; i<func.GetNpar(); ++i){
//		std::cout<<i<<"="<<func.GetParName(i)<<std::endl;
//	}
//	int par_number = func.GetParNumber(("amp_"+isotope).c_str()); // can do straight by name
	double v = static_cast<double>(func.GetParameter(("amp_"+isotope).c_str()));
	std::cout<<"amp of "<<isotope<<" was "<<v<<std::endl;
	fit_amps[isotope] = func.GetParameter(("amp_"+isotope).c_str());
}

// copy fit result amplitude value back into a component TF1
void PurewaterLi9Plots::PushFitAmp(double amp, std::string isotope){
	fit_amps[isotope] = amp;
}

// copy amplitude value from a component TF1 into a sum function for fitting
void PurewaterLi9Plots::PullFitAmp(TF1& func, std::string isotope, bool fix){
	std::string parname = (isotope.substr(0,5)=="const") ? "const" : "amp_"+isotope;
	int par_number = func.GetParNumber(parname.c_str());
	if(fix){
		func.FixParameter(par_number,fit_amps.at(isotope));
	} else {
		func.SetParameter(par_number,fit_amps.at(isotope));
	}
}

void PurewaterLi9Plots::PullPaperAmp(TF1& func, std::string isotope, bool threshold_scaling, double fixed_scaling){
	double paperval = papervals.at(isotope);
	double scaling = 1.;
	if(isotope!="const"){
		if(isotope.find("_")==std::string::npos){
			paperval *= (lifetimes.at(isotope) ); ///0.006);
			// not a pairing - can look up efficiency directly
			if(threshold_scaling) scaling = reco_effs_scaling.at(isotope);
			scaling *= fixed_scaling;
		} else {
			std::string first_isotope = isotope.substr(0,isotope.find_first_of("_"));
			std::string second_isotope = isotope.substr(isotope.find_first_of("_")+1,std::string::npos);
			paperval *= (((lifetimes.at(first_isotope)+lifetimes.at(second_isotope))*0.5)); ///0.006);
			// how do we handle efficiency for pairs?
			// best we can do is assume half for each and average the efficiency i think....
			if(threshold_scaling){
				scaling = 0.5*(reco_effs_scaling.at(first_isotope)
										 + reco_effs_scaling.at(second_isotope));
			}
			scaling *= fixed_scaling;
		}
	}
	// scale to match our max
	//std::cout<<"fixed scaling is "<<fixed_scaling<<std::endl;
	//std::cout<<"efficiency for "<<isotope<<" is "<<(scaling/fixed_scaling)<<std::endl;
	paperval *= scaling;
	std::string parname = (isotope=="const") ? "const" : "amp_"+isotope;
	func.SetParameter(parname.c_str(), paperval);
}

bool PurewaterLi9Plots::PlotLi9BetaEnergy(){
	// events passing dlt cut, li9 energy and lifetime cuts, and with a tagged neutron passing BDT cut
	TH1F li9_e_hist("li9_e_hist","Li9 Candidate Beta Energy",7,6,li9_endpoint);
	// as we do this we need to add rest mass as reconstructed energy is only kinetic...
	double e_rest_mass = 0.511; // [MeV] TODO replace this with TParticleDatabase lookup
	for(auto&& aval : li9_e_vals) li9_e_hist.Fill(aval+e_rest_mass); // FIXME weight by # post mus & num neutrons
	li9_e_hist.Write();
	
	// again to overlay the expected li9 and background plots, we need to take the beta spectra
	// of each isotope and propagate it through the efficiency chain TODO
	
	return true;
}

// =========================================================================
// Li9 ncapture fits
// =========================================================================

bool PurewaterLi9Plots::PlotNcaptureDt(){
	// events passing dlt cut, li9 energy and lifetime cuts, and with a tagged neutron passing BDT cut
	
	// make a histogram to bin the data
	std::cout<<"making li9 lowe-ncapture dt histogram"<<std::endl;
	TH1F li9_ncap_dt_hist("li9_ncap_dt_hist","Beta to ncapture dt for Li9 triplets",
	                       21,0,500E-6);
	
	// XXX debug investigation - no sign of expl decay of ncapture times... time span is correct,
	// there's no data beyond 500us. What's going on?
	std::cout<<"first 100 ncapture times were: {";
	int ncpi=0;
	for(auto&& aval : li9_ntag_dt_vals){
		// XXX FIXME REMOVE AFTER REPROCESSING IN ANALYSE XXX XXX XXX XXX XXX XXX 
		double ncap_time_adjusted = aval < 50000 ? aval : aval - 65000;
		ncap_time_adjusted /= 1E9;
		if(ncpi<100) std::cout<<ncap_time_adjusted<<", "; ++ncpi;
		li9_ncap_dt_hist.Fill(ncap_time_adjusted);  // FIXME weight by num_post_muons and num neutrons
	}
	std::cout<<"}"<<std::endl;
	std::cout<<"saving to file"<<std::endl;
	// for comparison to the paper, scale the x axis up to us
	li9_ncap_dt_hist.GetXaxis()->SetLimits(0,500);  // changes axis labels but doesn't affect binning
	li9_ncap_dt_hist.Write();
	// set it back for analysis
	li9_ncap_dt_hist.GetXaxis()->SetLimits(0,500E-6);  // XXX this is bad practice
	
	// TODO get expected number of background events
	// this is also needed for the expected background E spectrum, where we also need
	// the energies of the background, so calculate that first and we can just take the number
	
	// fit the ncapture lifetime, Fig 5 from the paper.
	// do we do chi2 fit or binned likelihood fit? why? FIXME how do we do binned likelihood?
	
	// perform an unbinned chi2 fit
	std::cout<<"doing binned chi2 fit"<<std::endl;
	double binned_estimate = BinnedNcapDtChi2Fit(&li9_ncap_dt_hist);
	
	// perform an unbinned likelihood fit
	std::cout<<"doing unbinned likelihood fit"<<std::endl;
	UnbinnedNcapDtLogLikeFit(&li9_ncap_dt_hist, binned_estimate);
	
	return true;
}

double PurewaterLi9Plots::BinnedNcapDtChi2Fit(TH1F* li9_ncap_dt_hist){
	// this fits the lifetime of ncapture to extract the amount of exponential and constant
	std::cout<<"making TF1 for binned chi2 fit with "<<li9_ncap_dt_hist->GetEntries()<<" values"<<std::endl;
	
	// number of neutrons left after time t = N = N0*exp(-t/τ)
	// rate of change of num neutrons, i.e. rate of observed ncapture events = dN/dt = -N0*τ*exp(-t/τ)
	// we can ignore the sign which just says N is decreasing.
	// simple chi2 fit of 'y = C + A*τ*exp(-dt/τ)' for all values of C, A, τ
	TF1 ncap_dt_func("ncap_dt_func","[0]+[1]*[2]*exp(-x/[2])",ncap_dtmin,ncap_dtmax);
	
	// set parameter names
	ncap_dt_func.SetParName(0,"rate of bg events");
	ncap_dt_func.SetParName(1,"num of Li9+n events");
	ncap_dt_func.SetParName(2,"ncapture lifetime");
	
	// set starting values - we can fix the ncapture lifetime
	//ncap_dt_func.SetParameters(allparams.data());  // pass an array or set individually
	ncap_dt_func.SetParameter(0,0);  // TODO estimate from accidental rate of ntag from MC and num events?
	ncap_dt_func.SetParameter(1,li9_ncap_dt_hist->GetBinContent(1));
	ncap_dt_func.FixParameter(2,lifetimes.at("ncapture"));
	
	// set num TF1 points for drawing
	ncap_dt_func.SetNpx(1000);
	
	// DO THE FIT
	std::cout<<"invoking li9 ncapture binned chi2 fit"<<std::endl;
	TFitResultPtr fitresult = li9_ncap_dt_hist->Fit("ncap_dt_func","MRS");
	// options here are  M: better fit, R: use range, S: get resultsptr
	
	// print result
	std::cout<<"ncapture dt binned chi2 fit parameters = {";
	for(int i=0; i<(ncap_dt_func.GetNpar()); i++){
		std::cout<<ncap_dt_func.GetParameter(i);
		(i<(ncap_dt_func.GetNpar()-1)) ? std::cout<<", " : std::cout<<"};\n";
	}
	//float fitchi2 = fitresult->Chi2();                 // same as below, which
	float fitchi2 = ncap_dt_func.GetChisquare();         // doesn't need fitresultptr
	Log(toolName+" li9 lowe->ncap dt fit chi2 was "+toString(fitchi2),v_message,verbosity);
	
	// draw result
	li9_ncap_dt_hist->Draw();
	gPad->WaitPrimitive();
	gPad->Clear();
//	//fit->Draw("lsame");   // not necessary, fit is added to histogram's list of functions and drawn automatically
	li9_ncap_dt_hist->GetListOfFunctions()->Clear();
	
	std::cout<<"binned ncapture chi2 fit done"<<std::endl;
	return ncap_dt_func.GetParameter(1);
}

bool PurewaterLi9Plots::UnbinnedNcapDtLogLikeFit(TH1F* li9_ncap_dt_hist, double num_li9_events){
	
	std::cout<<"doing unbinned likelihood fit with "<<li9_ntag_dt_vals.size()<<" values"<<std::endl;
	
	// TF1 of likelihood distribution
	// ROOT provides some basic pdfs: https://root.cern.ch/doc/v610/group__PdfFunc.html
	// but we need our likelihood to be normalized and I don't see how to do it suitably with these
//	TF1 ncapture_dt_unbinned_like("ncapture_dt_unbinned_like",
//	TString::Format("[0]*ROOT::Math::uniform_pdf (x,%.2f,%.2f,0) "
//	                " + [1]*ROOT::Math::exponential_pdf(x,[1],0)",xmin,xmax),xmin,xmax);
	
	// so let's make our own
	int fit_n_pars = 1;
	int func_n_dims = 1;
	std::cout<<"making TF1"<<std::endl;
	TF1 ncapture_dt_unbinned_like("ncapture_dt_unbinned_like",this,&PurewaterLi9Plots::ncap_lifetime_loglike,
	        ncap_dtmin, ncap_dtmax, fit_n_pars, "PurewaterLi9Plots", "ncap_lifetime_loglike");
	
	// set parameter names
	ncapture_dt_unbinned_like.SetParName(0,"fraction of background events");
	ncapture_dt_unbinned_like.SetParLimits(0,0,1.);
	
	//ncapture_dt_unbinned_like.SetParameters(allparams.data());  // pass an array or set individually
	// set starting values to fit val from binned fit
	std::cout<<"setting starting value of background fraction to "
			 <<(1.-(num_li9_events/li9_ntag_dt_vals.size()))<<std::endl;
	ncapture_dt_unbinned_like.SetParameter(0,(1.-(num_li9_events/li9_ntag_dt_vals.size())));
	
	int nbins=10;
	for(int i=nbins; i>0; --i){
		std::vector<double> pars{ncapture_dt_unbinned_like.GetParameter(0)};
		double xval = ncap_dtmin+double(i)*((ncap_dtmax-ncap_dtmin)/double(nbins-1));
		std::cout<<"ncapture unbinned function eval at "<<xval<<" = "
				 <<ncapture_dt_unbinned_like.Eval(xval)<<std::endl
				 <<", base function = "<<ncap_lifetime_loglike(&xval,pars.data())
				 <<std::endl;
	}
	
	// draw for comparison. Since our likelihood function is normalised we need to normalise the data too
	TH1F* li9_ncap_dt_hist_normalised = (TH1F*)li9_ncap_dt_hist->Clone();
	li9_ncap_dt_hist_normalised->Scale(1./li9_ncap_dt_hist->Integral());
	li9_ncap_dt_hist_normalised->Draw();
	ncapture_dt_unbinned_like.Draw("same");
	gPad->WaitPrimitive();
	return true;
	
	// convert the TF1 for use in unbinned fit
	std::cout<<"making TF1 wrapper"<<std::endl;
	ROOT::Math::WrappedMultiTF1 ncapture_dt_unbinned_func(ncapture_dt_unbinned_like,
	                                                      ncapture_dt_unbinned_like.GetNdim());
	
	// build a fitter
	// 'false' says let the fitter calculate derivatives of the function itself (recommended)
	std::cout<<"making fitter"<<std::endl;
	ROOT::Fit::Fitter fitter;
	fitter.SetFunction(ncapture_dt_unbinned_func);  // no bool in ROOT 5 - what's the default?
	// important note: the fit function should be normalized, but simply calculating and dividing by
	// the integral before returning may not be enough as this introduces a correlation
	// between the fit parameters. Apparently the best way is to fit N-1 parameters, then derive
	// the value of the remaining one from the others and a fixed normalization condition
	// https://root-forum.cern.ch/t/root-fit-unbindata-example-is-not-working/28462/3
	// in cases where it is sufficient to simply divide by the integral (which??) 
	// a wrapper is given in https://root-forum.cern.ch/t/unbinned-log-likelihood-fits/18981/2
	
	// can we do it just by fixing the parameter?
	// https://web.pa.msu.edu/people/brock/file_sharing/ATLAS/root/math/mathcore/test/fit/testRooFit.cxx
	//fitter.Config().ParSettings(0).Fix();
	
	// configure the fitter
	// there is a fitter.Config().ParSettings(#) for each parameter #
	// these are available only after calling ROOT::Fit::Fitter::SetFunction
	// the configurable options are:
	// the initial values of the parameters;   (can also be set by input model func)
	// the parameter step sizes;
	// eventual parameter bounds;
	// the minimizer library and the particular algorithm to use;
	// different minimization options (print level, tolerance, max iterations, etc.)
	// the type of parameter errors to compute (parabolic error, Minos errors, re-normalized errors using fitted chi2 values)
	
	// set initial values. Think it also inherits initial values from TF1 on construction?
	//double initialParams[] = {2,1};
	//fitter.Config().SetParamsSettings(2,initialParams);
	
	// set the minimizer algorithm to use
	// https://root.cern.ch/root/htmldoc/guides/users-guide/FittingHistograms.html
	// https://root.cern.ch/root/htmldoc/ROOT__Fit__FitConfig.html
	fitter.Config().SetMinimizer("Minuit2","Migrad");
//	fitter.Config().SetUpdateAfterFit();  // not in root 5
	
	// minimizer options include:
	// Minimizer type (Minuit, Fumili, GSLMultiMin...)
	// Minimizer algorithm (Migrad, Simplex, Scan...)
	// Strategy - Minuit default 1 is to only compute full Hessian matrix after minimization(?)
	// Print level (verbosity 0...)
	// Tolerance to control iterations
	// Max function calls
	// Max iterations (not used by Minuit)
	// Parameter errors value (??) - default 1 for chi2, 0.5 for log-likelihood
	// Precision in evaluation of minimization. Default double. Should we use float?
	
	// MinimizerOptions::SetMaxFunctionCalls(int )
	// MinimizerOptions::SetTolerance(double )
	// MinimizerOptions::SetPrecision(double )
	//ROOT::Math::MinimizerOptions opt = fitter.Config().MinimizerOptions();
	//opt.SetMaxFunctionCalls(5000);  // XXX like this?
	
//	ROOT::Fit::DataOptions opt;
//	ROOT::Fit::DataRange range(0,5);
//	ROOT::Fit::UnBinData ncap_dt_data(opt, range, x.size());
	
	// convert the data into a suitable 'UnBinData' object
	std::cout<<"making unbinned dataset"<<std::endl;
	ROOT::Fit::UnBinData ncap_dt_data(li9_ntag_dt_vals.size());
	for(auto aval : li9_ntag_dt_vals){   // note: use auto NOT auto&&
		ncap_dt_data.Add(aval);   // can we introduce weights here? FIXME
	}
	
	// DO THE FIT
	std::cout<<"doing the fit"<<std::endl;
	fitter.LikelihoodFit(ncap_dt_data);
	
	// print the results
	ROOT::Fit::FitResult r=fitter.Result();
	r.Print(std::cout);
	
	std::cout<<"unbinned likelihood fit done"<<std::endl;
	
	// TODO err, extract and return parameters for drawing
	li9_ncap_dt_hist_normalised->Draw();
	ncapture_dt_unbinned_like.Draw("same");
	gPad->WaitPrimitive();
	gPad->Clear();
	return true;
}

double PurewaterLi9Plots::ncap_lifetime_loglike(double* x, double* par){
	// shouldn't be trying to fit times outside our range
	if(((*x)<ncap_dtmin) || ((*x)>ncap_dtmax)) return 1e10;
	// shouldn't be trying to fit parameter values outside the valid range
	if((par[0]<0.) || (par[0]>1.)) return 1e10;
	
	// likelihood of an event occurring at a given time, using the specified background fraction
	// being that this is a probability the return value should be normalised 0-1
	// since the integral is fixed, all we can vary here is is the fraction
	// of events that are signal vs background
	
	// rate = C + A*exp(-dt/τ)
	// τ = neutron capture lifetime
	// C = rate of accidental backgrounds
	// A = rate of signal events
	
	// we need to make a PDF, so need to normalise. Calculating the integral we get:
	// integral = C*Δt + A*τ*[exp(-tmax/τ) - exp(-tmin/τ)]
	// the normalized function is obtained by scaling by this integral, turning C->c, A->a
	// a = A / integral, c = C / integral, giving:
	// c*Δt + a*τ*[exp(-tmax/τ) - exp(-tmin/τ)] = 1
	// we can relate a and c:
	// a = (c*Δt - 1 ) / τ*[exp(-tmax/τ) - exp(-tmin/τ)]
	// so now c*Δt is the fraction of background events
	// for a to remain poaitive we should restrict c to the range 0 to (1/Δt)
	// if we define par[0] = c' = c*Δt, then we can fit it directly.
	
	// Now, we have been given a value par[0], so
	// calculate a based on our value of c':
	double timespan = ncap_dtmax-ncap_dtmin;
	double ncap_lifetime = lifetimes.at("9Li");
	double a =  (par[0] - 1.) / (ncap_lifetime*(exp(-ncap_dtmax/ncap_lifetime) - exp(-ncap_dtmin/ncap_lifetime)));
	
	// calculate likelihood of an event at time (*x) given a background fraction of c'
	// i.e. evaluate the normalised function, c + a*exp(-dt/τ), at this time
	double thelikelihood = (par[0]/timespan) + a*exp(-(*x)/ncap_lifetime);
	
	// return log-likelihood
	return thelikelihood;
	//return log(thelikelihood);    // doesn't this mess with our normalization?
}


// =========================================================================
// Li9 lifetime fits
// =========================================================================

bool PurewaterLi9Plots::PlotLi9LifetimeDt(){
	
	// make a histogram to bin the data
	std::cout<<"making li9 mu-lowe dt histogram"<<std::endl;
	TH1F li9_muon_dt_hist("li9_muon_dt_hist","Muon to Low-E dt for Li9 triplets",
	                       15,li9_lifetime_dtmin,li9_lifetime_dtmax);
	for(auto&& aval : li9_muon_dt_vals){
		li9_muon_dt_hist.Fill(aval);
	}
	std::cout<<"saving to file"<<std::endl;
	li9_muon_dt_hist.Write();
	
	/*
	// TODO calculate the efficiency of background selection and scale down the number
	// of spallation events to get the expected number of accidental triplets
	// efficiency is based on background acceptance of BDT and of preceding cuts i guess
	// this is subtracted off to estimate the number of li9 events
	
	// get ntag signal and background efficiency
	// ------------------------------------------
	// read the BDT ROC and populate the splines
	fill_roc_ntag();
	// evaluate the splines at the working point used to retrieve the efficiencies
	double ntag_sig_eff = sig->Eval(ncut_li9);
	double ntag_bg_eff = bg->Eval(ncut_li9);
	
	// total background selection efficiency is product of all cuts specific to Li9
	// * dt in li9 lifetime range 0.05-0.5
	// * lowe energy in li9 beta energy range
	// * no other muons within 1ms
	// * ntag bdt fom > 0.995
	// how do we determine the fraction of non-li9 events that pass these? MC? paper doesn't say...
	double li9_lifetime_bg_eff = 1;
	double li9_energy_bg_eff = 1;
	double no_mu_1ms_bg_eff = 1;
	double li9_bg_eff = li9_lifetime_bg_eff * li9_energy_bg_eff * ntag_bg_eff * no_mu_1ms_bg_eff;
	// we multiply this by the number of spallation events to get the total bg contamination
	double num_spall_events = dt_mu_lowe_vals.size();
	double total_li9_bg = num_spall_events * li9_bg_eff;
	// and then we get the number of li9 events as the remainder
	// fig 4 plots the spectrum of li9 energy, and of accidentals - presumably from MC of all bgs?
	// and the expected li9 spectrum - presumably from some suitable reference
	*/
	
	// fit the Li9 lifetime, Fig 6 from the paper. NOT neutron capture times (plot 5!)
	// this fits a combination of 8 exponentials (including li9) with fixed abundances,
	// and uses it to extract the lifetime of li9. The paper only does a binned fit.
	// do we do chi2 fit or binned likelihood fit? why? FIXME how do we do binned likelihood?
	
	// perform a binned chi2 fit
	std::cout<<"doing Li9 lifetime binned chi2 fit"<<std::endl;
	double binned_estimate = BinnedLi9DtChi2Fit(&li9_muon_dt_hist);
	
	// TODO (un)binned extended likelihood fit to also extract number of Li9 events?
	
	return true;
}

double PurewaterLi9Plots::BinnedLi9DtChi2Fit(TH1F* li9_muon_dt_hist){
	// this fits the lifetime of li9 as a cross-check, by fitting 8 exponentials
	// based on backgrounds from other isotopes. The amount of other isotopes is based
	// on the previous global fit to muon-lowe dt, and the fraction of those that pass
	// all cuts, from first red, third, lt, li9 dt, li9 energy and ntag
	// To do this we need the efficiency of all of those cuts for all isotopes!
	// finally we fit (unbinned? binned?) the distribution of mu-lowe times leaving only li9 lifetime floating
	
	// currently the above paragraph does not describe what this code does FIXME
	std::cout<<"making TF1 for binned chi2 fit with "<<li9_muon_dt_hist->GetEntries()<<" values"<<std::endl;
	
	// number of li9 left after time t = N = N0*exp(-t/τ)
	// rate of change of li9, i.e. rate of observed li9 events = dN/dt = -N0*τ*exp(-t/τ)
	// we can ignore the sign which just says N is decreasing.
	// simple chi2 fit of 'y = C + A*τ*exp(-dt/τ)' for all values of C, A, τ
	TF1 li9_muon_dt_func("li9_muon_dt_func","[0]+[1]*[2]*exp(-x/[2])",
	                      li9_lifetime_dtmin,li9_lifetime_dtmax);
	
	// set parameter names
	li9_muon_dt_func.SetParName(0,"num bg events");
	li9_muon_dt_func.SetParName(1,"num Li9 events");
	li9_muon_dt_func.SetParName(2,"Li9 lifetime");
	
	// set starting values - we can fix the lifetime
	//li9_muon_dt_func.SetParameters(allparams.data());  // pass an array or set individually
	li9_muon_dt_func.SetParameter(0,0);  // TODO estimate from accidental rate of ntag from MC and num events?
	li9_muon_dt_func.SetParameter(1,li9_muon_dt_hist->GetBinContent(1));
	li9_muon_dt_func.FixParameter(2,lifetimes.at("9Li"));
	
	// set parameter limits. For minuit at least this is strongly discouraged
	// Maybe ok for chi2 fits but not likelihood fits?
	//li9_muon_dt_func.SetParLimits(4, 0., 200.);
	
	// set num TF1 points for drawing
	li9_muon_dt_func.SetNpx(1000);
	
	// DO THE FIT
	std::cout<<"invoking li9 lifetime chi2 fit"<<std::endl;
	TFitResultPtr fitresult = li9_muon_dt_hist->Fit("li9_muon_dt_func","MRS");
	// options here are  M: better fit, R: use range, S: get resultsptr
	
	// print result
	std::cout<<"li9 lifetime binned chi2 fit parameters = {";
	for(int i=0; i<(li9_muon_dt_func.GetNpar()); i++){
		std::cout<<li9_muon_dt_func.GetParameter(i);
		(i<(li9_muon_dt_func.GetNpar()-1)) ? std::cout<<", " : std::cout<<"};\n";
	}
	//float fitchi2 = fitresult->Chi2();                               // same as below, which
	float fitchi2 = li9_muon_dt_func.GetChisquare();         // doesn't need fitresultptr
	Log(toolName+" li9 mu->lowe dt fit chi2 was "+toString(fitchi2),v_message,verbosity);
	
	// draw result
	/*
	li9_muon_dt_hist->Draw();
	gPad->WaitPrimitive();
	gPad->Clear();
//	//fit->Draw("lsame");   // not necessary, fit is added to histogram's list of functions and drawn automatically
	*/
	li9_muon_dt_hist->GetListOfFunctions()->Clear();
	
	std::cout<<"li9 lifetime binned chi2 fit done"<<std::endl;
	return li9_muon_dt_func.GetParameter(1);
}

// =========================================================================
// =========================================================================

bool PurewaterLi9Plots::GetEnergyCutEfficiencies(){
	// 2015 paper had a low threshold of 6MeV, newer data has (as of now) 8MeV
	// so to compare yeilds we need to scale the paper values by the fraction
	// of decays that would be above the respective thresholds
	
	std::map<std::string,int> true_events_below_8MeV;
	std::map<std::string,int> true_events_above_8MeV;
	std::map<std::string,int> true_events_below_6MeV;
	std::map<std::string,int> true_events_above_6MeV;
	
	// same but using energy from bonsai
	std::map<std::string,int> reco_events_below_8MeV;
	std::map<std::string,int> reco_events_above_8MeV;
	std::map<std::string,int> reco_events_below_6MeV;
	std::map<std::string,int> reco_events_above_6MeV;
	
	BoostStore efficiencyStore(true,constants::BOOST_STORE_BINARY_FORMAT);
	efficiencyStore.Initialise(efficienciesFile);
	efficiencyStore.Get("true_events_below_8MeV",true_events_below_8MeV);
	efficiencyStore.Get("true_events_above_8MeV",true_events_above_8MeV);
	efficiencyStore.Get("true_events_below_6MeV",true_events_below_6MeV);
	efficiencyStore.Get("true_events_above_6MeV",true_events_above_6MeV);
	
	efficiencyStore.Get("reco_events_below_8MeV",reco_events_below_8MeV);
	efficiencyStore.Get("reco_events_above_8MeV",reco_events_above_8MeV);
	efficiencyStore.Get("reco_events_below_6MeV",reco_events_below_6MeV);
	efficiencyStore.Get("reco_events_above_6MeV",reco_events_above_6MeV);
	
	for(auto&& anisotope : true_events_below_6MeV){
		std::string isotope = anisotope.first;
		double true_below6 = true_events_below_6MeV.at(isotope);
		double true_above6 = true_events_above_6MeV.at(isotope);
		double true_eff6 = true_above6/(true_below6+true_above6);
		double true_below8 = true_events_below_8MeV.at(isotope);
		double true_above8 = true_events_above_8MeV.at(isotope);
		double true_eff8 = true_above8/(true_below8+true_above8);
		std::cout<<"TRUE ENERGY:"<<std::endl;
		std::cout<<"Isotope "<<anisotope.first<<" had "<<true_below6<<" events below 6MeV vs "
				 <<true_above6<<" above 6 MeV corresponding to an efficiency of "<<(true_eff6*100.)<<"%"
				 <<std::endl<<"compared to "<<true_below8<<" events below 8MeV and "<<true_above8
				 <<" events above 8 MeV corresponding to an efficiency of "<<(true_eff8*100.)<<"%"<<std::endl;
		
		true_effs_6mev.emplace(isotope,true_eff6);
		true_effs_8mev.emplace(isotope,true_eff8);
		true_effs_scaling.emplace(isotope,(true_eff8/true_eff6));
		
		// same with reconstructed energies - the spectra are VERY distorted...!?
		double reco_below6 = reco_events_below_6MeV.at(isotope);
		double reco_above6 = reco_events_above_6MeV.at(isotope);
		double reco_eff6 = reco_above6/(reco_below6+reco_above6);
		double reco_below8 = reco_events_below_8MeV.at(isotope);
		double reco_above8 = reco_events_above_8MeV.at(isotope);
		double reco_eff8 = reco_above8/(reco_below8+reco_above8);
		std::cout<<"RECONSTRUCTED ENERGY:"<<std::endl;
		std::cout<<"Isotope "<<anisotope.first<<" had "<<reco_below6<<" events below 6MeV vs "
				 <<reco_above6<<" above 6 MeV corresponding to an efficiency of "<<(reco_eff6*100.)<<"%"
				 <<std::endl<<"compared to "<<reco_below8<<" events below 8MeV and "<<reco_above8
				 <<" events above 8 MeV corresponding to an efficiency of "<<(reco_eff8*100.)<<"%"<<std::endl;
		
		reco_effs_6mev.emplace(isotope,reco_eff6);
		reco_effs_8mev.emplace(isotope,reco_eff8);
		reco_effs_scaling.emplace(isotope,(reco_eff8/reco_eff6));
	}
	return true;
}
















bool PurewaterLi9Plots::PlotPaperDt(THStack& ourplots){
	// get the first of our plots to use as a base to define binning
	TH1F* basehist = (TH1F*)ourplots.GetHists()->At(0);
	
	TH1F* h_paper_dt_mu_nonfitted = (TH1F*)basehist->Clone("paper_dt_nonfitted");
	h_paper_dt_mu_nonfitted->Reset();
	// digitized from paper. I didn't do all series' as most overlapped.
	std::vector<double> paper_dt_mu_nonfitted
		{ 0.5458, 0.3111, 0.1014, 0.0403, 0.0014, 0.0806, 0.0694, 0.0111, 0.0014, 0.0056};
	for(int i=0; i<paper_dt_mu_nonfitted.size(); ++i){
		h_paper_dt_mu_nonfitted->SetBinContent(i+1, paper_dt_mu_nonfitted.at(i));
	}
	h_paper_dt_mu_nonfitted->SetLineColor(kBlue);
	h_paper_dt_mu_nonfitted->SetLineStyle(2);
	
	TH1F* h_paper_dt_mu_single = (TH1F*)basehist->Clone("paper_dt_single");
	h_paper_dt_mu_single->Reset();
	std::vector<double> paper_dt_mu_single
		{ 0.5806, 0.2153, 0.0903, 0.0403, 0.0236, 0.0153, 0.0111, 0.0083, 0.0097, 0.0069 };
	for(int i=0; i<paper_dt_mu_single.size(); ++i){
		h_paper_dt_mu_single->SetBinContent(i+1, paper_dt_mu_single.at(i));
	}
	h_paper_dt_mu_single->SetLineColor(kBlack);
	h_paper_dt_mu_single->SetLineStyle(2);
	
	/*
	// put them into a stack
	THStack paper_dts;
	paper_dts.Add(h_paper_dt_mu_nonfitted);
	paper_dts.Add(h_paper_dt_mu_single);
	// write them out to file
	paper_dts.Write("spall_dts_paper");
	// cleanup
	paper_dts.GetHists()->Delete();
	*/
	
	// add to our stack for easy combined plotting
	ourplots.Add(h_paper_dt_mu_nonfitted);
	ourplots.Add(h_paper_dt_mu_single);
	
	return true;
}

bool PurewaterLi9Plots::PlotPaperDlt(THStack& ourplots){
	
	// get first one as a base for defining axes
	TH1F* basehist = (TH1F*)ourplots.GetHists()->At(0);
	
	// digitized from paper
	std::vector<double> paper_dl_multiple
		{ 0.1414, 0.2303, 0.1980, 0.1465, 0.1121, 0.0798, 0.0596, 0.0404 };
	TH1F* h_paper_dl_multiple = (TH1F*)basehist->Clone("paper_dl_multiple");
	h_paper_dl_multiple->Reset();
	for(int i=0; i<paper_dl_multiple.size(); ++i){
		h_paper_dl_multiple->SetBinContent(i+1,paper_dl_multiple.at(i));
	}
	h_paper_dl_multiple->SetLineColor(kRed);
	h_paper_dl_multiple->SetLineStyle(2);
	
	std::vector<double> paper_dl_single
		{ 0.3636, 0.3707, 0.1485, 0.0636, 0.0283, 0.0141, 0.0121, 0.0071 };
	TH1F* h_paper_dl_single = (TH1F*)basehist->Clone("paper_dl_single");
	h_paper_dl_single->Reset();
	for(int i=0; i<paper_dl_single.size(); ++i){
		h_paper_dl_single->SetBinContent(i+1,paper_dl_single.at(i));
	}
	h_paper_dl_single->SetLineColor(kBlack);
	h_paper_dl_single->SetLineStyle(2);
	
	std::vector<double> paper_dl_stopping
		{ 0.4161, 0.3242, 0.0455, 0.0414, 0.0707, 0.0434, 0.0263, 0.0303 };
	TH1F* h_paper_dl_stopping = (TH1F*)basehist->Clone("paper_dl_stopping");
	h_paper_dl_stopping->Reset();
	for(int i=0; i<paper_dl_stopping.size(); ++i){
		h_paper_dl_stopping->SetBinContent(i+1,paper_dl_stopping.at(i));
	}
	h_paper_dl_stopping->SetLineColor(kMagenta);
	h_paper_dl_stopping->SetLineStyle(2);
	
	/*
	// put them into a stack
	THStack paper_dls;
	paper_dls.Add(h_paper_dl_multiple);
	paper_dls.Add(h_paper_dl_single);
	paper_dls.Add(h_paper_dl_stopping);
	// write them out to file
	paper_dls.Write("spall_dls_paper");
	// cleanup
	paper_dls.GetHists()->Delete();
	*/
	
	// add them to our stack for combined plotting
	ourplots.Add(h_paper_dl_multiple);
	ourplots.Add(h_paper_dl_single);
	ourplots.Add(h_paper_dl_stopping);
	
	return true;
}

void PurewaterLi9Plots::fill_roc_ntag(){
	ifstream fin(bdt_outfile);
	const int splmax = 5000;
	double cut, sg, bkg, sgsys, bkgsys, dum1, dum2;
	std::vector<double> v_cut, v_sg, v_bkg, v_sgsys, v_bkgsys;
	int k = 0;
	while(fin){
		std::vector<double> res;
		fin >> cut >> sg >> bkg >> sgsys >> bkgsys >> dum1 >> dum2;
		if (cut > cut_lowth) break;
		v_cut.push_back(cut);
		v_sg.push_back(sg);
		v_bkg.push_back(bkg);
		v_sgsys.push_back(sgsys);
		v_bkgsys.push_back(bkgsys);
	}
	sig = new TSpline3("signal", v_cut.data(), v_sg.data(),v_cut.size(), 0, 1, 0);
	bg = new TSpline3("background", v_cut.data(), v_bkg.data(), v_cut.size(), 0, 1, 0);
	sigsys = new TSpline3("signal systematic", v_cut.data(), v_sgsys.data(), v_cut.size(), 0, 0, 0);
	bgsys = new TSpline3("background systematic", v_cut.data(), v_bkgsys.data(), v_cut.size(), 0, 0, 0);
}

