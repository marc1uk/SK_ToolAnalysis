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
	catch(std::exception& e){
		Log(toolName+" encountered error "+e.what()+" during Analyse() at "+__FILE__+"::"
			+std::to_string(__LINE__),v_error,verbosity);
		// __FILE__ gives the current file name
		// __func__ gives current function name
		// __LINE__ gives current line
		// unfortunately these do not refer to the throw site, but rather the site at which they're used
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
	for(size_t mu_i : pre_muboy_first_muons){
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
	for(size_t mu_i : post_muboy_first_muons){
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
	
	// pass out maps for Pre- and Post-Muon Dt / Dlt distributions
	// PlotMuonDtDlt tool will create plots to compare to paper Fig 2
	m_data->vars.Set("dlt_vals_pre",&dlt_vals_pre);
	m_data->CStore.Set("dt_vals_pre",&dt_vals_pre,false);
	m_data->CStore.Set("dlt_vals_post",&dlt_vals_post,false);
	m_data->CStore.Set("dt_vals_post",&dt_vals_post,false);
	
	// pass out mu-lowe time diffs for li9 triplet candidates
	// FitLi9Lifetime tool will produce plot to compare to paper Fig 6
	m_data->CStore.Set("li9_muon_dt_vals",&li9_muon_dt_vals,false);
	
	// pass out the mu-lowe time distribution for spallation
	// FitSpallationDt tool will create plots to compare to paper Fig 3
	m_data->CStore.Set("dt_mu_lowe_vals",&dt_mu_lowe_vals,false);
	m_data->CStore.Set("livetime",livetime);
	
	// pass out li9->ntag dt values
	// FitPurewaterLi9NcaptureDt tool will create plots to compare to paper Fig 5
	m_data->CStore.Set("li9_ntag_dt_vals",&li9_ntag_dt_vals,false);
	
	// TODO move this to another tool.
	// XXX tool will create plots to compare to paper Fig 4
	Log(toolName+" making plots of Li9 candidate beta spectrum",v_debug,verbosity);
	PlotLi9BetaEnergy();
	
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

void PurewaterLi9Plots::fill_roc_ntag(){
	ifstream fin(bdt_outfile);
	double cut_lowth = 0.999; // i dunno, comes from cut_third.cpp in li9 stuff FIXME make a config variable...
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

