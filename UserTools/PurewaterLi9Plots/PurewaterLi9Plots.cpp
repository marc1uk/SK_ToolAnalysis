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
#include "TMath.h"
#include "TString.h"
#include "TSpline.h"
#include "TEntryList.h"

// For Fitting
#include "Fit/Fitter.h"
//#include "Fit/BinData.h"
#include "Fit/UnBinData.h"
//#include "Fit/Chi2FCN.h"
#include "Fit/FitResult.h"
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
	//* unused when an upstream tool is driving the toolchain
	//* when using values from a BoostStore the ROOT files do not need to be read
	
	// cut thresholds
	m_variables.Get("li9_lifetime_dtmin",li9_lifetime_dtmin);
	m_variables.Get("li9_lifetime_dtmax",li9_lifetime_dtmax);
	
	
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
				li9_ntag_dt_vals.push_back(dt_lowe_n[neutron_i]);
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
	
	// this is for candidates that pass the dlt cut
	Log(toolName+" making plots of spallation Dt distributions",v_debug,verbosity);
	PlotSpallationDt();
	
	// mu-lowe time diff for li9 triplet candidates: compare to paper Fig 6
	Log(toolName+" making plot of Li9 candidates Dt distribution",v_debug,verbosity);
	PlotLi9TripletsDt();
	
	Log(toolName+" making plots of Li9 candidate beta spectrum",v_debug,verbosity);
	PlotLi9BetaEnergy();
	
	Log(toolName+" making plot of Li9 ncapture Dt distribution",v_debug,verbosity);
	PlotLi9NtagDt();
	
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

bool PurewaterLi9Plots::PlotMuonDlt(){
	// make histograms of transverse distance to muon for all pre- and post-muons
	// and their difference to extract the spallation distributions
	for(int mu_class_i=0; mu_class_i<6; ++mu_class_i){  // we have 5 muboy classifications
		TH1F ahist_pre(TString::Format("dlt_pre_%d",mu_class_i),
					     "All Pre-Muon to Low-E Transverse Distances",20,0,400);
		for(auto&& aval : dlt_vals_pre.at(mu_class_i)) ahist_pre.Fill(aval);
		ahist_pre.Write();
		
		TH1F ahist_post(TString::Format("dlt_post_%d",mu_class_i),
					     "All Post-Muon to Low-E Transverse Distances",20,0,400);
		for(auto&& aval : dlt_vals_post.at(mu_class_i)) ahist_post.Fill(aval);
		ahist_post.Write();
		
		TH1F* dlt_hist_spall = (TH1F*)ahist_pre.Clone(TString::Format("dlt_spall_%d",mu_class_i));
		dlt_hist_spall->Reset(); // clear entries
		// subtract post from pre
		dlt_hist_spall->Add(&ahist_pre,&ahist_post,1,-1);
		std::cout<<"subtracting "<<ahist_post.GetEntries()<<" post entries from "<<ahist_pre.GetEntries()
				 <<" pre entries results in "<<dlt_hist_spall->GetEntries()<<" ("
				 <<(ahist_post.GetEntries()-ahist_pre.GetEntries())<<") entries"<<std::endl;
		// XXX we seem to lose half our entries here....??
		dlt_hist_spall->Write();
		dlt_hist_spall->SetDirectory(0);
		delete dlt_hist_spall;
	}
	
	return true;
}

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
	
	// write the stack to file so we can easily see all histos
	my_spall_dts.Write("spall_dts_stack");
	
	// also create the paper ones for overlay
	TH1F* ahist = my_spall_dts.GetHists()->At(0);
	PlotPaperDt(ahist);
	PlotPaperDlt(ahist);
	
	// TH1::Clone creates a copy that we own, so are responsible for cleanup.
	// THStacks do not take ownership of their histos, and we can't make them.
	// The way to cleanup is therefore:
	my_spall_dts.GetHists()->Delete();
	
	return true;
}

bool PurewaterLi9Plots::PlotSpallationDt(){
	// dt distribution of muons passing dlt<200cm cut
	TH1F dt_mu_lowe_hist("dt_mu_lowe_hist","Spallation Muon to Low-E Time Differences",500,0,30);
	TH1F dt_mu_lowe_hist_short("dt_mu_lowe_hist_short","Spallation Muon to Low-E Time Differences",500,0,0.25);
	for(auto&& aval : dt_mu_lowe_vals){
		dt_mu_lowe_hist.Fill(fabs(aval));
		dt_mu_lowe_hist_short.Fill(fabs(aval));
	}
	dt_mu_lowe_hist.Write();
	dt_mu_lowe_hist_short.Write();
	
	// XXX fit this with all the rates....
	
	
	return true;
}

bool PurewaterLi9Plots::PlotLi9NtagDt(){
	// events passing dlt cut, li9 energy and lifetime cuts, and with a tagged neutron passing BDT cut:
	TH1F li9_ntag_dt_hist("li9_ntag_dt_hist","Neutron Capture Candidate dt",20,0,550);
	for(auto&& aval : li9_ntag_dt_vals) li9_ntag_dt_hist.Fill(aval); // FIXME weight by num_post_muons and num neutrons
	li9_ntag_dt_hist.Write();
	
	return true;
}

bool PurewaterLi9Plots::PlotLi9BetaEnergy(){
	// events passing dlt cut, li9 energy and lifetime cuts, and with a tagged neutron passing BDT cut:
	TH1F li9_e_hist("li9_e_hist","Li9 Candidate Beta Energy",7,6,li9_endpoint);
	for(auto&& aval : li9_e_vals) li9_e_hist.Fill(aval); // FIXME weight by # post mus & num neutrons
	li9_e_hist.Write();
	
	return true;
}

bool PurewaterLi9Plots::PlotLi9TripletsDt(){
	
	// make a histogram to bin the data
	std::cout<<"making li9 mu-lowe dt histogram"<<std::endl;
	TH1F li9_muon_dt_hist("li9_muon_dt_hist","Muon to Low-E dt for Li9 triplets",
	                       15,li9_lifetime_dtmin,li9_lifetime_dtmax);
	for(auto&& aval : li9_muon_dt_vals){
		li9_muon_dt_hist.Fill(aval);
	}
	std::cout<<"saving to file"<<std::endl;
	li9_muon_dt_hist.Write();
	
	// looks like they've drawn a TGraphErrors: convert our histogram
	std::cout<<"converting to TGraphErrors"<<std::endl;
	int n_bins = li9_muon_dt_hist.GetNbinsX();
	std::vector<float> yvals(n_bins);
	std::vector<float> xvals(n_bins);
	std::vector<float> yerrs(n_bins);
	std::vector<float> xerrs(n_bins);
	for(int bini=0; bini<n_bins; ++bini){
		xvals[bini]=li9_muon_dt_hist.GetBinCenter(bini);
		yvals[bini]=li9_muon_dt_hist.GetBinContent(bini);
		xerrs[bini]=li9_muon_dt_hist.GetBinWidth(1)/2.;
		yerrs[bini]= (yvals[bini]==0) ? 0 : 1./sqrt(yvals[bini]); // FIXME how to handle errors on bin counts 0? 
	}
	TGraphErrors li9_dt_tgraph(n_bins,xvals.data(),yvals.data(),xerrs.data(),yerrs.data());
	li9_dt_tgraph.SetTitle(li9_muon_dt_hist.GetTitle());
	li9_dt_tgraph.SetName("li9_muon_dt_tgraph");
	li9_dt_tgraph.Write();
	
	// do we do chi2 fit or binned likelihood fit? why? FIXME how do we do binned likelihood?
	
	// perform an unbinned chi2 fit
	std::cout<<"doing binned chi2 fit"<<std::endl;
	BinnedLi9DtChi2Fit(&li9_muon_dt_hist);
	
	// perform an unbinned likelihood fit
	// TODO
	
	// perform an unbinned likelihood fit
	std::cout<<"doing unbinned likelihood fit"<<std::endl;
	UnbinnedLi9DtLogLikeFit();
	
	// TODO unbinned extended likelihood fit to also extract number of Li9 events?
	
	return true;
}

bool PurewaterLi9Plots::BinnedLi9DtChi2Fit(TH1F* li9_muon_dt_hist){
	std::cout<<"making TF1 for binned chi2 fit with "<<li9_muon_dt_hist->GetEntries()<<" values"<<std::endl;
	
	// simple chi2 fit of 'y = C + A*exp(-dt/τ)' for all values of C, A, τ
	TF1 li9_muon_dt_func("li9_muon_dt_func","[0]+[1]*exp(-x/[2])",
	                      li9_lifetime_dtmin,li9_lifetime_dtmax);
	
	// set parameter names
	li9_muon_dt_func.SetParName(0,"num bg events");
	li9_muon_dt_func.SetParName(1,"num Li9 events");
	li9_muon_dt_func.SetParName(2,"Li9 lifetime");
	
	// set starting values
	//li9_muon_dt_func.SetParameters(allparams.data());  // pass an array or set individually
	li9_muon_dt_func.SetParameter(0,0);     // TODO estimate from accidental rate of ntag from MC and num events?
	li9_muon_dt_func.SetParameter(1,li9_muon_dt_hist->GetEntries());
	li9_muon_dt_func.SetParameter(2,0.26);  // li9 lifetime 0.26 seconds from paper
	
	// set parameter limits. For minuit at least this is strongly discouraged
	// Maybe ok for chi2 fits but not likelihood fits?
	//li9_muon_dt_func.SetParLimits(4, 0., 200.);
	
	// set num TF1 points for drawing
	li9_muon_dt_func.SetNpx(1000);
	
	// DO THE FIT
	std::cout<<"doing binned chi2 fit"<<std::endl;
	TFitResultPtr fitresult = li9_muon_dt_hist->Fit("li9_muon_dt_func","MRS");
	// options here are  M: better fit, R: use range, S: get resultsptr
	
	// print result
	std::cout<<"li9 dt binned chi2 fit parameters = {";
	for(int i=0; i<(li9_muon_dt_func.GetNpar()); i++){
		std::cout<<li9_muon_dt_func.GetParameter(i);
		(i<(li9_muon_dt_func.GetNpar()-1)) ? std::cout<<", " : std::cout<<"};\n";
	}
	//float fitchi2 = fitresult->Chi2();                               // same as below, which
	float fitchi2 = li9_muon_dt_func.GetChisquare();         // doesn't need fitresultptr
	Log(toolName+" li9 mu->lowe dt fit chi2 was "+toString(fitchi2),v_message,verbosity);
	
	// draw result
//	li9_muon_dt_hist->Draw();
//	//fit->Draw("lsame");   // not necessary, fit is added to histogram's list of functions and drawn automatically
	
	std::cout<<"binned chi2 fit done"<<std::endl;
	return true;
}

bool PurewaterLi9Plots::UnbinnedLi9DtLogLikeFit(){
	
	std::cout<<"doing unbinned likelihood fit with "<<li9_muon_dt_vals.size()<<" values"<<std::endl;
	
	// TF1 of likelihood distribution
	// ROOT provides some basic pdfs: https://root.cern.ch/doc/v610/group__PdfFunc.html
	// but we need our likelihood to be normalized and I don't see how to do it suitably with these
//	TF1 li9_muon_dt_unbinned_like("li9_muon_dt_unbinned_like",
//	TString::Format("[0]*ROOT::Math::uniform_pdf (x,%.2f,%.2f,0) "
//	                " + [1]*ROOT::Math::exponential_pdf(x,[1],0)",xmin,xmax),xmin,xmax);
	
	int fit_n_pars = 2;
	int func_n_dims = 1;
	std::cout<<"making TF1"<<std::endl;
	TF1 li9_muon_dt_unbinned_like("li9_muon_dt_unbinned_like",this,&PurewaterLi9Plots::li9_lifetime_loglike,
	                       li9_lifetime_dtmin, li9_lifetime_dtmax, fit_n_pars,
	                       "PurewaterLi9Plots","li9_lifetime_loglike");
	
	// set parameter names
	li9_muon_dt_unbinned_like.SetParName(0,"fraction of Li9 events");
	li9_muon_dt_unbinned_like.SetParName(1,"Li9 lifetime");
	
	// set starting values
	//li9_muon_dt_unbinned_like.SetParameters(allparams.data());  // pass an array or set individually
	li9_muon_dt_unbinned_like.SetParameter(0,0.8);   // TODO better way to estimate initial val?
	li9_muon_dt_unbinned_like.SetParameter(1,0.26);  // li9 lifetime 0.26 seconds from paper
	
	// convert the TF1 for use in unbinned fit
	std::cout<<"making TF1 wrapper"<<std::endl;
	ROOT::Math::WrappedMultiTF1 li9_mu_dt_unbinned_func(li9_muon_dt_unbinned_like,
	                                                    li9_muon_dt_unbinned_like.GetNdim());
	
	// build a fitter
	// 'false' says let the fitter calculate derivatives of the function itself (recommended)
	std::cout<<"making fitter"<<std::endl;
	ROOT::Fit::Fitter fitter;
	fitter.SetFunction(li9_mu_dt_unbinned_func);  // no bool in ROOT 5 - what's the default?
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
//	ROOT::Fit::UnBinData li9_dt_data(opt, range, x.size());
	
	// convert the data into a suitable 'UnBinData' object
	std::cout<<"making unbinned dataset"<<std::endl;
	ROOT::Fit::UnBinData li9_dt_data(li9_muon_dt_vals.size());
	for(auto aval : li9_muon_dt_vals){   // note: use auto NOT auto&&
		li9_dt_data.Add(aval);   // can we introduce weights here? FIXME
	}
	
	// DO THE FIT
	std::cout<<"doing the fit"<<std::endl;
	fitter.LikelihoodFit(li9_dt_data);
	
	// print the results
	ROOT::Fit::FitResult r=fitter.Result();
	r.Print(std::cout);
	
	std::cout<<"unbinned likelihood fit done"<<std::endl;
	
	// TODO err, extract and return parameters for drawing
	return true;
}

double PurewaterLi9Plots::li9_lifetime_loglike(double* x, double* par){
	
	if(*x<0) return 1e10;  // shouldn't be trying to fit negative times
	
	// like = C + A*exp(-dt/τ)
	//         C = # of accidental backgrounds (constant per unit time)
	//par[0] = A = # li9 decay events
	//par[1] = τ = li9 decay lifetime
	
	// calculate area of expl decay curve over the fitted range for normalization
	double expintegral = par[0]*par[1]*(exp(-li9_lifetime_dtmin/par[1]) - exp(-li9_lifetime_dtmax/par[1]));
	// to normalize this pdf the (constant part * x range) must make up the rest
	double bgfrac = (1. - expintegral)/(li9_lifetime_dtmax-li9_lifetime_dtmin);
	
	// bgfrac must be >0, so if bgfrac<0 return something that rapidly goes to very large values
	if(bgfrac<0) return exp(bgfrac*100.);
	
	// otherwise calculate and return log-likelihood of value x given this set of parameters
	// log(like/integral[like]) = log(like) - log(integral[like])
	return log(bgfrac + par[0]*exp(-(*x)/par[1]));
}



















bool PurewaterLi9Plots::PlotPaperDt(TH1F* basehist){
	TH1F* h_paper_dt_mu_nonfitted = basehist->Clone("paper_dt_nonfitted");
	h_paper_dt_mu_nonfitted->Reset();
	// digitized from paper. I didn't do all series' as most overlapped.
	std::vector<double> paper_dt_mu_nonfitted
		{ 0.5458, 0.3111, 0.1014, 0.0403, 0.0014, 0.0806, 0.0694, 0.0111, 0.0014, 0.0056};
	for(int i=0; i<paper_dt_mu_nonfitted.size(); ++i){
		h_paper_dt_mu_nonfitted.SetBinContent(i, paper_dt_mu_nonfitted[i]);
	}
	h_paper_dt_mu_nonfitted.SetLineColor(kBlue);
	h_paper_dt_mu_nonfitted.SetLineStyle(2);
	
	TH1F* h_paper_dt_mu_single = basehist->Clone("paper_dt_single");
	h_paper_dt_mu_single->Reset();
	std::vector<double> paper_dt_mu_single
		{ 0.5806, 0.2153, 0.0903, 0.0403, 0.0236, 0.0153, 0.0111, 0.0083, 0.0097, 0.0069 };
	for(int i=0; i<paper_dt_mu_single.size(); ++i){
		h_paper_dt_mu_single.SetBinContent(i, paper_dt_mu_single[i]);
	}
	h_paper_dt_mu_single.SetLineColor(kBlack);
	h_paper_dt_mu_single.SetLineStyle(2);
	
	// multigraph for easy combined plotting
	THStack paper_dts;
	paper_dts.Add(&h_paper_dt_mu_nonfitted);
	paper_dts.Add(&h_paper_dt_mu_single);
	
	// write them out to file
	paper_dts.Write("spall_dts_paper");
	
	// clean up
	paper_dts.GetHists()->Delete();
	
	return true;
}

bool PurewaterLi9Plots::PlotPaperDlt(TH1F* basehist){
	// digitized from paper
	std::vector<double> paper_dl_multiple
		{ 0.1414, 0.2303, 0.1980, 0.1465, 0.1121, 0.0798, 0.0596, 0.0404 };
	TH1F* h_paper_dl_multiple = basehist->Clone("paper_dl_multiple");
	h_paper_dl_multiple->Reset();
	for(int i=0; i<paper_dl_multiple.size(); ++i){
		h_paper_dl_multiple.SetBinContent(i,paper_dl_multiple[i]);
	}
	h_paper_dl_multiple.SetLineColor(kRed);
	h_paper_dl_multiple.SetLineStyle(2);
	
	std::vector<double> paper_dl_single
		{ 0.3636, 0.3707, 0.1485, 0.0636, 0.0283, 0.0141, 0.0121, 0.0071 };
	TH1F* h_paper_dl_single = basehist->Clone("paper_dl_single");
	h_paper_dl_single->Reset();
	for(auto&& apair : paper_dl_single){
		h_paper_dl_single.SetPoint(h_paper_dl_single.GetN(), apair.first, apair.second);
	}
	h_paper_dl_single.SetLineColor(kBlack);
	h_paper_dl_single.SetLineStyle(2);
	
	std::vector<double> paper_dl_stopping
		{ 0.4161, 0.3242, 0.0455, 0.0414, 0.0707, 0.0434, 0.0263, 0.0303 };
	TH1F* h_paper_dl_stopping = basehist->Clone("paper_dl_stopping");
	h_paper_dl_stopping->Reset();
	for(int i=0; i<paper_dl_stopping.size(); ++i){
		h_paper_dl_stopping.SetBinContent(i,paper_dl_stopping[i]);
	}
	h_paper_dl_stopping.SetLineColor(kMagenta);
	h_paper_dl_stopping.SetLineStyle(2);
	
	THStack paper_dls;
	paper_dls.Add(&h_paper_dl_multiple);
	paper_dls.Add(&h_paper_dl_single);
	paper_dls.Add(&h_paper_dl_stopping);
	
	// write them out to file
	paper_dls.Write("spall_dts_paper");
	
	paper_dls.GetHists()->Delete();
	
	return true;
}

















































bool PurewaterLi9Plots::PlotPaperDt(){
	// digitized from paper. I didn't do all series' as most overlapped.
	std::vector<std::pair<double,double>> paper_dt_mu_nonfitted{
		{0.012168141592920359,0.5458333333333334},
		{0.03761061946902655,0.3111111111111111},
		{0.062499999999999986,0.10138888888888886},
		{0.08794247787610619,0.040277777777777746},
		{0.11283185840707964,-0.001388888888888884},
		{0.1377212389380531,0.08055555555555549},
		{0.16261061946902655,0.06944444444444442},
		{0.1875,0.011111111111111072},
		{0.21294247787610615,0.001388888888888884},
		{0.2372787610619469,0.005555555555555536}
	};
	TGraph g_paper_dt_mu_nonfitted;
	for(auto&& apair : paper_dt_mu_nonfitted){
		g_paper_dt_mu_nonfitted.SetPoint(g_paper_dt_mu_nonfitted.GetN(), apair.first, apair.second);
	}
	g_paper_dt_mu_nonfitted.SetName("paper_dt_nonfitted");
	g_paper_dt_mu_nonfitted.SetMarkerColor(kBlue);
	g_paper_dt_mu_nonfitted.SetMarkerStyle(2);
	g_paper_dt_mu_nonfitted.SetLineWidth(0);
	std::vector<std::pair<double,double>> paper_dt_mu_single{
		{0.012168141592920359,0.5805555555555555},
		{0.03761061946902655,0.2152777777777778},
		{0.062499999999999986,0.09027777777777779},
		{0.08794247787610619,0.040277777777777746},
		{0.11283185840707964,0.023611111111111138},
		{0.1366150442477876,0.015277777777777724},
		{0.16261061946902655,0.011111111111111072},
		{0.1875,0.008333333333333304},
		{0.21294247787610615,0.009722222222222188},
		{0.23783185840707965,0.00694444444444442}
	};
	TGraph g_paper_dt_mu_single;
	for(auto&& apair : paper_dt_mu_single){
		g_paper_dt_mu_single.SetPoint(g_paper_dt_mu_single.GetN(), apair.first, apair.second);
	}
	g_paper_dt_mu_single.SetName("paper_dt_single");
	g_paper_dt_mu_single.SetMarkerColor(kBlack);
	g_paper_dt_mu_single.SetMarkerStyle(2);
	g_paper_dt_mu_single.SetLineWidth(0);
	
	// multigraph for easy combined plotting
	TMultiGraph paper_dts;
	paper_dts.Add(&g_paper_dt_mu_nonfitted);
	paper_dts.Add(&g_paper_dt_mu_single);
	
	// write them out to file
	paper_dts.Write("spall_dts_paper");
	
	return true;
}

bool PurewaterLi9Plots::PlotPaperDlt(){
	// colour scheme of paper
	std::map<std::string, EColor> class_colours{
		{"misfit",kBlue},
		{"single_thru_going",kBlack},
		{"single_stopping",kMagenta},
		{"multiple_mu",kRed},
		{"also_multiple_mu",kRed},
		{"corner_clipper",kWhite}
	
	// digitized from paper
	std::vector<std::pair<double,double>> paper_dl_multiple{
		{24.724061810154467,0.14141414141414138},
		{75.05518763796908,0.23030303030303023},
		{125.38631346578359,0.19797979797979792},
		{173.95143487858718,0.1464646464646464},
		{225.16556291390725,0.11212121212121207},
		{274.6136865342163,0.07979797979797976},
		{324.06181015452535,0.059595959595959536},
		{374.39293598233985,0.04040404040404033}
	};
	TGraph g_paper_dl_multiple;
	for(auto&& apair : paper_dl_multiple){
		g_paper_dl_multiple.SetPoint(g_paper_dl_multiple.GetN(), apair.first, apair.second);
	}
	g_paper_dl_multiple.SetName("paper_dl_multiple");
	g_paper_dl_multiple.SetMarkerColor(kBlack);
	g_paper_dl_multiple.SetMarkerStyle(2);
	g_paper_dl_multiple.SetLineWidth(0);
	std::vector<std::pair<double,double>> paper_dl_single{
		{23.841059602648897,0.3636363636363636},
		{75.05518763796908,0.37070707070707065},
		{124.50331125827813,0.14848484848484844},
		{174.83443708609263,0.0636363636363636},
		{225.16556291390725,0.028282828282828243},
		{274.6136865342163,0.014141414141414121},
		{324.9448123620308,0.012121212121212088},
		{374.39293598233985,0.007070707070707005}
	};
	TGraph g_paper_dl_single;
	for(auto&& apair : paper_dl_single){
		g_paper_dl_single.SetPoint(g_paper_dl_single.GetN(), apair.first, apair.second);
	}
	g_paper_dl_single.SetName("paper_dl_single");
	g_paper_dl_single.SetMarkerColor(kBlack);
	g_paper_dl_single.SetMarkerStyle(2);
	g_paper_dl_single.SetLineWidth(0);
	std::vector<std::pair<double,double>> paper_dl_stopping{
		24.724061810154467,0.4161616161616161},
		74.17218543046351,0.3242424242424242},
		124.50331125827813,0.045454545454545414},
		174.83443708609263,0.04141414141414135},
		225.16556291390725,0.07070707070707066},
		274.6136865342163,0.04343434343434338},
		323.1788079470198,0.02626262626262621},
		374.39293598233985,0.030303030303030276}
	};
	TGraph g_paper_dl_stopping;
	for(auto&& apair : paper_dl_stopping){
		g_paper_dl_stopping.SetPoint(g_paper_dl_stopping.GetN(), apair.first, apair.second);
	}
	g_paper_dl_stopping.SetName("paper_dl_stopping");
	g_paper_dl_stopping.SetMarkerColor(kBlack);
	g_paper_dl_stopping.SetMarkerStyle(2);
	g_paper_dl_stopping.SetLineWidth(0);
	
	TMultiGraph paper_dls;
	paper_dls.Add(&g_paper_dl_multiple
	paper_dls.Add(&g_paper_dl_single
	paper_dls.Add(&g_paper_dl_stopping);
	
	// write them out to file
	paper_dls.Write("spall_dts_paper");
	
	return true;
}

