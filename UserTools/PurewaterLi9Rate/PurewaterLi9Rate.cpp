#include "PurewaterLi9Rate.h"

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TTree.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TRolke.h"
#include "TMinuit.h"

PurewaterLi9Rate::PurewaterLi9Rate():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

const double li9_endpoint = 14.5; // MeV

bool PurewaterLi9Rate::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	m_data= &data;
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",verbosity);            // how verbose to be
	m_variables.Get("inputFile",inputFile);            // a single specific input file
	m_variables.Get("outputFile",outputFile);          // output file to write
	m_variables.Get("maxEntries",MAX_ENTRIES);         // terminate after processing at most this many events
	m_variables.Get("writeFrequency",WRITE_FREQUENCY); // how many events to TTree::Fill between TTree::Writes
	
	m_variables.Get("max_closest_muon_dt",max_closest_muon_dt); // cut lowe events too close to a mu (pre OR post!)
	m_variables.Get("max_closest_lowe_dx",max_closest_lowe_dx); // cut lowe events too close to another within 60s
	
	m_variables.Get("num_bdt_cut_bins",num_bdt_cut_bins);
	m_variables.Get("bdt_cut_min",bdt_cut_min);
	m_variables.Get("bdt_cut_max",bdt_cut_max);
	m_variables.Get("ntag_FOM_threshold",ntag_FOM_threshold);
	
	m_variables.Get("run_min",run_min);
	m_variables.Get("run_max",run_max);
	
	// open the input TFile and TTree
	// ------------------------------
	get_ok = myTreeReader.Load(inputFile, "data");
	DisableUnusedBranches();
	if(get_ok) ReadEntryNtuple(0);
	
	// create the output TFile and TTree
	// ---------------------------------
	CreateOutputFile(outputFile);
	
	// read ntag ROC curve file and make splines for efficiencies and systematics (functions of ntag cut)
//	fill_ntag_roc(&sig, &bg, &sigsys, &bgsys); // TODO enable when implemented
	// create a vector of BDT FOM cutoffs
//	make_BDT_bins();                           // TODO enable when implemented
	
	// create the histograms to fill
	// -----------------------------
	
	// Spallation and control sample dt and dlt distributions, for each class of muon
	//  - we search for up to 30s prior, but Fig 2 only plots up to 0.25s! XXX
	for(int mu_class_i=0; mu_class_i<5; ++mu_class_i){  // we have 5 muboy classifications
		dlt_hists_pre.emplace(mu_class_i,
							  new TH1F(TString::Format("dlt_hists_pre_%d",mu_class_i),
							           "All Pre-Muon to Low-E Transverse Distances",100,0,400));
		dt_hists_pre.emplace(mu_class_i,
							  new TH1F(TString::Format("dt_hists_pre_%d",mu_class_i),
							           "All Pre-Muon to Low-E Time Differences",100,-30,0));
		dlt_hists_post.emplace(mu_class_i,
							  new TH1F(TString::Format("dlt_hists_post_%d",mu_class_i),
							           "All Post-Muon to Low-E Transverse Distances",100,0,400));
		dt_hists_post.emplace(mu_class_i,
							  new TH1F(TString::Format("dt_hists_post_%d",mu_class_i),
							           "All Post Muon to Low-E Time Differences",100,-30,0));
	}
	
	// all spallation passing dlt cut - we search for up to 30s prior, but Fig 2 only plots up to 0.25s! XXX
	dt_mu_lowe_hist = new TH1F("dt_mu_lowe_hist","Spallation Muon to Low-E Time Differences",100,0,30);
	
	// events passing dlt cut, li9 energy and lifetime cuts, and with a tagged neutron passing BDT cut:
	li9_muon_dt_hist = new TH1F("li9_muon_dt_hist","Muon to Low-E dt for Li9 triplets",100,-30,0);
	li9_ntag_dt_hist = new TH1F("li9_ntag_dt_hist","Neutron Capture Candidate dt",100,0,550e-6);
	li9_e_hist = new TH1F("li9_e_hist","Li9 Candidate Beta Energy",100,6,li9_endpoint);
	
	// Pre-populate event count tracker with all cuts to get the ordering of reduction right.
	// FIXME find a better way. Maybe we should split lowe events // mu-lowe pairs // mu-lowe-ntag triplets,
	// as otherwise it need not be monotonic. Easiest way for now is:
	// `grep IncrementEventCount UserTools/PurewaterLi9Rate/PurewaterLi9Rate.cpp`
	// TODO build cut names to reflect actual cuts
	cut_tracker.emplace("all",0);
	cut_tracker.emplace("68671<run<73031",0);
	cut_tracker.emplace("SNR>0.5",0);
	cut_tracker.emplace("dwall>200cm",0);
	cut_tracker.emplace("dt_mu_lowe>50us",0);
//	cut_tracker.emplace("thirdred",0);
//	cut_tracker.emplace("max_hits_200ns_AFT<50",0);
//	cut_tracker.emplace("min(dt_mu_lowe)>1ms",0);
//	cut_tracker.emplace("nearest_other_lowe>490cm",0);
	cut_tracker.emplace("lowe_energy>6MeV",0);
	cut_tracker.emplace("mu_lowe_pairs",0);
	cut_tracker.emplace("muboy_index==0",0);
	cut_tracker.emplace("dlt_mu_lowe>200cm",0);
	cut_tracker.emplace("lowe_energy_in_li9_range",0);
	cut_tracker.emplace("dt_mu_lowe_in_li9_range",0);
	cut_tracker.emplace("lowest_other_mu_dt>1ms",0);
	cut_tracker.emplace("ntag_FOM>0.995",0);
	cut_tracker.emplace("mu_lowe_ntag_triplets",0);
	
	return true;
}

bool PurewaterLi9Rate::Execute(){
	
	Log(toolName+" processing entry "+toString(entry_number),v_debug,verbosity);
	
	// clear output vectors so we don't carry anything over
	Log(toolName+" clearing output vectors",v_debug,verbosity);
	ClearOutputTreeBranches();
	
	// Do Li9 analysis
	Log(toolName+" doing analysis",v_debug,verbosity);
	Analyse();
	
	// Fill the output tree
	Log(toolName+" filling output TTree entry",v_debug,verbosity);
	outtree->Fill();
	
	// update the output file so we don't lose everything if we crash
	if((entry_number%WRITE_FREQUENCY)==0) WriteTree();
	
	// stop at user-defined limit to the number of events to process
	++entry_number;
	if((MAX_ENTRIES>0)&&(entry_number>=MAX_ENTRIES)){
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

bool PurewaterLi9Rate::Finalise(){
	
	Log(toolName+" event counts trace: ",v_warning,verbosity);
	for(auto&& event_count : cut_tracker){
		std::cout<<event_count.first<<" : "<<event_count.second<<"\n";
	}
	
	// subtract distribution of post-muons from pre-muons to obtain lt and dt distributions
//	PlotSpallationDtDlt(); // using dlt_hists_pre and dlt_hists_post  TODO
	
	// measure dlt cut systematic TODO
//	MeasureDltSystematic(); // using variation in dlt_systematic_dt_cuts in bin corresponding to dlt=200cm
	// for each entry in dlt_systematic_dt_cuts_pre, dlt_systematic_dt_cuts_post, subtract post from pre.
	// this gives a set of efficiencies of spallation for varying dts.
	// compare efficiency across various dts in each dlt bin to get the dlt systematic error.
	
	// fit time distribution of muon-lowe pairs with the nominal dlt cut to obtain isotope counts
//	MeasureIsotopeRates();  // TODO
	// split into sum of contributions from each isotope by doing piecewise fit to time distribution:
	// F(t) = sum_i (N_i/τ_i) exp(-dt/τ_i) + const
	
	// ensure everything is written to the output file
	// -----------------------------------------------
	get_ok = WriteTree();
	if(not get_ok){
		Log(toolName+" Error writing output TTree!",v_error,verbosity);
	}
	
	// Close and delete the file handle
	// --------------------------------
	CloseFile();
	
	return true;
}

// #####################################################################

// main body of the tool
bool PurewaterLi9Rate::Analyse(){
	// This gets called for each Execute iteration, to process one lowe event
	
	Log(toolName+" entry "+toString(entry_number)+" run "+toString(HEADER->nrunsk),v_debug,verbosity);
	
	IncrementEventCount("all");
	if(HEADER->nrunsk < run_min) return false;    // SHE trigger threshold lowered to 8 MeV
	if(HEADER->nrunsk > run_max){                // Yang Zhang's time range
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	//if (HEADER->nrunsk > 74781) return false;  // WIT started after this run. What's the significance of this?
	Log(toolName+" entry "+toString(entry_number)+" run "+toString(HEADER->nrunsk)+" ok",v_debug,verbosity);
	IncrementEventCount("68671<run<73031");
	
	// find lowe events ✅
	// for reference, paper says 54,963 beta events... though not clear after which cuts
	
	// n_hits_with_Q_lt_0.5pe / n_hits_total > 0.55
	if((double)thirdredvars->q50 / (double)(LOWE->bsn50) > 2.) return false;   // is this the right cut? XXX
	IncrementEventCount("SNR>0.5");
	
	// dwall > 2m
	if(thirdredvars->dwall < 200.) return false;
	IncrementEventCount("dwall>200cm");
	
	// dt_muon_lowe > 50us
	// check the closest preceding muon. 
	// n.b. if muboy found multiple, they all have the same time, so we don't neeed to scan
	if(dt_mu_lowe[num_pre_muons-1] < 50e-6) return false;
	IncrementEventCount("dt_mu_lowe>50us");
	
	Log(toolName+" entry "+toString(entry_number)+" passed 1st reduction",v_debug,verbosity);
	
//	// new set of cuts from atmospheric analysis TODO enable?
//	if (apply_third_reduction(third, LOWE)) return false;
//	Log(toolName+" entry "+toString(entry_number)+" passed 3rd reduction",v_debug,verbosity);
//	IncrementEventCount("thirdred");
	
//	// new cut of events with a low energy muon in the ncapture window TODO enable?
//	if (max_hits_200ns_AFT > 50) return false;
//	Log(toolName+" entry "+toString(entry_number)+" passed muon in AFT cut",v_debug,verbosity);
//	IncrementEventCount("max_hits_200ns_AFT<50");
	
//	// new cut of events within 1ms of nearest muon (pre- or post) cut TODO enable?
//	// how would this work? We're looking for spallation: we want events in proximity to a muon...
//	bool fail_1ms = false;
//	for ( int i = 0; i < num_pre_muons; i++ ) {
//		if ( fabs( dt_mu_lowe[ i ] ) < 0.001 ){
//			fail_1ms = true;
//			break;
//		}
//	}
//	if(fail_1ms) return false;
//	Log(toolName+" entry "+toString(entry_number)+" passed >1ms to nearest pre-muon cut",v_debug,verbosity);
//	IncrementEventCount("min(dt_mu_lowe)>1ms");
	
//	// new cut to remove all lowe events in close proximity to another lowe event TODO enable?
//	if (closest_lowe_60s[0] < max_closest_lowe_dx) return false;    // not applied by Zhang 
//	Log(toolName+" entry "+toString(entry_number)+" passed >60s to nearest lowe event cut",v_debug,verbosity);
//	IncrementEventCount("nearest_other_lowe>490cm");
	
	// cut all lowe events with energy < 6 MeV (mostly non-spallation)
	if( LOWE->bsenergy < 6.f ) return false;
	Log(toolName+" entry "+toString(entry_number)+" passed LOWE energy >6 MeV cut",v_debug,verbosity);
	IncrementEventCount("lowe_energy>6MeV");
	
	// find muons within 30 prior to lowe events ( muons identified by: >1000pe in ID ) ✅
	// find muons within 30s post lowe events    ( muons identified by: >1000pe in ID ) ✅
	// pair all lowe events with all pre-muons ✅
	// pair all lowe events with all post-muons ✅
	
	// plot lt, dt distributions as a function of muon class, for both pre- and post-muons
	Log(toolName+" looping over "+toString(num_pre_muons)+" preceding muons and "+toString(num_post_muons)
				+" following muons to fill spallation and control dl, dt histograms",v_debug,verbosity);
	// pre muons
	for(int mu_i=0; mu_i<num_pre_muons; ++mu_i){
		// only consider first muboy muon (only for multi-mu events?)
		if(mu_index[mu_i] > 0) return false;
		
		dlt_hists_pre.at(mu_class[mu_i])->Fill(dlt_mu_lowe[mu_i]);   // FIXME weight by num_pre_muons
		dt_hists_pre.at(mu_class[mu_i])->Fill(dt_mu_lowe[mu_i]);     // FIXME weight by num_pre_muons
		
		// to evaluate systematic on lt cut, apply various dt cuts and see how the lt cut efficiency varies
		// since we're interested in the effect on the spallation sample, which is given by
		// the total - post-muon sample, record both pre- and post- muon samples with various dt cuts
		for(auto&& a_dt_cut : dlt_systematic_dt_cuts_pre){
			if(dt_mu_lowe[mu_i] < a_dt_cut.first) a_dt_cut.second->Fill(dt_mu_lowe[mu_i]);
		}
	}
	// post muons
	for(int mu_i=num_pre_muons; mu_i<(num_pre_muons+num_post_muons); ++mu_i){
		if (mu_index[mu_i] > 0) continue;
		
		dlt_hists_post.at(mu_class[mu_i])->Fill(dlt_mu_lowe[mu_i]);  // FIXME weight by num_post_muons
		dt_hists_post.at(mu_class[mu_i])->Fill(dt_mu_lowe[mu_i]);    // FIXME weight by num_post_muons
		
		for(auto&& a_dt_cut : dlt_systematic_dt_cuts_post){
			if(dt_mu_lowe[mu_i] < a_dt_cut.first) a_dt_cut.second->Fill(dt_mu_lowe[mu_i]);
		}
	}
	// in Finalise we'll substract the two to get dt and dlt distributions for spallation only.
	// we'll also compare across various dt cuts to get the systematic error on the spallation dlt cut.
	
	// the following cuts are based on muon-lowe pair variables, so loop over muon-lowe pairs
	Log(toolName+" Looping over "+toString(num_pre_muons)
				+" preceding muons to look for spallation events",v_debug,verbosity);
	for(int mu_i=0; mu_i<num_pre_muons; ++mu_i){
		IncrementEventCount("mu_lowe_pairs");
		
		// safety check: should not consider muons after lowe event
		if (dt_mu_lowe[mu_i] >= 0) break;
		IncrementEventCount("muboy_index==0");
		
		// Apply nominal lt cut, unless muon type was misfit or a poorly fit single muon
		if(not (dlt_mu_lowe[mu_i] < 200 || mu_class[mu_i] == muboy_classes::misfit || 
			   (mu_class[mu_i] == muboy_classes::single_thru_going && mu_fit_goodness[mu_i] < 0.4))) continue;
		Log(toolName+" entry "+toString(entry_number)+", muon "+toString(mu_i)
					+" passed dlt cut",v_debug,verbosity);
		IncrementEventCount("dlt_mu_lowe>200cm");
		
		// record the distribution of dt_mu_lowe
		dt_mu_lowe_hist->Fill(dt_mu_lowe[mu_i]);  // FIXME weight by num_pre_muons
		// in finalise we'll fit this distribution to estimate the number of events from each isotope
		
		// That's all for assessing the amount of general spallation isotopes
		// --------------------------------------------------------------
		// for Li9 we have a number of additional cuts
		
		// apply Li9 energy range cut
		if ( LOWE->bsenergy <= 7.5 || LOWE->bsenergy >= li9_endpoint ) continue;
		Log(toolName+" entry "+toString(entry_number)+", muon "+toString(mu_i)
					+" passed Li9 energy range cut",v_debug,verbosity);
		IncrementEventCount("lowe_energy_in_li9_range");
		
		// apply Li9 lifetime cut
		if( dt_mu_lowe[mu_i] < 0.05 || dt_mu_lowe[mu_i] > 0.5 ) continue;
		Log(toolName+" entry "+toString(entry_number)+", muon "+toString(mu_i)
					+" passed Li9 lifetime cut",v_debug,verbosity);
		IncrementEventCount("dt_mu_lowe_in_li9_range");
		
		// no other mu within 1ms of this lowe event
		bool other_muon_within_1ms = false;  // XXX check we're interpreting this cut right
		for(int othermu_i=0; othermu_i<(num_pre_muons+num_post_muons); ++othermu_i){
			if(othermu_i==mu_i) continue; // looking for muons other than the current one
			if(fabs(dt_mu_lowe[othermu_i])<max_closest_muon_dt){
				other_muon_within_1ms = true;
				break;
			}
		}
		if(other_muon_within_1ms) continue;
		Log(toolName+" entry "+toString(entry_number)+", muon "+toString(mu_i)
					+" passed >1ms to another muon cut",v_debug,verbosity);
		IncrementEventCount("lowest_other_mu_dt>1ms");
		
		// search for ncapture candidates: >7 hits within 10ns T-TOF in 50ns-535us after lowe events ✅
		
		// calculate neutron FOM and cut failing ones
		if( (num_neutron_candidates==0) || 
			(*std::max_element(ntag_FOM.begin(), ntag_FOM.end())<ntag_FOM_threshold)) continue;
		Log(toolName+" entry "+toString(entry_number)+", muon "+toString(mu_i)
					+" passed ntag cut",v_debug,verbosity);
		// XXX as reference, we should have 116 remaining candidate events here
		IncrementEventCount("ntag_FOM>0.995");
		
		// apparently in Zhang study no events had multiple ntag candidates
		if(num_neutron_candidates>1){
			Log(toolName+" event with "+toString(num_neutron_candidates)
				+" neutron candidates!",v_debug,verbosity);
		}
		
		// plot distribution of beta energies from passing triplets, compare to fig 4
		li9_e_hist->Fill(LOWE->bsenergy); // FIXME weight by num_post_muons
		
		// plot distirbution of mu->beta   dt from passing triplets, compare to fig 6
		li9_muon_dt_hist->Fill(dt_mu_lowe[mu_i]); // FIXME weight by num_post_muons
		
		// plot distribution of beta->ntag dt from passing triplets, compare to fig 5
		// Zhang had no events with >1 ntag candidate: should we only take the first? XXX
		for(int neutron_i=0; neutron_i<num_neutron_candidates; ++neutron_i){
			li9_ntag_dt_hist->Fill(dt_lowe_n[neutron_i]); // FIXME weight by num_post_muons and num neutrons
			IncrementEventCount("mu_lowe_ntag_triplets");
		}
	
	} // end loop over muons
	
	return true;
}

bool PurewaterLi9Rate::SoniasAnalyse(){
	// neatened version of sonia's analyse code. Doesn't reproduce all histograms,
	// but reproduces the general series of cuts. Bugs included. ;)
	
	//if (HEADER->nrunsk > 73031) return false;             // Yang Zhang's time range
	//if (HEADER->nrunsk < 68671) return false;             // SHE trigger threshold lowered to 8 MeV
	
	if (max_hits_200ns_AFT > 50) return false;              // Remove low energy muons in AFT trigger
	if (apply_third_reduction(thirdredvars, LOWE)) return false;   // Apply third reduction cuts
	
	// remove all lowe events within 1ms of a (pre or post) muon
	bool fail_1ms = false;
	for ( int i = 0; i < num_pre_muons; i++ ) {
		if ( fabs( dt_mu_lowe[ i ] ) < 0.001 ){
			fail_1ms = true;
			break;
		}
	}
	if(fail_1ms) return false;
	
	// remove all lowe events in close proximity to another lowe event
	//if (closest_lowe_60s[0] < cut_multi) return false;    // not applied by Zhang
	
	// loop over ncapture candidates for this lowe event
	for (Int_t j=0; j<num_neutron_candidates; j++) {
		
		// Remove 10mus time shift for AFT events and convert lowe->neutron time to microseconds
		float dtt = dt_lowe_n[j] < 50000 ? dt_lowe_n[j]/1000 : dt_lowe_n[j]/1000 - 65;  // XXX ???
		
		// truncate BDT FOM to histogrammed range
		float neutj = ntag_FOM[j];
		if (ntag_FOM[j] > bdt_cut_thresholds[num_bdt_cut_bins])
			neutj = (bdt_cut_thresholds[num_bdt_cut_bins]+bdt_cut_thresholds[num_bdt_cut_bins - 1])/2.;
		if (ntag_FOM[j] < bdt_cut_thresholds[0])
			neutj = (bdt_cut_thresholds[0]+bdt_cut_thresholds[1])/2.;
		
		// cut ncapture events <18us from lowe event "Same range as for MC studies" (??)
		if (dt_lowe_n[j] < 18000) continue;
		
		// "preselection" ...?
		if (ntag_FOM[j] < -1) continue;
		
		// the following 2 loops are basically the same, with the following differences:
		// the first:
		//      uses a fixed BDT FOM cut (ncut_li9)
		//      excludes misfit muons from lt cut
		//      applies an energy < li9 endpoint cut
		//      stores ncapture times in the array dtimes_li9 (ignoring muon pairings)
		// the second:
		//      scans various BDT FOM cuts (x-axis of ndt histogram)
		//      does not exclude misfit muons from lt cut
		//      does not apply li9 endpoint cut
		//      stores ncapture times in the 2D histogram li9dt (ignoring muon pairings)
		//
		// both of these loop over all muon pairings, but each time record dtt,
		// which is the neutron capture time, which will be the same for each loop???
		// only nspatot doesn't get double counted ... even in the full code
		
		
		// loop over muon pairings for this lowe event
		bool first_passing_muon = true;
		for ( int i = 0; i < num_pre_muons; i++ ) {
			
			// do not consider muons after the lowe event
			if (dt_mu_lowe[i] >= 0)  break;
			
			// only consider first muboy muon (only for multi-mu events?)
			if (mu_index[i] > 0) continue;
			
			// Apply lt cut, unless muon type was misfit (0), or a poorly fit single muon
			if(not (dlt_mu_lowe[i] < 200 || mu_class[i] == 0 ||
				   (mu_class[i] == 1 && mu_fit_goodness[i] < 0.4))) continue;
			
			// apply Li9 lifetime cut
			if(fabs(dt_mu_lowe[i]) <= 0.05 || fabs(dt_mu_lowe[i]) >= 0.5) continue;
			
			// apply beta energy < Li9 endpoint cut
			if( LOWE->bsenergy >= li9_endpoint ) continue;
			
			// track how many lowe events pass li9 candidate selection before ntag cut
			// only increment for the first muon pairing of this lowe event to avoid double-counting
			if (first_passing_muon){ nspatot++; first_passing_muon = false; }
			
			// apply BDT FOM cut
			if (neutj <= ntag_FOM_threshold ) continue;
			
			// remaining events are Li9 candidates! note time from mu->lowe event (li9 decay time)
			dtimes_li9_ncap.push_back(dtt);
			
		}
		
		// likelihood fit of dt_neutron, similar to Yang Zhang 
		int max_bdtbin = 0; //ndt->GetXaxis()->FindBin(neutj); // Bin corresponding to the ntag BDT value TODO
		// we explore a range of BDT cutoffs
		// add this dt to all collections for which the BDT FOM passes the cutoff
		for(int q = 0; q <= max_bdtbin; q++){
			
			// add this lowe-neutron dt to the distribution for this bdt cut threshold
			//Float_t neutcenter = ndt->GetXaxis()->GetBinCenter(q);  TODO implement
			//ndt->Fill(neutcenter,dtt);   // all lowe->ntag times, no cuts, no muon consideration
			
			// Loop over muon pairings for this lowe event
			for ( int i = 0; i < num_pre_muons; i++ ) {
				
				// do not consider muons after the lowe event
				if (dt_mu_lowe[i] >= 0) break;
				
				// only consider first muboy muon (only for multi-mu events?)
				if (mu_index[i] > 0) continue;
				
				// Apply lt cut
				if(dlt_mu_lowe[i] >= 200) continue;  // *** no special treatment for poorly fit muons ***
				
				// apply 9Li lifetime cut
				if(fabs(dt_mu_lowe[i]) <= 0.05 || fabs(dt_mu_lowe[i]) >= 0.5) continue;
				
				// apply beta energy < Li9 endpoint cut
				//if( LOWE->bsenergy >= li9_endpoint ) continue;   // ??? not applied here?
				
				// save neutron capture time
				//li9dt->Fill(neutcenter,dtt);    // poorly named histogram! TODO implement
				
			}   // loop over muon pairings
		}       // loop over ntag thresholds
	}           // loop over neutron candidates
	
	return true;
}

bool PurewaterLi9Rate::apply_third_reduction(const ThirdRed *th, const LoweInfo *LOWE){
	Log(toolName+" applying third reduction",v_debug,verbosity);
	if (th->maxpre >= 12)             return true;                // remove events with high pre-activity
	if (th->nmue > 0)                 return true;                // remove events ...with a muon...? how? XXX
	if (th->dwall < 200.)             return true;                // remove events within 2m from wall
	if (th->effwall < 500.)           return true;                // remove events with back projection <5m to wall
	if (th->pilike >= 0.36)           return true;                // remove pion-like events
	if ((double)th->q50 / (double)(LOWE->bsn50) > 2) return true; // remove events with low signal-to-noise
	if (th->angle<38 || th->angle>50) return true;                // remove events with Cherenkov angle != 42°
	if (LOWE->bsenergy<8)             return true;                // remove events with energy < 8 MeV
	if ( LOWE->bsgood[1] < 0.5 )      return true;                // remove events with ... poor vertex fit? XXX
	if (th->ovaq < 0.25 )             return true;                // remove events poor event quality
	return false;
}

//// Read in the 'receiver operating characteristic' (ROC) curves for the BDT.
//void PurewaterLi9Rate::fill_ntag_roc(TSpline3 **sig, TSpline3 **bg, TSpline3 **sigsys, TSpline3 **bgsys){
//	// these specify signal acceptance vs background acceptance over a span of FOM cutoff values.
//	// Fit a spline to the ROC so that we can determine the efficiency and background acceptances
//	// for an arbitrary FOM cutoff.
//	Log(toolName+" reading BDT ROC from file "+bdt_outfile,v_message,verbosity);
//	ifstream fin(bdt_outfile);
//	const int splmax = 5000;
//	double cut[splmax], sg[splmax], bkg[splmax], sgsys[splmax], bkgsys[splmax], dum1[splmax], dum2[splmax];
//	int k = 0;
//	while(fin){
//		fin >> cut[k] >> sg[k] >> bkg[k] >> sgsys[k] >> bkgsys[k] >> dum1[k] >> dum2[k];
//		if (cut[k] > bdt_cut_lowth) break; // XXX what is bdt_cut_lowth??? 0.999 in cut_third.h
//		k++;
//	}
//	Log(toolName+" creating spline fit to "+toString(k)+" BDT threshold values",v_message,verbosity);
//	*sig = new TSpline3("signal", cut, sg, k, 0, 1, 0);
//	*bg = new TSpline3("signal", cut, bkg, k, 0, 1, 0);
//	*sigsys = new TSpline3("signal", cut, sgsys, k, 0, 0, 0);
//	*bgsys = new TSpline3("signal", cut, bkgsys, k, 0, 0, 0);
//}

// make an array of BDT FOM thresholds. Make (num_bdt_cut_bins+1) thresholds between
// bdt_cut_min and bdt_cut_max, corresponding to (num_bdt_cut_bins) bins in a histogram.
void PurewaterLi9Rate::make_BDT_bins(){
	Log(toolName+" building array of "+toString(num_bdt_cut_bins)+" BDT cutoff values",v_debug,verbosity);
	bdt_cut_thresholds.resize(num_bdt_cut_bins);
	// convert to log scale
	int cutmin = log(1 - bdt_cut_min)/log(10) - 0.5;   // For cut edges (dec log)
	int cutmax = log(1 - bdt_cut_max)/log(10) + 0.5; // For cut edges (dec log)
	
	for(int i = num_bdt_cut_bins; i >=0; i--){
		float cutval = cutmin + i * 1./num_bdt_cut_bins * (cutmax - cutmin);
		bdt_cut_thresholds[num_bdt_cut_bins - i] = 1 - pow(10, cutval);
	}
	printVals(bdt_cut_thresholds, v_debug, verbosity, "BDT cutoff values: ");
}

void PurewaterLi9Rate::IncrementEventCount(std::string cutname){
	// track remaining numbers of events after each cut
	if(cut_tracker.count(cutname)==0){
		cut_tracker.emplace(cutname,1);
	} else {
		cut_tracker.at(cutname)++;
	}
}

// #####################################################################

int PurewaterLi9Rate::ReadEntryNtuple(long entry_number){
	// get next entry from TreeReader
	int bytesread = myTreeReader.GetEntry(entry_number);
	if(bytesread<=0) return bytesread;
	
	// retrieve variables from branches
	// for efficiency, add all used branches to DisableUnusedBranches
	int success = 
	(myTreeReader.GetBranchValue("HEADER", HEADER)) &&
	(myTreeReader.GetBranchValue("LOWE", LOWE)) &&
	(myTreeReader.GetBranchValue("ThirdRed", thirdredvars)) &&
	(myTreeReader.GetBranchValue("np", num_neutron_candidates)) &&
	(myTreeReader.GetBranchValue("N200M", max_hits_200ns_AFT)) &&
	(myTreeReader.GetBranchValue("neutron5", ntag_FOM)) &&
	(myTreeReader.GetBranchValue("dt", dt_lowe_n)) &&
	(myTreeReader.GetBranchValue("nmusave_pre", num_pre_muons)) &&
	(myTreeReader.GetBranchValue("nmusave_post", num_post_muons)) &&
	(myTreeReader.GetBranchValue("mubstatus", mu_class)) &&
	(myTreeReader.GetBranchValue("mubitrack", mu_index)) &&
	(myTreeReader.GetBranchValue("mubgood", mu_fit_goodness)) &&
	(myTreeReader.GetBranchValue("spadt", dt_mu_lowe)) &&
	(myTreeReader.GetBranchValue("spadlt", dlt_mu_lowe)) &&
	(myTreeReader.GetBranchValue("multispa_dist", closest_lowe_60s));
	
	return success;
}

int PurewaterLi9Rate::DisableUnusedBranches(){
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
	
	return myTreeReader.OnlyEnableBranches(used_branches);
}

int PurewaterLi9Rate::CreateOutputFile(std::string outputFile){
	Log(toolName+": Creating output file "+outputFile,v_message,verbosity);
	outfile = new TFile(outputFile.c_str(), "RECREATE");
	outtree = new TTree("eventtree", "Li9 Spallation Analysis");
	
	// make output branches
	// XXX FILL ME
	
	return 0;
}

void PurewaterLi9Rate::ClearOutputTreeBranches(){
	// XXX FILL ME
}

int PurewaterLi9Rate::WriteTree(){
	Log(toolName+" writing TTree",v_debug,verbosity);
	outfile->cd();
	// TObject::Write returns the total number of bytes written to the file.
	// It returns 0 if the object cannot be written.
	int bytes = outtree->Write("",TObject::kOverwrite);
	if(bytes<=0){
		Log(toolName+" Error writing TTree!",v_error,verbosity);
	} else if(verbosity>2){
		Log(toolName+ " Wrote "+toString(get_ok)+" bytes",v_debug,verbosity);
	}
	return bytes;
};

void PurewaterLi9Rate::CloseFile(){
	outtree->ResetBranchAddresses();
	outfile->Write("*",TObject::kOverwrite);
	outfile->Close();
	delete outfile;
	outfile = nullptr;
};

