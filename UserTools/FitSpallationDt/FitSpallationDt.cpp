/* vim:set noexpandtab tabstop=4 wrap */
#include "FitSpallationDt.h"

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"

#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include "TH1.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

FitSpallationDt::FitSpallationDt():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

const double fiducial_vol = 22.5;               // kton, from paper
const double paper_livetime = 1890;             // days, from paper
const double paper_first_reduction_eff = 99.;   // fraction, from paper
const double paper_dlt_cut_eff = 78.8;          // fraction, from paper

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

const std::map<std::string,double> papervals{  // from 2015 paper table II, Rate in events/kton/day
	// NOTE: these have been corrected for the efficiency of 6MeV energy threshold,
	// so are the rate of ALL DECAYS, not just the rate of decays that produce betas with >6MeV!
	{"11Be",5.7},   // read off plot, 16.9 upper limit in table
	{"16N",39.7},
	{"15C",3.5},    // read off plot, 6.7 upper limit in table
	{"8Li",3.9},    // read off plot, combined in table but plot splits them
	{"8B",4.9},     // read off plot, combined in table but plot splits them
//	{"8Li_8B",8.3},
	{"16C",0},
	{"9Li",0.9},
	//{"9C",0.1},   // read off plot, 0.7 upper limit in table
	//{"8He",0.2},  // read off plot, 0.7 upper limit in table (combined in table but plotted separately?)
	{"8He_9C",0.3}, // read off plot, 1.4 upper limit
	{"12Be",0},
	{"12B",19.8},
	{"13B",0},
	{"14B",0},
	{"12N",2.8},
	{"13O",0},
	{"11Li",0},
	{"const",330}   // not in table, use value from plot, uhhhhh units? Events/0.006s/FV?
};

const std::map<std::string,double> papereffs{
	{"11Be",48.84},
	{"16N",57.68},
	{"15C",40.76},
	{"8Li",54.86},
	{"8B",65.76},
	{"16C",0},
	{"9Li",50.25},
	{"9C",64.35},
	{"8He",28.84},
	{"12Be",0},
	{"12B",58.32},
	{"13B",0},
	{"14B",0},
	{"12N",72.04},
	{"13O",0},
	{"11Li",0}
};

bool FitSpallationDt::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",verbosity);            // how verbose to be
	m_variables.Get("outputFile",outputFile);          // where to save data. If empty, current TFile
	// hacky way to scale the total number of events when comparing to the paper plots
	m_variables.Get("paper_scaling",paper_scaling);
	
	// normally this Tool obtains its data from PurewaterLi9Plots, which fills the vector
	// of spallation muon to lowe event time differences (dt_mu_lowe_vals)
	// and places a pointer to it in the CStore. However we also support filling that
	// vector from spallation files directly
	m_variables.Get("inputFile",inputFile);            // input spallation files for dt values
	m_variables.Get("treeName",treeName);              // name of input tree in spallation files
	m_variables.Get("maxEvents",maxEvents);            // user limit to number of events to process
	m_variables.Get("livetime",livetime);              // if not provided by upstream tool
	
	// energy threshold efficiencies, from FLUKA
	m_variables.Get("efficienciesFile",efficienciesFile);
	
	// read efficiencies of the different thresholds of old vs new data
	GetEnergyCutEfficiencies();
	
	// quick hack to load data from spallation files by laura for validation
	if(inputFile!=""){
		get_ok = myTreeReader.Load(inputFile, treeName);
		if(not get_ok){
			Log(toolName+" failed to open reader on tree "+treeName+" in file "+inputFile,v_error,verbosity);
			return false;
		}
		DisableUnusedBranches();  // for efficiency of reading, only enable used branches
		
		// normally populated by upstream tools
		m_data->CStore.Set("livetime",livetime);
		m_data->CStore.Set("dt_mu_lowe_vals",&my_dt_mu_lowe_vals);
	}
	
	return true;
}


bool FitSpallationDt::Execute(){
	
	Log(toolName+" getting entry "+toString(entrynum),v_debug,verbosity);
	
	// retrieve desired branches
	get_ok = GetBranches();
	
	// process the data
	try{
		Analyse();
	}
	catch(std::exception& e){
		// catch any exceptions to ensure we always increment the event number
		// and load the next entry. This prevents us getting stuck in a loop
		// forever processing the same broken entry!
		Log(toolName+" encountered error "+e.what()+" during Analyse()",v_error,verbosity);
	}
	
	// move to next entry
	entrynum++;
	// check if we've hit the user-requested entry limit
	if((maxEvents>0)&&(entrynum==maxEvents)){
		Log(toolName+" hit max events, setting StopLoop",v_error,verbosity);
		m_data->vars.Set("StopLoop",1);
		return 1;
	}
	
	// pre-load the next ttree entry
	get_ok = ReadEntry(entrynum);
	if(get_ok==0){
		return 1; // end of file
	} else if (get_ok<0){
		return 0; // read error
	}
	
	return true;
}

bool FitSpallationDt::Analyse(){
	// we do our analysis (fitting the distribution) in finalise
	// all we need to do here is add the next dt value
	my_dt_mu_lowe_vals.push_back(dt);
	return true;
}

int FitSpallationDt::ReadEntry(long entry_number){
	// load next entry data from TTree
	int bytesread = myTreeReader.GetEntry(entry_number);
	
	// stop loop if we ran off the end of the tree
	if(bytesread==0){
		Log(toolName+" hit end of input file, stopping loop",v_message,verbosity);
		m_data->vars.Set("StopLoop",1);
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

int FitSpallationDt::GetBranches(){
	int success = (
		(myTreeReader.GetBranchValue("dt",dt))
	);
	return success;
}

int FitSpallationDt::DisableUnusedBranches(){
	std::vector<std::string> used_branches{
		// list used branches here
		"dt"
	};
	return myTreeReader.OnlyEnableBranches(used_branches);
}


bool FitSpallationDt::Finalise(){
	
	// make a new file if given a filename, or if blank check there is a valid file open
	TFile* fout = nullptr;
	if(outputFile!=""){
		fout = new TFile(outputFile.c_str(),"RECREATE");
		fout->cd();
	} else {
		if(gDirectory->GetFile()==nullptr){
			Log(toolName+" Error! No output file given and no file open!",v_error,verbosity);
			return true;
		}
	}
	
	Log(toolName+" Fitting spallation dt distribution",v_debug,verbosity);
	PlotSpallationDt();
	
	if(fout!=nullptr){
		fout->Close();
		delete fout;
		fout=nullptr;
	}
	
	return true;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// =====================================================================

bool FitSpallationDt::PlotSpallationDt(){
	
	// get the data from the CStore
	std::vector<float>* dt_mu_lowe_vals; // input data
	get_ok = m_data->CStore.Get("dt_mu_lowe_vals",dt_mu_lowe_vals);
	if(not get_ok){
		Log(toolName+" Error! No dt_mu_lowe_vals in CStore!",v_error,verbosity);
		return true;
	}
	std::cout<<"fitting "<<dt_mu_lowe_vals->size()<<" spallation dt values"<<std::endl;
	
	// we also need the livetime to convert number of events to rates
	get_ok = m_data->CStore.Get("livetime",livetime);
	if(not get_ok){
		Log(toolName+" Error! No livetime in CStore!",v_error,verbosity);
		return true;
	}
	livetime = paper_livetime; // XXX XXX REMOVE WHEN IMPLEMENTED
	
	// dt distribution of muons passing dlt<200cm cut
	TH1F dt_mu_lowe_hist("dt_mu_lowe_hist","Spallation Muon to Low-E Time Differences",5000,0,30);
	TH1F dt_mu_lowe_hist_short("dt_mu_lowe_hist_short","Spallation Muon to Low-E Time Differences",500,0,0.25);
	for(auto&& aval : (*dt_mu_lowe_vals)){
		if(aval>0) continue;  // only fit preceding muons (true spallation)
		dt_mu_lowe_hist.Fill(fabs(aval));
		dt_mu_lowe_hist_short.Fill(fabs(aval));
	}
	dt_mu_lowe_hist.Write();
	dt_mu_lowe_hist_short.Write();
	
	// make a histo with equally spaced bins in log scale
	// we need to fit the same histogram (with the same binning) for the intermediate fits,
	// otherwise the bin widths change, the contents change, and the fit parameters aren't the same.
	// 5000 bins (30/0.006) gives chi2/NDOF that matches the paper.
	int nbins=5000;
	std::vector<double> binedges = MakeLogBins(0.001, 30, nbins+1);
	TH1F dt_mu_lowe_hist_log("dt_mu_lowe_hist_log","Data;dt(s);Events/0.006 s",
							 nbins, binedges.data());
	for(auto&& aval : (*dt_mu_lowe_vals)){ if(aval>0) continue; dt_mu_lowe_hist_log.Fill(fabs(aval)); }
	// scale each bin by its width to correct back to number of events per equal time interval
	for(int bini=1; bini<dt_mu_lowe_hist_log.GetNbinsX()+1; ++bini){
		// numbers here are from binning of dt_mu_lowe_hist
		dt_mu_lowe_hist_log.SetBinContent(bini,
			dt_mu_lowe_hist_log.GetBinContent(bini)/dt_mu_lowe_hist_log.GetBinWidth(bini));
	}
	// for some ungodly reason they plot things in events / 0.006s
	// so scale down by events / second to events / 0.006s by * 0.006s
	dt_mu_lowe_hist_log.Scale(0.006);
//	std::cout<<"this plot will be used for all the time distribution fits:"<<std::endl;
//	dt_mu_lowe_hist_log.Draw();
//	gPad->WaitPrimitive();
	
	// fit this with all the rates....
	// we do the fitting in 5 stages, initially fitting sub-ranges of the distribution
	for(int i=0; i<5; ++i) FitDtDistribution(dt_mu_lowe_hist, dt_mu_lowe_hist_log, i);
	
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
		if(isotope.find("_")==std::string::npos){
			double energy_cut_eff = reco_effs_8mev.at(anisotope.first)/100.;
			// FIXME for now just assume the same 1st reduction and dlt cut efficiency
			double efficiency = energy_cut_eff * (paper_first_reduction_eff/100.) * (paper_dlt_cut_eff/100.);
			rates[anisotope.first] = anisotope.second/(efficiency * fiducial_vol * livetime);
			std::cout<<"Num of "<<anisotope.first<<" events is "<<anisotope.second
					 <<", energy cut efficiency is "<<energy_cut_eff<<", total efficiency is "<<efficiency
					 <<" giving a total rate of "<<rates[anisotope.first]<<"/kton/day"<<std::endl;
		} else {
			std::string first_isotope = isotope.substr(0,isotope.find_first_of("_"));
			std::string second_isotope = isotope.substr(isotope.find_first_of("_")+1,std::string::npos);
			double e_cut_eff_1 = reco_effs_8mev.at(first_isotope)/100.;
			double e_cut_eff_2 = reco_effs_8mev.at(second_isotope)/100.;
			// FIXME for now just assume the same 1st reduction and dlt cut efficiency
			double first_efficiency = e_cut_eff_1 * (paper_first_reduction_eff/100.) * (paper_dlt_cut_eff/100.);
			double second_efficiency = e_cut_eff_2 * (paper_first_reduction_eff/100.) * (paper_dlt_cut_eff/100.);
			rates[first_isotope] = 0.5*anisotope.second/(first_efficiency * fiducial_vol * livetime);
			rates[second_isotope] = 0.5*anisotope.second/(second_efficiency * fiducial_vol * livetime);
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

std::vector<double> FitSpallationDt::MakeLogBins(double xmin, double xmax, int nbins){
	std::vector<double> binedges(nbins+1);
	double xxmin=log10(xmin);
	double xxmax = log10(xmax);
	for(int i=0; i<nbins; ++i){
		binedges[i] = pow(10,xxmin + (double(i)/(nbins-1.))*(xxmax-xxmin));
	}
	binedges[nbins+1] = xmax; // required
	return binedges;
}

bool FitSpallationDt::FitDtDistribution(TH1F& dt_mu_lowe_hist, TH1F& dt_mu_lowe_hist_log, int rangenum){
	// do sub-range fitting based on section B of the 2015 paper
	// formula 1 from the paper, with a sub-range as described in section B2
	// we have 4 sub-ranges to fit
	switch (rangenum){
	case 0:{
		// fit range 50us -> 0.1s with 12B + 12N only
		// make the function which we'll fit
		std::cout<<"calling BuildFunction for case "<<rangenum<<", 12B+12N"<<std::endl;
		TF1 func_sum = BuildFunction({"12B","12N"},50e-6,0.1);
		// do the fit
		std::cout<<"fitting"<<std::endl;
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
		dt_mu_lowe_hist_log.Draw();
		//gPad->WaitPrimitive();
		dt_mu_lowe_hist_log.GetListOfFunctions()->Clear();
		break;
		}
	case 1:{
		// fit range 6-30s with 16N + 11B only
		std::cout<<"calling BuildFunction for case "<<rangenum<<", 16N+11Be"<<std::endl;
		TF1 func_sum = BuildFunction({"16N","11Be"},6,30);
		std::cout<<"fitting"<<std::endl;
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
		dt_mu_lowe_hist_log.Draw();
		//gPad->WaitPrimitive();
		dt_mu_lowe_hist_log.GetListOfFunctions()->Clear();
		break;
		}
	case 2:{
		// fit the range 0.1-0.8s with the components previously fit now fixed,
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
		// allowing additional components 15C + 16N
		std::cout<<"calling BuildFunction for case "<<rangenum<<", 15C+16N+8Li_8B"<<std::endl;
		TF1 func_sum = BuildFunction({"15C","16N","8Li_8B"},0.8,6);
		std::cout<<"retrieving past results"<<std::endl;
		PullFitAmp(func_sum,"8Li_8B");
		// fit the new components
		std::cout<<"fitting"<<std::endl;
		dt_mu_lowe_hist_log.Fit(&func_sum,"R","",0.8,6);
		// record the results for the next step
		std::cout<<"recording fit results"<<std::endl;
		PushFitAmp(func_sum,"15C");
		PushFitAmp(func_sum,"16N");
		fit_amps["const_3"] = func_sum.GetParameter("const");
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
		
		std::cout<<"fitting"<<std::endl;
		TFitResultPtr fitresult =  dt_mu_lowe_hist_log.Fit(&func_sum,"RS","",0,30);
		
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
		dt_mu_lowe_hist_log.Draw();
		gPad->SetLogx();
		gPad->SetLogy();
		
		// since we're having issues, let's see how well the old results fit our data
		TF1 func_paper = BuildFunction2({"12B","12N","16N","11Be","9Li","8He_9C","8Li","8B","15C"},0,30);
		// scale the paper down to be around the same num events as us, to ease comparison
		std::cout<<"retrieving paper results"<<std::endl;
		bool correct_energy_threshold = false;
		PullPaperAmp(func_paper,"12B",correct_energy_threshold,paper_scaling);
		PullPaperAmp(func_paper,"12N",correct_energy_threshold,paper_scaling);
		PullPaperAmp(func_paper,"16N",correct_energy_threshold,paper_scaling);
		PullPaperAmp(func_paper,"11Be",correct_energy_threshold,paper_scaling);
		PullPaperAmp(func_paper,"9Li",correct_energy_threshold,paper_scaling);
		PullPaperAmp(func_paper,"8He_9C",correct_energy_threshold,paper_scaling);
		//PullPaperAmp(func_paper,"8Li_8B",correct_energy_threshold,paper_scaling);
		PullPaperAmp(func_paper,"8Li",correct_energy_threshold,paper_scaling);
		PullPaperAmp(func_paper,"8B",correct_energy_threshold,paper_scaling);
		PullPaperAmp(func_paper,"15C",correct_energy_threshold,paper_scaling);
		PullPaperAmp(func_paper,"const",correct_energy_threshold,paper_scaling);
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

void FitSpallationDt::BuildPaperPlot(){
	// add the individual isotopic contributions to reproduce the old plot
	std::vector<TF1> indiv_funcs;
	std::vector<std::string> isotope_names;
	indiv_funcs.reserve(papervals.size());
	for(auto&& theisotope : papervals){
		if(theisotope.first.substr(0,5)=="const") continue; // not a real isotope
		if(theisotope.second==0) continue; // do not add to plot isotopes with no abundance
		
		std::cout<<"****** buildpaperplot rate of isotope "
				 <<theisotope.first<<" is "<<theisotope.second<<std::endl;
		std::string anisotope = theisotope.first;
		TF1 next_func = BuildFunction({anisotope},0,30);
		
		PullPaperAmp(next_func,anisotope,false,paper_scaling);
		next_func.SetLineColor(colourwheel.GetNextColour());
		
		indiv_funcs.push_back(next_func);
		indiv_funcs.back().SetName(anisotope.c_str());
		indiv_funcs.back().SetTitle(anisotope.c_str());
		isotope_names.push_back(anisotope);
	}
	// and last but not least the constant background
	TF1 constfunc("constfunc","[0]",0,30);
	constfunc.SetParameter(0,papervals.at("const"));
	indiv_funcs.push_back(constfunc);
	indiv_funcs.back().SetName("const");
	indiv_funcs.back().SetTitle("const");
	
	// now the sum of everything
	TF1 func_paper = BuildFunction(isotope_names,0.001,30);
	func_paper.SetName("total");
	func_paper.SetTitle("total");
	// set the amplitudes
	for(auto&& theisotope : isotope_names){
		PullPaperAmp(func_paper,theisotope,false,paper_scaling);
	}
	func_paper.SetParameter("const",papervals.at("const"));
	
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
TF1 FitSpallationDt::BuildFunction2(std::vector<std::string> isotopes, double func_min, double func_max){
	std::string total_func="";
	std::string func_name="";
	std::map<std::string,int> parameter_posns;
	int next_par_index=0;
	for(std::string& anisotope : isotopes){
		func_name += anisotope+"_";
		//std::cout<<"adding isotope "<<anisotope<<std::endl;
		// first, this 'isotope' may be a degenerate pair, so try to split it apart
		if(anisotope.find("_")==std::string::npos){
			// not a pair
			//std::cout<<"not a pair"<<std::endl;
			int first_index=next_par_index;
			// "([0]/[1])*exp(-x/[1])"
			
			// Or that's the function you'd expect. Fig 3 plots from 0--30s, with x-axis in seconds;
			// which means that F(t) = dN/dt is also in seconds ... but Fig 3's y-axis is in events / 0.006s!!
			// The bin width is variable (it's a log-log plot), so we already need to scale our bin counts
			// by the bin width to get consistent units, so there's no reason not to use events/second.
			// Still, to make a comparable plot, we could scale our histogram bin counts up using TH1::Scale,
			// but then our fit values will be off unless our fit function accounts for it.
			// (this also ensures the paper plot overlay comparison has the correct scaling).
			
			std::string scalestring = "0.006*";
			std::string this_func =  scalestring+"(abs(["+toString(next_par_index)
									+"])/["+toString(next_par_index+1)
									+"])*exp(-x/["+toString(next_par_index+1)+"])";
			std::cout<<"function is "<<this_func<<std::endl;
			next_par_index +=2;
			// add this function to the total function string
			total_func = (total_func=="") ? this_func : total_func+" + "+this_func;
			// add the parameter names to our map
			parameter_posns.emplace("amp_"+anisotope,first_index);
			parameter_posns.emplace("lifetime_"+anisotope,first_index+1);
		} else {
			// it's a pair. For now, only support two isotopes at a time.
			//std::cout<<"a pair, splitting into ";
			std::string first_isotope = anisotope.substr(0,anisotope.find_first_of("_"));
			std::string second_isotope = anisotope.substr(anisotope.find_first_of("_")+1,std::string::npos);
			//std::cout<<first_isotope<<" and "<<second_isotope<<std::endl;
			// the fit function isn't just the sum of two single isotope functions
			// as they share an amplitude and constant
			//"[0]*0.5*(exp(-x/[1])/[1]+exp(-x/[2])/[2])
			// as above, add in the additional scaling factor to get Y units of events/0.006s
			std::string scalestring = "0.006*";
			int first_index=next_par_index;
			std::string this_func =  scalestring+"abs(["+toString(next_par_index)+"])*0.5*"
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
		std::cout<<"setting parname "<<next_par.second<<" to "<<next_par.first<<std::endl;
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
TF1 FitSpallationDt::BuildFunction(std::vector<std::string> isotopes, double func_min, double func_max){
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
			
			// get paper amplitude, corrected for energy efficiency and scaled by config file scaling
			double paperval = GetPaperAmp(anisotope,true,paper_scaling);
			std::string paperstring = std::to_string(paperval/2.); // fix half the paper val, fit the rest
			// "((x.xx + abs([0]))/[1])*exp(-x/[1])"
			std::string scalestring = "0.006*";
			std::string this_func =  scalestring+"(("+paperstring+"+abs(["+toString(next_par_index)
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
			
			// get paper amplitude, corrected for energy efficiency and scaled by config file scaling
			std::cout<<"getting paper amp"<<std::endl;
			double paperval = GetPaperAmp(anisotope,true,paper_scaling);
			std::cout<<"paperval= "<<paperval<<std::endl;
			
			// fix the abundance to at least half the paper val, fit the rest
			std::string paperstring = std::to_string(paperval/2.);
			std::string scalestring = "0.006*";
			std::string this_func =  scalestring+"("+paperstring+"+abs(["+toString(next_par_index)+"]))*0.5*"
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
void FitSpallationDt::FixLifetime(TF1& func, std::string isotope){
	int par_number = func.GetParNumber(("lifetime_"+isotope).c_str());
	func.FixParameter(par_number,lifetimes.at(isotope));
}

// copy fit result amplitude value back into a component TF1
void FitSpallationDt::PushFitAmp(TF1& func, std::string isotope){
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
void FitSpallationDt::PushFitAmp(double amp, std::string isotope){
	fit_amps[isotope] = amp;
}

// copy amplitude value from a component TF1 into a sum function for fitting
void FitSpallationDt::PullFitAmp(TF1& func, std::string isotope, bool fix){
	std::string parname = (isotope.substr(0,5)=="const") ? "const" : "amp_"+isotope;
	int par_number = func.GetParNumber(parname.c_str());
	if(fix){
		func.FixParameter(par_number,fit_amps.at(isotope));
	} else {
		func.SetParameter(par_number,fit_amps.at(isotope));
	}
}

double FitSpallationDt::GetPaperAmp(std::string isotope, bool threshold_scaling, double fixed_scaling){
	std::cout<<"getting paper amp for isotope "<<isotope<<std::endl;
	double paperval;
	if(isotope!="const"){
		// the values in papervals are rates in events / kton / day, integrated over all beta energies.
		// To obtain Ni, the initial number of observable events over our live time,
		// we need to multiply by the livetime, the fiducial volume, and the efficiency of observation
		// (i.e. total efficiency all cuts)
		// Everything except the low energy cut is isotope independent.
		// the efficiency of the energy cut we need to look up for each isotope
		// The threshold_scaling bool defines whether we use the 6MeV efficiency from the paper,
		// or an 8MeV efficiency from FLUKA + skdetsim. For the same livetime this would give the number
		// of events we would expect to see with a higher E threshold so that we can compare.
		double energy_cut_eff = 1.;
		if(isotope.find("_")==std::string::npos){
			// not a pairing - can look up efficiency directly
			if(threshold_scaling){
				energy_cut_eff = reco_effs_8mev.at(isotope);
			} else {
				energy_cut_eff = papereffs.at(isotope);
				std::cout<<" 6MeV cut efficiency is "<<energy_cut_eff<<std::endl;
			}
			paperval = papervals.at(isotope);
		} else {
			std::string first_isotope = isotope.substr(0,isotope.find_first_of("_"));
			std::string second_isotope = isotope.substr(isotope.find_first_of("_")+1,std::string::npos);
			// how do we handle efficiency for pairs?
			// best we can do is assume half for each and average the efficiency i think....
			if(threshold_scaling){
				energy_cut_eff = 0.5*(reco_effs_8mev.at(first_isotope) + reco_effs_8mev.at(second_isotope));
			} else {
				energy_cut_eff = 0.5*(papereffs.at(first_isotope) + papereffs.at(second_isotope));
			}
			// for a pair of isotopes we may have one or both
			if(papervals.count(isotope)){
				// if we have a combined amplitude, use that
				paperval = papervals.at(isotope);
			} else {
				// otherwise take the sum
				double paperval1 = papervals.at(first_isotope);
				double paperval2 = papervals.at(second_isotope);
				paperval = (paperval1+paperval2);
			}
		}
		// ok, convert Ri to Ni
		paperval *= fiducial_vol * paper_livetime * (paper_first_reduction_eff/100.)
					 * (paper_dlt_cut_eff/100.) * (energy_cut_eff/100.);
		std::cout<<"paper Ni is "<<paperval<<std::endl;
	} else {
		// the constant term isn't a rate, it's read straight off the plot so doesn't need conversion.
		// But, if we're comparing across energy thresholds, it should also be scaled
		// down by the relative rate of random backgrounds above 8MeV vs above 6MeV
		// TODO obtain that number...
		paperval = papervals.at(isotope); // for now just use the value read off the plot, no scaling
	}
	
	paperval *= fixed_scaling; // superfluous, just in case
	return paperval;
}

void FitSpallationDt::PullPaperAmp(TF1& func, std::string isotope, bool threshold_scaling, double fixed_scaling){
	double paperval = GetPaperAmp(isotope, threshold_scaling, fixed_scaling);
	std::string parname = (isotope=="const") ? "const" : "amp_"+isotope;
	func.SetParameter(parname.c_str(), paperval);
}

// =========================================================================
// =========================================================================

bool FitSpallationDt::GetEnergyCutEfficiencies(){
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
		
		reco_effs_6mev.emplace(isotope,reco_eff6*100.);
		reco_effs_8mev.emplace(isotope,reco_eff8*100.);
		reco_effs_scaling.emplace(isotope,(reco_eff8/reco_eff6));
	}
	return true;
}
