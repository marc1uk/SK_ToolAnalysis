/* vim:set noexpandtab tabstop=4 wrap */
#include "PlotNeutronCaptures.h"

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"

#include <memory>  // unique_ptr
#include <map>

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPie.h"

// TODO move to DataModel RootAlgorithms or something
std::unique_ptr<TPie> GeneratePieFromHisto(TH1F* histo, int verbose=0);
std::unique_ptr<TPie> GeneratePieFromHisto(std::string histoname, int verbose=0);

PlotNeutronCaptures::PlotNeutronCaptures():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

bool PlotNeutronCaptures::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="") m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",verbosity);            // how verbose to be
	m_variables.Get("inputFile",inputFile);            // a single specific input file
	m_variables.Get("outputFile",outputFile);          // output file to write
	m_variables.Get("drawPlots",DrawPlots);            // show root plots while working? otherwise just save
	m_variables.Get("maxEvents",maxEvents);            // user limit to number of events to process
	
	// open the input TFile and TTree
	// ------------------------------
	get_ok = myTreeReader.Load(inputFile, "eventtree");
	intree = myTreeReader.GetTree();
	
	// open the output TFile and TTree
	// -------------------------------
	outfile = new TFile(outputFile.c_str(),"RECREATE");
	friendTree = new TTree("ntree","Process Variables");
	friendTree->Branch("neutrino_momentum",&neutrino_momentump,32000,0);
	friendTree->Branch("muon_momentum",&muon_momentump,32000,0);
	// travel distance components relative to neutrino dir
	friendTree->Branch("neutron_longitudinal_travel",&neutron_longitudinal_travel,32000,0);
	friendTree->Branch("neutron_perpendicular_travel",&neutron_perpendicular_travel,32000,0);
	friendTree->Branch("total_gamma_energy",&total_gamma_energy,32000,0);
	
	if(DrawPlots){
		// Only one TApplication may exist. Get it, or make it if there isn't one
		int myargc=0;
		intptr_t tapp_ptr=0;
		get_ok = m_data->CStore.Get("RootTApplication",tapp_ptr);
		if(not get_ok){
			if(verbosity>2) std::cout<<toolName<<": making global TApplication"<<std::endl;
			rootTApp = new TApplication("rootTApp",&myargc,0);
			tapp_ptr = reinterpret_cast<intptr_t>(rootTApp);
			m_data->CStore.Set("RootTApplication",tapp_ptr);
		} else {
			if(verbosity>2) std::cout<<toolName<<": Retrieving global TApplication"<<std::endl;
			rootTApp = reinterpret_cast<TApplication*>(tapp_ptr);
		}
		int tapplicationusers;
		get_ok = m_data->CStore.Get("RootTApplicationUsers",tapplicationusers);
		if(not get_ok) tapplicationusers=1;
		else tapplicationusers++;
		m_data->CStore.Set("RootTApplicationUsers",tapplicationusers);
	}
	
	return true;
}


bool PlotNeutronCaptures::Execute(){
	
	Log(toolName+" processing entry "+toString(entrynum),v_debug,verbosity);
	
	// ttree entry is already loaded so just retrieve the desired branches
	Log(toolName+" getting data from input branches",v_debug,verbosity);
	get_ok = GetBranches();
	
	// process the data
	Log(toolName+" calculating processed data for output tree",v_debug,verbosity);
	get_ok = FillFriend();
	
	// move to next entry
	Log(toolName+" checking if stopping the toolchain",v_debug,verbosity);
	entrynum++;
	// check if we've hit the user-requested entry limit
	if((maxEvents>0)&&(entrynum==maxEvents)){
		Log(toolName+" hit max events, setting StopLoop",v_message,verbosity);
		m_data->vars.Set("StopLoop",1);
		return 1;
	}
	
	// pre-load the next ttree entry
	Log(toolName+" pre-loading entry "+toString(entrynum),v_debug,verbosity);
	get_ok = ReadEntry(entrynum);
	if(get_ok==0){
		return 1; // end of file
	} else if (get_ok<0){
		return 0; // read error
	}
	
	return get_ok;
}

int PlotNeutronCaptures::FillFriend(){
	// don't carry anything over
	Log(toolName+" clearing output ttree variables",v_debug,verbosity);
	ClearOutputTreeBranches();
	
	// loop over primaries and extract the neutrino and primary muon,
	// since we want their momenta for later derived values
	Log(toolName+" getting primary nu/mu momenta",v_debug,verbosity);
	int neutrino_pdg = 12;    // recall that skdetsim has just one neutrino type, which gets saved as ν-e
	int muon_pdg = 13;
	for(int primary_i=0; primary_i<primary_pdg->size(); ++primary_i){
		int next_primary_pdg = primary_pdg->at(primary_i);
		if(next_primary_pdg==neutrino_pdg){
			neutrino_momentum = primary_start_mom->at(primary_i);
		} else if(next_primary_pdg==muon_pdg){
			muon_momentum = primary_start_mom->at(primary_i);
		}
	}
	
	// loop over neutrons in this entry and build the auxilliary info for the friend tree
	Log(toolName+" calculating neutron travel components",v_debug,verbosity);
	for(int neutron_i=0; neutron_i<neutron_start_pos->size(); ++neutron_i){
		// longitudinal distance = (neutron_travel_vector).(neutrino_direction_vector)
		TVector3 neutron_travel_vector = 
			neutron_end_pos->at(neutron_i).Vect() - neutron_start_pos->at(neutron_i).Vect();
		double next_neutron_longitudinal_travel = neutron_travel_vector.Dot(neutrino_momentum.Unit());
		double next_neutron_perpendicular_travel = neutron_travel_vector.Mag()-next_neutron_longitudinal_travel;
		neutron_longitudinal_travel.push_back(next_neutron_longitudinal_travel);
		neutron_perpendicular_travel.push_back(next_neutron_perpendicular_travel);
		
		double total_gamma_E=0;
		for(auto&& agamma : gamma_energy->at(neutron_i)){
			total_gamma_E+=agamma;
		}
		total_gamma_energy.push_back(total_gamma_E);
		
		// keep a map with capture nuclides to num capture events
		int capture_nuclide_pdg = nuclide_daughter_pdg->at(neutron_i);
		std::string capture_nuclide_name = PdgToString(nuclide_daughter_pdg->at(neutron_i));
		if(capture_nuclide_vs_count.count(PdgToString(nuclide_daughter_pdg->at(neutron_i)))){
			capture_nuclide_vs_count.at(PdgToString(nuclide_daughter_pdg->at(neutron_i)))++;
		} else {
			capture_nuclide_vs_count.emplace(PdgToString(nuclide_daughter_pdg->at(neutron_i)),1);
		}
	}
	
	// XXX any further event-wise info we want to add to the friend tree?
	Log(toolName+" filling friendTree",v_debug,verbosity);
	friendTree->Fill();
	
//	// angle between vector 'd' (from nu intx vertex and neutron capture)
//	// and "inferred neutron momentum (direction)", 'p', calculated somehow from CCQE assumption...?
//	// expect neutron and muon to have sum of transverse momentum = 0,
//	// -> neutron should be emitted in plane of muon/neutrino tracks,
//	// if we know neutrino direction and muon momentum, we can know muon pt,
//	// which should be balanced by neutron pt (CCQE assumption)
//	// .... but then what? How do we know pl to know expected neutron angle?
//	// is it just angle from neutrino direction??
//	theta_inferred = acos(d.p)/|d|
	
	// TODO:
	// check neutron travel distance as a function of concentration?
	// check neutron capture fraction on nuclei as a function of concentration?
	// compare travel distance to expected distribution for given energy distribution+concentration?
	
	// TODO stretch goals:
	// calculate Evis
	// calculate detection efficiency (efficiency of successful identification? reconstruction? ) (expect 90% on Gd, 25% on H)
	// capture vertex resolution
	// travel distance resolution
	// neutron multiplicity? vs muon pt? - need to revert TruthNeutronCaptures tool to version that keeps nuclides
	
	return 1;
}


bool PlotNeutronCaptures::Finalise(){
	
	// write out the friend tree
	Log(toolName+" writing output TTree",v_debug,verbosity);
	outfile->cd();
	friendTree->Write();
	
	// make and write out histograms
	Log(toolName+" making histograms",v_debug,verbosity);
	MakeHistos();
	
	if(friendTree) friendTree->ResetBranchAddresses();
	if(outfile){ outfile->Close(); delete outfile; outfile=nullptr; }
	
	// unregister ourselves as a user of the TApplication, and delete it if we're the last one
	if(DrawPlots){
		int tapplicationusers=0;
		get_ok = m_data->CStore.Get("RootTApplicationUsers",tapplicationusers);
		if(not get_ok || tapplicationusers==1){
			if(rootTApp){
				std::cout<<toolName<<": Deleting global TApplication"<<std::endl;
				delete rootTApp;
				rootTApp=nullptr;
			}
		} else if(tapplicationusers>1){
			m_data->CStore.Set("RootTApplicationUsers",tapplicationusers-1);
		}
	}
	
	return true;
}

int PlotNeutronCaptures::ReadEntry(long entry_number){
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

int PlotNeutronCaptures::MakeHistos(){
	// the lazy way (proper would be to call TH1::Fill during FillFriend)
	// ============
	outfile->cd();
	
	// FIXME Specify colours for plots added to THStacks
	// figure out how to save TLegend with THStack
	
	// ======================
	// cumulative plots
	// ======================
	Log(toolName+" making aggregate plots",v_debug,verbosity);
	// neutron energy
	TH1D hNeutronE("hNeutronE","Neutron Energy;Neutron Energy [MeV];Num Events",100,0,100);
	intree->Draw("neutron_start_energy>>hNeutronE");
	hNeutronE.Write();
	
	// neutron travel distance
	TH1D hNeutronTravelDist("hNeutronTravelDist","Neutron Travel Distance;Distance [cm];Num Events",100,0,2200);
	intree->Draw("n_travel_dist>>hNeutronTravelDist");
	hNeutronTravelDist.Write();
	
	// gamma mulitiplicity
	TH1D hGammaNum("hGammaNum","Gamma Multiplicity (All Nuclides);Num Gammas Emitted;Num Events",100,0,20);
	intree->Draw("neutron_n_daughters>>hGammaNum");
	hGammaNum.Write();
	
	// gamma energy
	TH1D hGammaE("hGammaE", "Gamma Energy (All Nuclides);Gamma Energy [MeV];Num Events",100,0,10);
	intree->Draw("gamma_energy>>hGammaE");
	hGammaE.Write();
	
	// total gamma energy from the neutron capture
	TH1D hSumGammaE("hSumGammaE","Total Emitted Gamma Energy (All Nuclides);Sum of Gamma Energy [MeV];Num Events",100,0,10);
	friendTree->Draw("total_gamma_energy>>hSumGammaE");
	hSumGammaE.Write();
	
	// gamma emission time (parent lifetime)
	TH1D hGammaT("hGammaT", "Gamma Emission Time (All Nuclides);Gamma Emission Time [ns];Num Events",100,0,100);
	intree->Draw("gamma_time>>hGammaT");
	hGammaT.Write();
	
	// pie chart of capture nuclei
	Log(toolName+" making pie chart",v_debug,verbosity);
	TPie pCaptureNuclidePdg = TPie("pCaptureNuclidePdg", "Captures by Nuclide",capture_nuclide_vs_count.size());
	int nuclide_i=0;
	for(auto&& anuclide : capture_nuclide_vs_count){
		pCaptureNuclidePdg.SetEntryLabel(nuclide_i,anuclide.first.c_str());
		pCaptureNuclidePdg.SetEntryVal(nuclide_i,anuclide.second);
		++nuclide_i;
	}
	// making it look nice
//	pCaptureNuclidePdg.SetAngularOffset(333);
	pCaptureNuclidePdg.SetLabelFormat("#splitline{%txt}{#splitline{%val}{(%perc)}}");
	pCaptureNuclidePdg.SetValueFormat("%4.0f");
	pCaptureNuclidePdg.SetPercentFormat("%3.0f");
	pCaptureNuclidePdg.SetCircle(0.5, 0.4702026, 0.3302274);
	pCaptureNuclidePdg.SetTextSize(0.03455766);
	// saving to file
	pCaptureNuclidePdg.Draw();
	pCaptureNuclidePdg.Write();
	
	// ==============================
	// broken down by capture nucleus
	// ==============================
	Log(toolName+" making total gamma energy stack",v_debug,verbosity);
	// Total Gamma Energy
	// ------------------
	// capture on H
	TH1D hTotGammaE_H("hTotGammaE_H", "Total Gamma Energy (Capture on H);Gamma Energy [MeV];Num Events",100,0,10);
	// friend the input tree so we can cut on nuclide_daughter_pdg
	friendTree->AddFriend(intree);
	friendTree->Draw("total_gamma_energy>>hTotGammaE_H","nuclide_daughter_pdg==100045");
	
	// capture on Gd-155 -> daughter nuclide Gd-156
	TH1D hTotGammaE_Gd_155("hTotGammaE_Gd_155", "Total Gamma Energy (Capture on Gd-155);Gamma Energy [MeV];Num Events",100,0,10);
	friendTree->Draw("total_gamma_energy>>hTotGammaE_Gd_155","nuclide_daughter_pdg==1000641560");
	
	// capture on Gd-157 -> daughter nuclide Gd-158
	TH1D hTotGammaE_Gd_157("hTotGammaE_Gd_157", "Total Gamma Energy (Capture on Gd-157);Gamma Energy [MeV];Num Events",100,0,10);
	friendTree->Draw("total_gamma_energy>>hTotGammaE_Gd_157","nuclide_daughter_pdg==1000641580");
	
	// Stack of all of them
	THStack hTotGammaE_Stack("hTotGammaE_Stack","Gamma Spectrum by Capture Nucleus;Gamma Energy [MeV];Num Events");
	hTotGammaE_Stack.Add(&hTotGammaE_H);
	hTotGammaE_Stack.Add(&hTotGammaE_Gd_155);
	hTotGammaE_Stack.Add(&hTotGammaE_Gd_157);
	TLegend StackLegend(0.65,0.7,0.88,0.88,NULL);
	StackLegend.SetFillStyle(0);
	StackLegend.SetLineStyle(0);
	StackLegend.AddEntry(&hTotGammaE_H,"Hydrogen","l");
	StackLegend.AddEntry(&hTotGammaE_Gd_155,"Gd-155","l");
	StackLegend.AddEntry(&hTotGammaE_Gd_157,"Gd-157","l");
	hTotGammaE_Stack.Draw();
	StackLegend.Draw();
	hTotGammaE_Stack.Write();
	
	// Gamma Multiplicity
	// -------------------
	Log(toolName+" making gamma multiplicity stack",v_debug,verbosity);
	// capture on H
	TH1D hGammaNum_H("hGammaNum_H", "Gamma Multiplicity (Capture on H);Num Gammas;Num Events",100,0,10);
	intree->Draw("neutron_n_daughters>>hGammaNum_H","nuclide_daughter_pdg==100045");
	
	// capture on Gd-155 -> daughter nuclide Gd-156
	TH1D hGammaNum_Gd_155("hGammaNum_Gd_155", "Gamma Multiplicity (Capture on Gd-155);Num Gammas;Num Events",100,0,10);
	intree->Draw("neutron_n_daughters>>hGammaNum_Gd_155","nuclide_daughter_pdg==1000641560");
	
	// capture on Gd-157 -> daughter nuclide Gd-158
	TH1D hGammaNum_Gd_157("hGammaNum_Gd_157", "Gamma Multiplicity (Capture on Gd-157);Num Gammas;Num Events",100,0,10);
	intree->Draw("neutron_n_daughters>>hGammaNum_Gd_157","nuclide_daughter_pdg==1000641580");
	
	// Stack of all of them
	THStack hGammaNum_Stack("hGammaNum_Stack","Gamma Multiplicity by Capture Nucleus;Num Gammas;Num Events");
	hGammaNum_Stack.Add(&hGammaNum_H);
	hGammaNum_Stack.Add(&hGammaNum_Gd_155);
	hGammaNum_Stack.Add(&hGammaNum_Gd_157);
	hGammaNum_Stack.Draw();
	StackLegend.Draw();
	hGammaNum_Stack.Write();
	
	// Gamma Energy Spectrum
	// ---------------------
	Log(toolName+" making gamma spectrum stack",v_debug,verbosity);
	// capture on H
	TH1D hGammaE_H("hGammaE_H", "Gamma Energy (Capture on H);Gamma Energy [MeV];Num Events",100,0,10);
	intree->Draw("gamma_energy>>hGammaE_H","nuclide_daughter_pdg==100045");
	
	// capture on Gd-155 -> daughter nuclide Gd-156
	TH1D hGammaE_Gd_155("hGammaE_Gd_155", "Gamma Energy (Capture on Gd-155);Gamma Energy [MeV];Num Events",100,0,10);
	intree->Draw("gamma_energy>>hGammaE_Gd_155","nuclide_daughter_pdg==1000641560");
	
	// capture on Gd-157 -> daughter nuclide Gd-158
	TH1D hGammaE_Gd_157("hGammaE_Gd_157", "Gamma Energy (Capture on Gd-157);Gamma Energy [MeV];Num Events",100,0,10);
	intree->Draw("gamma_energy>>hGammaE_Gd_157","nuclide_daughter_pdg==1000641580");
	
	// Stack of all of them
	THStack hGammaE_Stack("hGammaE_Stack","Gamma Spectrum by Capture Nucleus;Gamma Energy [MeV];Num Events");
	hGammaE_Stack.Add(&hGammaE_H);
	hGammaE_Stack.Add(&hGammaE_Gd_155);
	hGammaE_Stack.Add(&hGammaE_Gd_157);
	hGammaE_Stack.Draw();
	StackLegend.Draw();
	hGammaE_Stack.Write();
	
	return 1;
}


int PlotNeutronCaptures::GetBranches(){
	int success = (
//	(myTreeReader.GetBranchValue("filename",filename))                         &&
//	(myTreeReader.GetBranchValue("water_transparency",water_transparency))     &&
//	(myTreeReader.GetBranchValue("entry_number",entry_number))                 &&
//	(myTreeReader.GetBranchValue("subevent_num",subevent_number))              &&
	(myTreeReader.GetBranchValue("primary_pdg",primary_pdg))                   &&
//	(myTreeReader.GetBranchValue("primary_energy",primary_energy))             &&
	(myTreeReader.GetBranchValue("primary_start_mom",primary_start_mom))       &&
//	(myTreeReader.GetBranchValue("primary_start_pos",primary_start_pos))       &&
//	(myTreeReader.GetBranchValue("primary_end_pos",primary_end_pos))           &&
//	(myTreeReader.GetBranchValue("nuclide_parent_pdg",nuclide_parent_pdg))     &&
//	(myTreeReader.GetBranchValue("nuclide_creation_pos",nuclide_creation_pos)) &&
//	(myTreeReader.GetBranchValue("nuclide_decay_pos",nuclide_decay_pos))       &&
	(myTreeReader.GetBranchValue("nuclide_daughter_pdg",nuclide_daughter_pdg)) &&
	(myTreeReader.GetBranchValue("neutron_start_pos",neutron_start_pos))       &&
	(myTreeReader.GetBranchValue("neutron_end_pos",neutron_end_pos))           &&
//	(myTreeReader.GetBranchValue("neutron_start_energy",neutron_start_energy)) &&
//	(myTreeReader.GetBranchValue("neutron_end_energy",neutron_end_energy))     &&
//	(myTreeReader.GetBranchValue("neutron_end_process",neutron_end_process))   &&
//	(myTreeReader.GetBranchValue("neutron_n_daughters",neutron_ndaughters))    &&
	(myTreeReader.GetBranchValue("gamma_energy",gamma_energy))                 //&&
//	(myTreeReader.GetBranchValue("gamma_time",gamma_time))                     &&
	);
	
	return success;
}

void PlotNeutronCaptures::ClearOutputTreeBranches(){
	neutrino_momentum.SetXYZ(0,0,0);
	muon_momentum.SetXYZ(0,0,0);
	neutron_longitudinal_travel.clear();
	neutron_perpendicular_travel.clear();
	total_gamma_energy.clear();
}

// Produce pie chart of nuclei that captured neutrons
// ==================================================
std::unique_ptr<TPie> GeneratePieFromHisto(std::string histoname, int verbose){
	TH1F* histo = (TH1F*)gROOT->FindObject(histoname.c_str());
	if(histo==nullptr){
		std::cerr<<"GeneratePieFromHisto could not find histo "<<histoname<<std::endl;
		return nullptr;
	}
	return GeneratePieFromHisto(histo, verbose);
}

std::unique_ptr<TPie> GeneratePieFromHisto(TH1F* histo, int verbose){
	std::string histoname = std::string(histo->GetName());
	if(verbose) std::cout<<"creating pie chart from histo "<<histoname<<", which has "
						 <<histo->GetNbinsX()<<" bins with contents: "<<std::endl;
	std::vector< std::pair<std::string,float> > histbins;
	for(int bini=0; bini<histo->GetNbinsX(); bini++){
		TString binlabel = histo->GetXaxis()->GetBinLabel(bini+1);
		double binconts = histo->GetBinContent(bini+1);
		if(binconts<0.01) binconts = 0.0f;  // round floats. useful if the histo has been scaled.
		if(verbose && binconts!=0.0f) std::cout<<binlabel.Data()<<" : "<<binconts<<std::endl;
		if(binconts<0) std::cerr<<"error converting "<<histoname<<" to pie chart: bin "<<binlabel.Data()
			<<" has "<<binconts<<" entries!"<<std::endl;
		if(binconts!=0) histbins.emplace_back(binlabel.Data(),binconts);
	}
	
	auto thepie = std::unique_ptr<TPie>(new TPie(TString::Format("%sPie",histoname.c_str()), TString::Format("%s",histoname.c_str()), histbins.size()));
	
	for(int bini=0; bini<histbins.size(); bini++){
		std::pair<std::string,float> abin = histbins.at(bini);
		std::string thebinlabel = abin.first;
		float thebincontents = abin.second;
		if(thebincontents<0) std::cerr<<"error converting "<<histoname<<" to pie chart: bin "<<thebinlabel
			<<" has "<<thebincontents<<" entries!"<<std::endl;
		thepie->SetEntryVal(bini,thebincontents);  // NO +1 - TPie's have no underflow bin!
		thepie->SetEntryLabel(bini,thebinlabel.c_str());
	}
	return thepie;
}
