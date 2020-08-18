/* vim:set noexpandtab tabstop=4 wrap */
#include "PlotNeutrons.h"
#include "Algorithms.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "type_name_as_string.h"

PlotNeutrons::PlotNeutrons():Tool(){
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

bool PlotNeutrons::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	m_data= &data;
	
	// Get the Tool configuration variables
	// ====================================
	m_variables.Get("verbosity",verbosity);           // how verbose to be
	m_variables.Get("drawPlots",DrawPlots);           // show root plots while working? otherwise just save
	m_variables.Get("inputFile",inputFile);           // a single specific input file
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
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
	
	// set the files into the CStore
	if(inputFile==""){
		get_ok = m_data->CStore.Get("InputFileList", input_file_names);
		if(not get_ok){
			Log(toolName+" Error: No inputFile given and no InputFileList in CStore!",v_error,verbosity);
			return false;
		}
	}
	
	return true;
}


bool PlotNeutrons::Execute(){
	
	
	
	
	m_data->vars.Set("StopLoop",1);
	
	return true;
}


bool PlotNeutrons::Finalise(){
	
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

/*
int PlotNeutrons::SetupInputBranches(){
	
	// file level
//	skdetsim_version?
	intree->SetBranchAddress("wlen",&water_transparency);  // [cm]
	
	// event meta info
	intree->SetBranchAddress("nrun",&run_number);
	intree->SetBranchAddress("nsub",&subrun_number);
	intree->SetBranchAddress("nev",&event_number);
//	intree->SetBranchAddress("date",&event_date); // [year,month,day] array
//	intree->SetBranchAddress("time",&time);       // [hour,minute,second,???] array
	
	// event level detector info
	intree->SetBranchAddress("nhit",&N_hit_ID_PMTs);      // "nqisk"
	intree->SetBranchAddress("potot",&total_ID_pes);      // "qismsk"
	intree->SetBranchAddress("pomax",&max_ID_PMT_pes);    // "qimxsk", presumably max # PEs from a single ID PMT
	
	// neutrino interaction info - first primaries array includes neutrino and target (index 0 and 1)
//	intree->SetBranchAddress("mode",&nu_intx_mode);        // use neut_mode_to_string(mode) to decode
//	intree->SetBranchAddress("numnu",&tot_n_primaries);    // both ingoing and outgoing, size of subsequent arrays
	// following are arrays of size numnu
//	intree->SetBranchAddress("ipnu",&primary_pdg);         // see constants::numnu_code_to_string for index mapping
//	intree->SetBranchAddress("pnu",&primary_momentum);     // [GeV/c]
	
	// primary event - second primaries array includes more info
	intree->SetBranchAddress("posv",&primary_event_vertex);      // [cm]
	intree->SetBranchAddress("wallv",&primary_event_dist_from_wall); // [cm]
	intree->SetBranchAddress("npar",&n_outgoing_primaries);      // should be (tot_n_primaries - 2)...
	intree->SetBranchAddress("ipv",&primary_G3_code);            // use constants::g3_to_pdg to map to pdg code
	intree->SetBranchAddress("pmomv",&primary_start_mom);        // [units?] this and ipv are arrays of size npar
	
	// secondaries - first secondaries arrays...
	intree->SetBranchAddress("npar2",&n_secondaries_1);
	// following are arrays of size npar2
	intree->SetBranchAddress("ipv2",&secondary_G3_code_1);                // 
	intree->SetBranchAddress("posv2",&secondary_start_vertex_1);          // array of 3, [cm?] what about time?
	intree->SetBranchAddress("wallv2",&secondary_start_dist_from_wall_1); // [cm?]
	intree->SetBranchAddress("pmomv2",&secondary_start_mom_1);            // [units?]
	intree->SetBranchAddress("iorg",&secondary_origin_1);                 // what is "origin"?
	
	// secondaries - second secondaries array...
	intree->SetBranchAddress("nscndprt",&n_secondaries_2);
	// following are arrays of size nscndprt
	intree->SetBranchAddress("iprtscnd",&secondary_PDG_code_2);      //
	intree->SetBranchAddress("vtxscnd",&secondary_start_vertex_2);   // [units?]
	intree->SetBranchAddress("tscnd",&secondary_start_time_2);       // [ns]? relative to event start?
	intree->SetBranchAddress("pscnd",&secondary_start_mom_2);        // [units?]
	intree->SetBranchAddress("lmecscnd",&secondary_gen_process);     // use constants::G3_process_code_to_string
	intree->SetBranchAddress("nchilds",&secondary_n_daughters);
	intree->SetBranchAddress("ichildidx",&secondary_first_daugher_index); // if >0, 1-based index in this array
	intree->SetBranchAddress("iprntidx",&parent_index);                   // if >0, 1-based index in this array
	
	// further parentage information - still arrays of size nscndprt. Useful?
//	intree->SetBranchAddress("iprntprt",&parent_G3_code);            // or is it a PDG code?
//	intree->SetBranchAddress("pprnt",&parent_mom_at_sec_creation);   // use w/daughter γ to see n energy @ capture
//	intree->SetBranchAddress("vtxprnt",&parent_init_pos);            // [cm?]
//	intree->SetBranchAddress("pprntinit",&parent_init_mom);          // [units?]
//	intree->SetBranchAddress("itrkscnd",&parent_G3_trackid);         // how do we use this?
//	intree->SetBranchAddress("istakscnd",&parent_G3_stack_trackid);  // how do we use this?
//	intree->SetBranchAddress("iprnttrk",&parent_trackid);            // how do we use this?
//	intree->SetBranchAddress("iorgprt",&parent_track_pid_code);      // i'm so confused
	
}

int PlotNeutrons::CreateOutputFile(std::string filename){
	// create the output ROOT file and TTree for writing
	// =================================================
	outfile = new TFile(filename.c_str(), "RECREATE");
	outtree = new TTree("eventtree", "Events with Neutron Captures");
	
	// Set up all the output TTree branches
	// Each TTree entry will correspond to one event, which may have multiple neutrons
	// Each neutron capture may have multiple gammas
	// Each gamma may have multiple PMT hits
	
	// create branches
	// ---------------
	// file level
	outtree->Branch("filename",&filename);
//	outtree->Branch("skdetsim_version",&skdetsim_version);    // where?
//	outtree->Branch("tba_table_version",&tba_table_version);  // where?
	outtree->Branch("water_transparency",&water_transparency);
	
	// event level
	outtree->Branch("entry_number",&entry_number);
	
	// primary particle
//	outtree->Branch("primary_id",&primary_id);
	outtree->Branch("primary_pdg",&primary_pdg);
	outtree->Branch("primary_energy",&primary_energy);
	outtree->Branch("primary_start_pos",&primary_start_pos);
	outtree->Branch("primary_end_pos",&primary_end_pos);
	
	// parent nuclide - one for each neutron
//	outtree->Branch("parent_primary_id",&parent_primary_id);
//	outtree->Branch("nuclide_id",&nuclide_id);
	outtree->Branch("nuclide_pdg",&nuclide_pdg);
	outtree->Branch("nuclide_creation_pos",&nuclide_creation_pos);
	outtree->Branch("nuclide_decay_pos",&nuclide_decay_pos);
	
	// neutron
//	outtree->Branch("parent_nuclide_id",&parent_nuclide_id);
//	outtree->Branch("neutron_id",&neutron_id);
	outtree->Branch("neutron_start_pos",&neutron_start_pos);
	outtree->Branch("neutron_end_pos",&neutron_end_pos);
	outtree->Branch("neutron_start_energy",&neutron_start_energy);
	outtree->Branch("neutron_end_process",&neutron_end_process);
	
	// gamma
//	outtree->Branch("parent_neutron_id",&parent_neutron_id);
	outtree->Branch("gamma_energy",&gamma_energy);
	
	return 1;
}

int PlotNeutrons::WriteTree(){
	if(verbosity>2) std::cout<<"writing TTree"<<std::endl;
	outfile->cd();
	// TObject::Write returns the total number of bytes written to the file.
	// It returns 0 if the object cannot be written.
	int bytes = outtree->Write("",TObject::kOverwrite);
	if(bytes<=0){
		std::cerr<<"ShowNeutrons Error writing TTree!"<<std::endl;
	} else if(verbosity>2){
		std::cout<<"Wrote "<<get_ok<<" bytes"<<std::endl;
	}
	return bytes;
};

void PlotNeutrons::CloseFile(){
	outtree->ResetBranchAddresses();
	outfile->Write("*",TObject::kOverwrite);
	outfile->Close();
	delete outfile;
	outfile = nullptr;
};

void PlotNeutrons::UpdateStorageArrays(){
	// because some branches stores static arrays we need to make sure
	// that our internal arrays have enough space before we call GetEntry
	
	// there are two sets of branches that store information about primary particles
	// the sizes of the two sets of corresponding arrays are...
	TBranch* numnubranch = intree->GetBranch("numnu"); // incident and outgoing particles to primary interaction
	TBranch* nparbranch  = intree->GetBranch("npar");  // just outgoing particles from primary interaction
	// likewise there are two sets of branches that store information about secondary particles
	// the branches storing their array sizes are...
	TBranch* npar2branch = intree->GetBranch("npar2");         // array of secondaries 1
	TBranch* nscndprtbranch  = intree->GetBranch("nscndprt");  // array of secondaries 2
	
	// retrieve their values
	numnubranch->GetEntry(entry_number);       // numnu    is read into tot_n_primaries
	nparbranch->GetEntry(entry_number);        // npar     is read into n_outgoing_primaries
	npar2branch->GetEntry(entry_number);       // npar2    is read into n_secondaries_1
	nscndprtbranch->GetEntry(entry_number);    // nscndprt is read into n_secondaries_2
	
	// resize our vectors
	if(expand_primaries){
		// primary array set 1
		primary_pdg.resize(tot_n_primaries);
		primary_momentum.resize(tot_n_primaries);
		// primary array set 2
		primary_G3_code.resize(n_outgoing_primaries);
		primary_start_mom.resize(n_outgoing_primaries);
	}
	
	if(expand_secondaries){
		// secondary array set 1
		secondary_G3_code_1.resize(n_secondaries_1);
		secondary_start_vertex_1.resize(n_secondaries_1);
		secondary_start_dist_from_wall_1.resize(n_secondaries_1);
		secondary_start_mom_1.resize(n_secondaries_1);
		secondary_origin_1.resize(n_secondaries_1);
		
		// secondary array set 2
		secondary_PDG_code_2.resize(n_secondaries_2);
		secondary_start_vertex_2.resize(n_secondaries_2);
		secondary_start_time_2.resize(n_secondaries_2);
		secondary_start_mom_2.resize(n_secondaries_2);
		secondary_gen_process.resize(n_secondaries_2);
		secondary_n_daughters.resize(n_secondaries_2);
		secondary_first_daugher_index.resize(n_secondaries_2);
		parent_index.resize(n_secondaries_2);
//		parent_G3_code.resize(n_secondaries_2);
//		parent_mom_at_sec_creation.resize(n_secondaries_2);
//		parent_init_pos.resize(n_secondaries_2);
//		parent_init_mom.resize(n_secondaries_2);
//		parent_G3_trackid.resize(n_secondaries_2);
//		parent_G3_stack_trackid.resize(n_secondaries_2);
//		parent_trackid.resize(n_secondaries_2);
//		parent_track_pid_code.resize(n_secondaries_2);
	}
	
	
	// most of the time we can just resize our vectors, which is really just updating
	// the internal 'size' variable. But when the new size exceeds the old capacity
	// the internal array will be reallocated and we need to update the branch pointers
	// so compare our reference data pointer with the new one to see if we need to do this
	if(PRIMARY_PDG_DATAP_1 != primary_pdg.data()){
		PRIMARY_PDG_DATAP_1=primary_pdg.data();
		// 


 const char *TLeafObject::GetTypeName() const
 {
    return fTitle.Data();
 }



	// event meta info
	intree->SetBranchAddress("nrun",&run_number);
	intree->SetBranchAddress("nsub",&subrun_number);
	intree->SetBranchAddress("nev",&event_number);
//	intree->SetBranchAddress("date",&event_date); // [year,month,day] array
//	intree->SetBranchAddress("time",&time);       // [hour,minute,second,???] array
	
	// event level detector info
	intree->SetBranchAddress("nhit",&N_hit_ID_PMTs);      // "nqisk"
	intree->SetBranchAddress("potot",&total_ID_pes);      // "qismsk"
	intree->SetBranchAddress("pomax",&max_ID_PMT_pes);    // "qimxsk", presumably max # PEs from a single ID PMT
	
	// neutrino interaction info - first primaries array includes neutrino and target (index 0 and 1)
//	intree->SetBranchAddress("mode",&nu_intx_mode);        // use neut_mode_to_string(mode) to decode
//	intree->SetBranchAddress("numnu",&tot_n_primaries);    // both ingoing and outgoing, size of subsequent arrays
	// following are arrays of size numnu
//	intree->SetBranchAddress("ipnu",&primary_pdg);         // see constants::numnu_code_to_string for index mapping
//	intree->SetBranchAddress("pnu",&primary_momentum);     // [GeV/c]
	
	// primary event - second primaries array includes more info
	intree->SetBranchAddress("posv",&primary_event_vertex);      // [cm]
	intree->SetBranchAddress("wallv",&primary_event_dist_from_wall); // [cm]
	intree->SetBranchAddress("npar",&n_outgoing_primaries);      // should be (tot_n_primaries - 2)...
	intree->SetBranchAddress("ipv",&primary_G3_code);            // use constants::g3_to_pdg to map to pdg code
	intree->SetBranchAddress("pmomv",&primary_start_mom);        // [units?] this and ipv are arrays of size npar
	
	// secondaries - first secondaries arrays...
	intree->SetBranchAddress("npar2",&n_secondaries_1);
	// following are arrays of size npar2
	intree->SetBranchAddress("ipv2",&secondary_G3_code_1);                // 
	intree->SetBranchAddress("posv2",&secondary_start_vertex_1);          // array of 3, [cm?] what about time?
	intree->SetBranchAddress("wallv2",&secondary_start_dist_from_wall_1); // [cm?]
	intree->SetBranchAddress("pmomv2",&secondary_start_mom_1);            // [units?]
	intree->SetBranchAddress("iorg",&secondary_origin_1);                 // what is "origin"?
	
	// secondaries - second secondaries array...
	intree->SetBranchAddress("nscndprt",&n_secondaries_2);
	// following are arrays of size nscndprt
	intree->SetBranchAddress("iprtscnd",&secondary_PDG_code_2);      //
	intree->SetBranchAddress("vtxscnd",&secondary_start_vertex_2);   // [units?]
	intree->SetBranchAddress("tscnd",&secondary_start_time_2);       // [ns]? relative to event start?
	intree->SetBranchAddress("pscnd",&secondary_start_mom_2);        // [units?]
	intree->SetBranchAddress("lmecscnd",&secondary_gen_process);     // use constants::G3_process_code_to_string
	intree->SetBranchAddress("nchilds",&secondary_n_daughters);
	intree->SetBranchAddress("ichildidx",&secondary_first_daugher_index); // if >0, 1-based index in this array
	intree->SetBranchAddress("iprntidx",&parent_index);                   // if >0, 1-based index in this array
	
	// further parentage information - still arrays of size nscndprt. Useful?
//	intree->SetBranchAddress("iprntprt",&parent_G3_code);            // or is it a PDG code?
//	intree->SetBranchAddress("pprnt",&parent_mom_at_sec_creation);   // use w/daughter γ to see n energy @ capture
//	intree->SetBranchAddress("vtxprnt",&parent_init_pos);            // [cm?]
//	intree->SetBranchAddress("pprntinit",&parent_init_mom);          // [units?]
//	intree->SetBranchAddress("itrkscnd",&parent_G3_trackid);         // how do we use this?
//	intree->SetBranchAddress("istakscnd",&parent_G3_stack_trackid);  // how do we use this?
//	intree->SetBranchAddress("iprnttrk",&parent_trackid);            // how do we use this?
//	intree->SetBranchAddress("iorgprt",&parent_track_pid_code);      // i'm so confused









		
		
	
	if(reset_primary_branches_1){
		MAX_N_PRIMARIES = primary_pdg.capacity();
		// set branch addresses

*/
