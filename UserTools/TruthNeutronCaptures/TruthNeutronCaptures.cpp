/* vim:set noexpandtab tabstop=4 wrap */
#include "TruthNeutronCaptures.h"
#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"

#include <algorithm> // std::find

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TVector3.h"
#include "TLorentzVector.h"

TruthNeutronCaptures::TruthNeutronCaptures():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

bool TruthNeutronCaptures::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	m_data= &data;
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",verbosity);            // how verbose to be
	m_variables.Get("inputFile",inputFile);            // a single specific input file
	m_variables.Get("outputFile",outputFile);          // output file to write
	m_variables.Get("maxEvents",MAX_EVENTS);           // terminate after processing at most this many events
	m_variables.Get("writeFrequency",WRITE_FREQUENCY); // how many events to TTree::Fill between TTree::Writes
	
	// get the list of input files from the CStore
	// -------------------------------------------
	// filled if using LoadFileList tool
	// TODO yet to implement support for this in MTreeReader
	if(inputFile==""){
		get_ok = m_data->CStore.Get("InputFileList", input_file_names);
		if(not get_ok){
			Log(toolName+" Error: No inputFile given and no InputFileList in CStore!",v_error,verbosity);
			return false;
		}
	}
	
	// open the input TFile and TTree
	// ------------------------------
	get_ok = myTreeReader.Load(inputFile, "h1"); // official ntuple TTree is descriptively known as 'h1'
	if(get_ok) ReadEntryNtuple(0);
	
	// create the output TFile and TTree
	// ---------------------------------
	CreateOutputFile(outputFile);
	
	return true;
}

bool TruthNeutronCaptures::Execute(){
	
	Log(toolName+" processing entry "+toString(entry_number),v_debug,verbosity);
	
	// clear output vectors so we don't carry anything over
	Log(toolName+" clearing output vectors",v_debug,verbosity);
	ClearOutputTreeBranches();
	
	// Copy over directly transferred variables
	Log(toolName+" copying output vectors",v_debug,verbosity);
	CopyVariables();
	
	// Calculate derived variables
	Log(toolName+" calculating output variables",v_debug,verbosity);
	CalculateVariables();
	
	// print the current event
	if(verbosity>1) PrintBranches();
	
	// Fill the output tree
	Log(toolName+"filling output TTree entry",v_debug,verbosity);
	outtree->Fill();
	
	// update the output file so we don't lose everything if we crash
	if((entry_number%WRITE_FREQUENCY)==0) WriteTree();
	
	// stop at user-defined limit to the number of events to process
	++entry_number;
	if((MAX_EVENTS>0)&&(entry_number>=MAX_EVENTS)){
		Log(toolName+" reached MAX_EVENTS, setting StopLoop",v_error,verbosity);
		m_data->vars.Set("StopLoop",1);
	} else {
		// Pre-Load next input entry so we can stop the toolchain
		// if we're about to run off the end of the tree or encounter a read error
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


bool TruthNeutronCaptures::Finalise(){
	
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

void TruthNeutronCaptures::CopyVariables(){
	// copy over variables from the input tree
	// ---------------------------------------
	// those we want to keep in the output tree without modification
	
	out_filename = myTreeReader.GetTree()->GetCurrentFile()->GetName();
//	out_skdetsim_version =    ??
//	out_tba_table_version =   ??
	out_water_transparency = water_transparency;
	
	// event level
//	out_run_number = run_number;
//	out_subrun_number = subrun_number;
//	out_event_number = event_number;
	out_entry_number = entry_number;
	out_subevent_number = subevent_number;
	
}

int TruthNeutronCaptures::CalculateVariables(){
	// calculate remaining variables
	// -----------------------------
	// those we want to save to the output tree but need to derive
	
	int neutron_pdg = 2112; //StringToPdg("Neutron");
	int neutron_g3 = 13; //StringToG3ParticleCode("Neutron");
	int gamma_pdg = 22; //StringToPdg("Gamma");
	double neutron_mass = PdgToMass(neutron_pdg);
	
	// note the primary event vertex. Primary particles don't have individual start vertices,
	// but they *should* all start from the primary event vertex ... unless the geant3 event
	// had primaries at different locations. Which is entirely reasonable. In that case,
	// we may need to use the daughters' "parent vertex at birth" branch to find our primary
	// start vertices.... FIXME for events that should have all particles originating from a
	// single primary event vertex, these are not consistent... why not??
	TLorentzVector primary_vertex_tvector(primary_event_vertex.at(0),
										  primary_event_vertex.at(1),
										  primary_event_vertex.at(2),
										  0.f); // is this correct? does this define T=0?
	
	// loop over event particles and scan for neutrons
	// we want to nest them in their parent nuclides, but we need to find the neutrons first
	// (no way to tell if a nuclide has a daughter neutron)
	// so scan for the neutrons first, then group them by parent after.
	// This means storing the neutron info in temporay container for now
	
	// some helper vectors
	std::map<int,int> primary_n_ind_to_loc;
	std::map<int,int> secondary_n_ind_to_loc;
	std::vector<int> neutron_parent_nuclide_indices; // parent index within the secondaries array
	std::vector<bool> neutron_terminfo_unknown;      // whether we've yet to set the neutron stopping info
	// information that will end up in the output file
	std::vector<TLorentzVector> neutron_start_pos;
	std::vector<TLorentzVector> neutron_end_pos;
	std::vector<double> neutron_start_energy;
	std::vector<double> neutron_end_energy;
	std::vector<int> neutron_end_process;
	std::vector<int> neutron_ndaughters;  // unused. Which branch stores ndaughters for primaries??
	std::vector<std::vector<double> > gamma_energy;
	std::vector<std::vector<double> > gamma_time;
	// used as fall-back for parent nuclide information when the nuclide isn't stored itself
	std::vector<int> parent_nuclide_pdg;
	std::vector<TLorentzVector> parent_nuclide_creation_pos;
	std::vector<int> nuclide_daughter_pdg;
	
	// ==========================
	// SCAN FOR PRIMARY NEUTRONS
	// ==========================
	// loop over the primary particles first, because for IBD events we have a primary neutron
	Log(toolName+" event had "+toString(n_outgoing_primaries)+" primary particles",v_debug,verbosity);
	for(int primary_i=0; primary_i<n_outgoing_primaries; ++primary_i){
		Log(toolName+" primary "+toString(primary_i)+" had G3 code "+toString((int)primary_G3_code.at(primary_i))
				+" ("+G3ParticleCodeToString(primary_G3_code.at(primary_i))+")",v_debug,verbosity);
		if(primary_G3_code.at(primary_i)==neutron_g3){
			// we found a neutron!
			Log(toolName+" NEUTRON!", v_debug,verbosity);
			primary_n_ind_to_loc.emplace(primary_i,neutron_start_energy.size());
			// note neutron info
			neutron_start_pos.push_back(primary_vertex_tvector);
			double startE = pow(primary_start_mom.at(primary_i),2.) / (2.*neutron_mass);
			// FIXME primary momenta appear to be in different units to secondaries...? Need to convert to MeV
			neutron_start_energy.push_back(startE);
			neutron_ndaughters.push_back(0);  // XXX which branch???
			
			// note its parent nuclide information
			neutron_parent_nuclide_indices.push_back(-1); // XXX which branch???
			int parent_nuclide_index = -1;
			if(parent_nuclide_index>=0){
				parent_nuclide_pdg.push_back(secondary_PDG_code_2.at(parent_nuclide_index));
				TLorentzVector parentstartpos(secondary_start_vertex_2.at(parent_nuclide_index).at(0),
						secondary_start_vertex_2.at(parent_nuclide_index).at(1),
						secondary_start_vertex_2.at(parent_nuclide_index).at(2),
						secondary_start_time_2.at(parent_nuclide_index));
				parent_nuclide_creation_pos.push_back(parentstartpos);
			} else {
				// even if we don't have a parent index fill placeholders to keep vectors in sync
				parent_nuclide_pdg.push_back(-1);
				parent_nuclide_creation_pos.push_back(TLorentzVector(0,0,0,0));
			}
			
			// termination information isn't saved, so we'll have to use the creation
			// info of a matched secondary gamma to infer neutron termination information
			// for now, to keep synchronization, fill with placeholders
			neutron_end_pos.push_back(TLorentzVector(0,0,0,0));
			neutron_end_energy.push_back(0);
			neutron_end_process.push_back(0); // do we have any way to determine this?
			neutron_terminfo_unknown.push_back(true);
			
			// a neutron may have many daughters, so we'll have to find those in a subsequent scan
			gamma_energy.push_back(std::vector<double>{});
			gamma_time.push_back(std::vector<double>{});
//			neutron_firstdaughter_indices.push_back(secondary_first_daugher_index.at(secondary_i));
			// likewise we'll note the daughter nuclide from the capture in the subsequent scan
			nuclide_daughter_pdg.push_back(-1);
		}
	}
	Log(toolName+" found "+toString(neutron_start_energy.size())+" primary neutrons", v_debug,verbosity);
	
	// ===========================
	// SCAN FOR SECONDARY NEUTRONS
	// ===========================
	Log(toolName+" event had "+toString(n_secondaries_2)+" secondary particles",v_debug,verbosity);
	for(int secondary_i=0; secondary_i<n_secondaries_2; ++secondary_i){
		Log(toolName+" secondary "+toString(secondary_i)+" had pdg "+toString(secondary_PDG_code_2.at(secondary_i))
				+" ("+PdgToString(secondary_PDG_code_2.at(secondary_i))+")",v_debug,verbosity);
		if(secondary_PDG_code_2.at(secondary_i)==neutron_pdg){
			// we found a neutron!
			Log(toolName+" NEUTRON!", v_debug,verbosity);
			secondary_n_ind_to_loc.emplace(secondary_i,neutron_start_energy.size());
			// note neutron info
			TLorentzVector startpos(secondary_start_vertex_2.at(secondary_i).at(0),
									secondary_start_vertex_2.at(secondary_i).at(1),
									secondary_start_vertex_2.at(secondary_i).at(2),
									secondary_start_time_2.at(secondary_i));
			neutron_start_pos.push_back(startpos);
			TVector3 startmom(secondary_start_mom_2.at(secondary_i).at(0),
							  secondary_start_mom_2.at(secondary_i).at(1),
							  secondary_start_mom_2.at(secondary_i).at(2));
			double startE = startmom.Mag2() / (2.*neutron_mass);
			neutron_start_energy.push_back(startE);
			neutron_ndaughters.push_back(secondary_n_daughters.at(secondary_i));
			
			// note its parent nuclide information
			neutron_parent_nuclide_indices.push_back(parent_index.at(secondary_i)-1);
			int parent_nuclide_index = parent_index.at(secondary_i)-1; // 1-based, 0 if not known
			if(parent_nuclide_index>=0){
				parent_nuclide_pdg.push_back(secondary_PDG_code_2.at(parent_nuclide_index));
				TLorentzVector parentstartpos(secondary_start_vertex_2.at(parent_nuclide_index).at(0),
						secondary_start_vertex_2.at(parent_nuclide_index).at(1),
						secondary_start_vertex_2.at(parent_nuclide_index).at(2),
						secondary_start_time_2.at(parent_nuclide_index));
				parent_nuclide_creation_pos.push_back(parentstartpos);
			} else {
				// even if we don't have a parent index fill placeholders to keep vectors in sync
				parent_nuclide_pdg.push_back(-1);
				parent_nuclide_creation_pos.push_back(TLorentzVector(0,0,0,0));
			}
			
			// termination information isn't saved, so we'll have to use the creation
			// info of a matched secondary gamma to infer neutron termination information
			// for now, to keep synchronization, fill with placeholders
			neutron_end_pos.push_back(TLorentzVector(0,0,0,0));
			neutron_end_energy.push_back(0);
			neutron_end_process.push_back(0); // do we have any way to determine this?
			neutron_terminfo_unknown.push_back(true);
			
			// a neutron may have many daughters, so we'll have to find those in a subsequent scan
			gamma_energy.push_back(std::vector<double>{});
			gamma_time.push_back(std::vector<double>{});
//			neutron_firstdaughter_indices.push_back(secondary_first_daugher_index.at(secondary_i));
		}
	}
	Log(toolName+" found "+toString(neutron_start_energy.size())+" primary+secondary neutrons", v_debug,verbosity);
	
	// ==================================
	// SCAN FOR GAMMAS AND CAPTURE NUCLEI
	// ==================================
	// scan again, this time looking for gammas
	int n_gammas=0;  // for debug
	for(int secondary_i=0; secondary_i<n_secondaries_2; ++secondary_i){
		// ---------------
		// Scan for Gammas
		// ---------------
		if(secondary_PDG_code_2.at(secondary_i)==gamma_pdg){
			Log(toolName+" GAMMA!", v_debug,verbosity);
			n_gammas++;
			// we found a gamma! See if it came from neutron capture (G3 process 18)
			bool from_ncapture = (secondary_gen_process.at(secondary_i)==18);
			Log(toolName+" from ncapture="+toString(from_ncapture),v_debug,verbosity);
			//if(not from_ncapture) continue;  // if not, not interested, skip it
			// for now, as a sanity check, checks if its parent was an identified neutron first
			// parent may either be a primary particle or secondary particle
			// sanity check: it should have one or the other, but not both
			int neutron_parent_loc = -1;
			int primary_parent_index = parent_trackid.at(secondary_i);
			int secondary_parent_index = parent_index.at(secondary_i);
			Log(toolName+" primary parent index "+toString(primary_parent_index)
						+" secondary parent index "+toString(secondary_parent_index),v_debug,verbosity);
			// first check if it has a valid SECONDARY parent
			if(secondary_parent_index>0){
				// its parent was a secondary: check if it's in our list of secondary neutrons
				if(secondary_n_ind_to_loc.count(secondary_parent_index-1)){ // -1 as indices are 1-based
					neutron_parent_loc = secondary_n_ind_to_loc.at(secondary_parent_index-1);
				}
				// else its parent was a secondary, but not one we know
				else if(from_ncapture){
					// if it came from ncapture of a secondary neutron, why don't we know about that neutron?
					Log(toolName+" WARNING, GAMMA FROM NCAPTURE WITH UNKNOWN SECONDARY PARENT INDEX "
							+toString(secondary_parent_index),v_warning,verbosity);
					continue;
				}
			}
			// only fall-back to getting parent PRIMARY if parent secondary index = 0
			// this is because parent primary index is carried over, so daughters of secondaries
			// will have the same primary parent index
			else if(primary_parent_index>0){
				// its parent was a primary: check if it's in our list of primary neutrons
				if(primary_n_ind_to_loc.count(primary_parent_index-1)){  // -1 as indices are 1-based
					neutron_parent_loc = primary_n_ind_to_loc.at(primary_parent_index-1);
				}
				// else its parent was a primary, but not one we know
				else if(from_ncapture){
					// if it came from ncapture of a primary neutron, why don't we know about that neutron?
					Log(toolName+"WARNING, GAMMA FROM NCAPTURE WITH UNKNOWN PRIMARY PARENT INDEX "
							 +toString(primary_parent_index),v_warning,verbosity);
					continue;
				}
			}
			if((neutron_parent_loc<0)&&(from_ncapture)){
				// if we got here, it suggests this gamma is from ncapture
				// but had neither a valid primary or secondary parent index
				Log(toolName+" WARNING, GAMMA FROM NCAPTURE WITH NO PRIMARY OR SECONDARY PARENT INDEX ",
						v_warning,verbosity);
				continue;
			}
			if((neutron_parent_loc>=0)&&(not from_ncapture)){
				// gamma with a parent matched to one of our known neutrons, but not from ncapture?
				// maybe from fast neutron scattering? Leave these for later
				Log(toolName+" WARNING, GAMMA WITH NEUTRON PARENT BUT NOT FROM NCAPTURE",v_warning,verbosity);
				continue;
			}
			// one final check extracts the ones we're after
			if(neutron_parent_loc>=0){
				// neutron from ncapture with known neutron parent! hurray!
				if(neutron_terminfo_unknown.at(neutron_parent_loc)){
					// this is the first daughter gamma we've found for this neutron,
					// so use it's creation info to update the termination info for the neutron
					TLorentzVector gamma_start_pos(secondary_start_vertex_2.at(secondary_i).at(0),
												   secondary_start_vertex_2.at(secondary_i).at(1),
												   secondary_start_vertex_2.at(secondary_i).at(2),
												   secondary_start_time_2.at(secondary_i));
					neutron_end_pos.at(neutron_parent_loc) = gamma_start_pos;
					// get momentum of parent at time of creation of secondary (neutron E at capture)
					TVector3 neutron_end_mom(parent_mom_at_sec_creation.at(secondary_i).at(0),
											 parent_mom_at_sec_creation.at(secondary_i).at(1),
											 parent_mom_at_sec_creation.at(secondary_i).at(2));
					double endE = neutron_end_mom.Mag();
					neutron_end_energy.at(neutron_parent_loc) = endE;
					neutron_end_process.at(neutron_parent_loc) = secondary_gen_process.at(secondary_i);
					neutron_terminfo_unknown.at(neutron_parent_loc) = false;
					
					// double check - for primary neutrons we assume the neutron start pos is
					// the primary event vertex. Compare with the "parent position at creation"
					TVector3 parent_neutron_start_pos(parent_init_pos.at(secondary_i).at(0),
													  parent_init_pos.at(secondary_i).at(1),
													  parent_init_pos.at(secondary_i).at(2));
					if(parent_neutron_start_pos!=neutron_start_pos.at(neutron_parent_loc).Vect()){
						std::cout<<"WARNING, NEUTRON START LOC FROM PARENT POSITION AT BIRTH ("
								 <<parent_neutron_start_pos.X()<<", "<<parent_neutron_start_pos.Y()
								 <<", "<<parent_neutron_start_pos.Z()<<") "
								 <<" DIFFERS FROM NEUTRON START LOC FROM NEUTRON ITSELF ("
								 <<neutron_start_pos.at(neutron_parent_loc).X()<<", "
								 <<neutron_start_pos.at(neutron_parent_loc).Y()<<", "
								 <<neutron_start_pos.at(neutron_parent_loc).Z()<<") "
								 <<std::endl;
					}
					
					// do the same with initial momentum via parent_init_mom
					TVector3 parent_neutron_start_mom(parent_init_mom.at(secondary_i).at(0),
													  parent_init_mom.at(secondary_i).at(1),
													  parent_init_mom.at(secondary_i).at(2));
					double nStartE = parent_neutron_start_mom.Mag2() / (2.*neutron_mass);
					if(nStartE!=neutron_start_energy.at(neutron_parent_loc)){
						std::cout<<"WARNING, NEUTRON START ENERGY FROM PARENT AT BIRTH ("<<nStartE
								 <<") DIFFERS FROM NEUTRON START ENERGY FROM NEUTRON ITSELF ("
								 <<neutron_start_energy.at(neutron_parent_loc)<<") "<<std::endl;
					}
				}
				
				// ok, now add the gamma info
				TVector3 startmom(secondary_start_mom_2.at(secondary_i).at(0),
								  secondary_start_mom_2.at(secondary_i).at(1),
								  secondary_start_mom_2.at(secondary_i).at(2));
				double startE = startmom.Mag();
				double startT = secondary_start_time_2.at(secondary_i);
				gamma_energy.at(neutron_parent_loc).push_back(startE);
				
				gamma_time.at(neutron_parent_loc).push_back(startT);
			}
			// otherwise this is just a plain old gamma. Valid parent, primary or secondary but not both,
			// not from ncapture, and whose parent is not a neutron. Not interested.
		}
		// ------------------------
		// Scan for Daughter Nuclei
		// ------------------------
		else if((secondary_gen_process.at(secondary_i)==18)&&
				(secondary_PDG_code_2.at(secondary_i)!=neutron_pdg)){
			// if it's not a neutron but came from neutron capture it's the daughter nuclide
			// its parent will also be the captured neutron, same as the decay gamma.
			int neutron_parent_loc=-1;
			int primary_parent_index = parent_trackid.at(secondary_i);
			int secondary_parent_index = parent_index.at(secondary_i);
			if(secondary_parent_index>0){
				if(secondary_n_ind_to_loc.count(secondary_parent_index-1)){ // -1 as indices are 1-based
					neutron_parent_loc = secondary_n_ind_to_loc.at(secondary_parent_index-1);
				} else {
					// came from capture of a neutron we don't know?
					Log(toolName+" WARNING, "+PdgToString(secondary_PDG_code_2.at(secondary_i))
						+" FROM NCAPTURE WITH UNKNOWN SECONDARY PARENT (NEUTRON) INDEX "
							+toString(secondary_parent_index),v_warning,verbosity);
					continue;
				}
			} else if(primary_parent_index>0){
				if(primary_n_ind_to_loc.count(primary_parent_index-1)){  // -1 as indices are 1-based
					neutron_parent_loc = primary_n_ind_to_loc.at(primary_parent_index-1);
				} else {
					// came from capture of a neutron we don't know?
					Log(toolName+" WARNING, "+PdgToString(secondary_PDG_code_2.at(secondary_i))
						+" FROM NCAPTURE WITH UNKNOWN PRIMARY PARENT (NEUTRON) INDEX "
							+toString(primary_parent_index),v_warning,verbosity);
					continue;
				}
			}
			if(neutron_parent_loc<0){
				Log(toolName+" WARNING, "+PdgToString(secondary_PDG_code_2.at(secondary_i))
						+" FROM NCAPTURE WITH NO PARENT (NEUTRON) INDEX ",v_warning,verbosity);
				continue;
			} else {
				// nuclide from ncapture with known neutron parent! hurray!
				// Convert from secondary index (i.e. position in array of all secondaries)
				// into neutron index (i.e. position in our array of neutrons)
				if(nuclide_daughter_pdg.at(neutron_parent_loc)>0){
					Log(toolName+" WARNING, FOUND SECOND DAUGHTER NUCLIDE FROM NEUTRON CAPTURE."
						+" FIRST DAUGHTER PDG: "+toString(nuclide_daughter_pdg.at(neutron_parent_loc))
						+" ("+PdgToString(nuclide_daughter_pdg.at(neutron_parent_loc))+") "
						+" SECOND DAUGHTER PDG: "+toString(secondary_PDG_code_2.at(secondary_i))
						+" ("+PdgToString(secondary_PDG_code_2.at(secondary_i))+")",v_warning,verbosity);
					continue;
				}
				nuclide_daughter_pdg.at(neutron_parent_loc) = secondary_PDG_code_2.at(secondary_i);
			}
		} // end if from ncapture
	} // end scan for gammas/daughter nuclides
	Log(toolName+" found "+toString(n_gammas)+" gammas", v_debug,verbosity);
	
	// ==========================
	// TRANSFER TO OUTPUT VECTORS
	// ==========================
	// ok now we can create the desired output structure
	std::map<int, int> nuclide_id; // map known nuclide indices (in secondaries array) to index in nuclide vector
	for(int neutron_i=0; neutron_i<neutron_start_energy.size(); ++neutron_i){
		Log(toolName+" adding neutron "+toString(neutron_i)+" to output",v_debug,verbosity);
		int parent_nuclide_index = neutron_parent_nuclide_indices.at(neutron_i);
		Log(toolName+" it has parent nuclide index "+toString(parent_nuclide_index),v_debug,verbosity);
		// see if we've already noted this nuclide
		if(nuclide_id.count(parent_nuclide_index)==0){
			// new nuclide, create it
			Log(toolName+" new nuclide",v_debug,verbosity);
			nuclide_id.emplace(parent_nuclide_index,out_nuclide_pdg.size()); // 🚩
			out_nuclide_pdg.push_back(parent_nuclide_pdg.at(neutron_i));
			out_nuclide_creation_pos.push_back(parent_nuclide_creation_pos.at(neutron_i));
			out_nuclide_decay_pos.push_back(neutron_start_pos.at(neutron_i));
			out_nuclide_daughter_pdg.push_back(nuclide_daughter_pdg.at(neutron_i));
			// 🚩 pdg -1 is dummy nuclide to hold neutrons that did not have a valid parent nuclide index
			// TODO use daughter nuclide (which should not be a dummy) - daughter->parent pdg map lookup?
			
			// make the corresponding vectors for neutrons and gammas associated with it
			out_neutron_start_pos.push_back(std::vector<TLorentzVector>{});
			out_neutron_end_pos.push_back(std::vector<TLorentzVector>{});
			out_neutron_start_energy.push_back(std::vector<double>{});
			out_neutron_end_energy.push_back(std::vector<double>{});
			out_neutron_end_process.push_back(std::vector<int>{});
			out_neutron_ndaughters.push_back(std::vector<int>{});
			out_gamma_energy.push_back(std::vector<std::vector<double> >{});
			out_gamma_time.push_back(std::vector<std::vector<double> >{});
		}
		Log(toolName+" nuclide ready, adding neutron",v_debug,verbosity);
		// convert nuclide index in the secondaries array to index within our vector of nuclides
		int parent_nuclide_loc = nuclide_id.at(parent_nuclide_index);
		Log(toolName+" adding this neutron to nuclide "+toString(parent_nuclide_loc),v_debug,verbosity);
		// add this neutron and its daughter gammas to that parent nuclide
		out_neutron_start_pos.at(parent_nuclide_loc).push_back(neutron_start_pos.at(neutron_i));
		out_neutron_end_pos.at(parent_nuclide_loc).push_back(neutron_end_pos.at(neutron_i));
		out_neutron_start_energy.at(parent_nuclide_loc).push_back(neutron_start_energy.at(neutron_i));
		out_neutron_end_energy.at(parent_nuclide_loc).push_back(neutron_end_energy.at(neutron_i));
		out_neutron_end_process.at(parent_nuclide_loc).push_back(neutron_end_process.at(neutron_i));
		out_neutron_ndaughters.at(parent_nuclide_loc).push_back(gamma_energy.at(neutron_i).size());
		out_gamma_energy.at(parent_nuclide_loc).push_back(gamma_energy.at(neutron_i));
		out_gamma_time.at(parent_nuclide_loc).push_back(gamma_time.at(neutron_i));
	}
	Log(toolName+" logged "+toString(out_nuclide_pdg.size())+" nuclides", v_debug,verbosity);
	
	// record all primaries...? do we need this info?
	for(int primary_i=0; primary_i<n_outgoing_primaries; ++primary_i){
		int primary_pdg_code = G3ParticleCodeToPdg(primary_G3_code.at(primary_i));
		out_primary_pdg.push_back(primary_pdg_code);
		double mom_sq = pow(primary_start_mom.at(primary_i),2.);
		double mass = PdgToMass(primary_pdg_code);
		mass = (mass==0) ? 1 : mass;
		out_primary_energy.push_back(mom_sq/(2.*mass));
		out_primary_start_pos.push_back(primary_vertex_tvector); // not sure about the validity of this
		out_primary_end_pos.push_back(TLorentzVector(0,0,0,0)); // need to get this from a daughter
	}
	
	return 1;
}

int TruthNeutronCaptures::ReadEntryNtuple(long entry_number){
	int bytesread = myTreeReader.GetEntry(entry_number);
	if(bytesread<=0) return bytesread;
	
	int success = 
	// file level
	// simulation version?
	(myTreeReader.GetBranchValue("wlen",water_transparency)) &&  // [cm]
	
	// event meta info
	(myTreeReader.GetBranchValue("nrun",run_number))         &&
	(myTreeReader.GetBranchValue("nsub",subrun_number))      &&
	(myTreeReader.GetBranchValue("nev",event_number))        &&
	(myTreeReader.GetBranchValue("nsube",subevent_number))   &&  // how does this relate to after trigger?
//	(myTreeReader.GetBranchValue("date",event_date))         &&  // [year,month,day]
//	(myTreeReader.GetBranchValue("time",time))               &&  // [hour,minute,second,?]
	
	// event level detector info
	(myTreeReader.GetBranchValue("nhit",N_hit_ID_PMTs))      &&  // "nqisk"
	(myTreeReader.GetBranchValue("potot",total_ID_pes))      &&  // "qismsk"
	(myTreeReader.GetBranchValue("pomax",max_ID_PMT_pes))    &&  // "qimxsk", presumably max # PEs from an ID PMT?
	
	// numnu is 0 even when npar is >3...
//	// neutrino interaction info - first primaries array includes neutrino and target (index 0 and 1)
//	(myTreeReader.GetBranchValue("mode",nu_intx_mode))       &&  // see neut_mode_to_string(mode)
//	(myTreeReader.GetBranchValue("numnu",tot_n_primaries))   &&  // both ingoing and outgoing
//	
//	// following are arrays of size numnu
//	(myTreeReader.GetBranchValue("ipnu",primary_pdg))        &&  // see constants::numnu_code_to_string
//	(myTreeReader.GetBranchValue("pnu",primary_momentum))    &&  // [GeV/c]
	
	// primary event - second primaries array includes more info
	(myTreeReader.GetBranchValue("posv",primary_event_vertex))          &&  // [cm]
//	(myTreeReader.GetBranchValue("wallv",primary_event_dist_from_wall)) &&  // [cm]
	(myTreeReader.GetBranchValue("npar",n_outgoing_primaries))          &&  // should be (tot_n_primaries - 2)?
	
	// following are arrays of size npar
	(myTreeReader.GetBranchValue("ipv",primary_G3_code))                &&  // see constants::g3_to_pdg
	(myTreeReader.GetBranchValue("pmomv",primary_start_mom))            &&  // [units?]
	
//	// secondaries - first secondaries arrays...
//	(myTreeReader.GetBranchValue("npar2",n_secondaries_1))              &&
	// npar2 is 0 even when nscndprt is not???
//	
//	// following are arrays of size npar2
//	(myTreeReader.GetBranchValue("ipv2",secondary_G3_code_1))                &&  // 
//	(myTreeReader.GetBranchValue("posv2",secondary_start_vertex_1))          &&  // [cm?] what about time?
//	(myTreeReader.GetBranchValue("wallv2",secondary_start_dist_from_wall_1)) &&  // [cm?]
//	(myTreeReader.GetBranchValue("pmomv2",secondary_start_mom_1))            &&  // [units?]
//	(myTreeReader.GetBranchValue("iorg",secondary_origin_1))                 &&  // what is "origin"?
	
	// secondaries - second secondaries array...
	(myTreeReader.GetBranchValue("nscndprt",n_secondaries_2))          &&
	
	// following are arrays of size nscndprt
	(myTreeReader.GetBranchValue("iprtscnd",secondary_PDG_code_2))     &&  //
	(myTreeReader.GetBranchValue("vtxscnd",secondary_start_vertex_2))  &&  // [units?]
	(myTreeReader.GetBranchValue("tscnd",secondary_start_time_2))      &&  // [ns]? relative to event start?
	(myTreeReader.GetBranchValue("pscnd",secondary_start_mom_2))       &&  // [units?]
	(myTreeReader.GetBranchValue("lmecscnd",secondary_gen_process))    &&  // constants::G3_process_code_to_string
	(myTreeReader.GetBranchValue("nchilds",secondary_n_daughters))     &&  // 
	(myTreeReader.GetBranchValue("iprntidx",parent_index))             &&  // if >0, 1-based index in this array
//	(myTreeReader.GetBranchValue("ichildidx",secondary_first_daugher_index)) &&  // if >0, 1-based index in this
	
	// further parentage information - still arrays of size nscndprt. Useful?
//	(myTreeReader.GetBranchValue("iprntprt",parent_G3_code))           &&  // or is it a PDG code?
	(myTreeReader.GetBranchValue("pprnt",parent_mom_at_sec_creation))  &&  // use w/γ to get n energy @ capture
	(myTreeReader.GetBranchValue("vtxprnt",parent_init_pos))           &&  // [cm?] parent pos @ birth
	(myTreeReader.GetBranchValue("pprntinit",parent_init_mom))         &&  // [MeV?] parent mom @ birth
//	(myTreeReader.GetBranchValue("itrkscnd",parent_G3_trackid))        &&  // how do we use this?
//	(myTreeReader.GetBranchValue("istakscnd",parent_G3_stack_trackid)) &&  // how do we use this?
	(myTreeReader.GetBranchValue("iprnttrk",parent_trackid));        //&&  // relates secondaries to primaries
	// NOTE this is carried over to daughters of secondaries, so only use as parent if iprntidx==0
//	(myTreeReader.GetBranchValue("iorgprt",parent_track_pid_code))     &&  // i'm so confused
	
	return success;
}

int TruthNeutronCaptures::CreateOutputFile(std::string filename){
	// create the output ROOT file and TTree for writing
	// =================================================
	outfile = new TFile(filename.c_str(), "RECREATE");
	outtree = new TTree("eventtree", "Events with Neutron Captures");
	
	// create branches
	// ---------------
	// file level
	outtree->Branch("filename",&out_filename);
//	outtree->Branch("skdetsim_version",&out_skdetsim_version);    // where?
//	outtree->Branch("tba_table_version",&out_tba_table_version);  // where?
	outtree->Branch("water_transparency",&out_water_transparency);
	
	// event level
//	outtree->Branch("run_number",&out_run_number);
//	outtree->Branch("subrun_number",&out_subrun_number);
	outtree->Branch("entry_number",&out_entry_number);
	outtree->Branch("subevent_num",&out_subevent_number);
	
	// primary particle
	outtree->Branch("primary_pdg",&out_primary_pdg);
	outtree->Branch("primary_energy",&out_primary_energy);
	outtree->Branch("primary_start_pos",&out_primary_start_pos);
	outtree->Branch("primary_end_pos",&out_primary_end_pos);
	
	// parent nuclide - one for each neutron
	outtree->Branch("nuclide_pdg",&out_nuclide_pdg);
	outtree->Branch("nuclide_creation_pos",&out_nuclide_creation_pos);
	outtree->Branch("nuclide_decay_pos",&out_nuclide_decay_pos);
	outtree->Branch("nuclide_daughter_pdg",&out_nuclide_daughter_pdg);
	
	// neutron
	outtree->Branch("neutron_start_pos",&out_neutron_start_pos);
	outtree->Branch("neutron_end_pos",&out_neutron_end_pos);
	outtree->Branch("neutron_start_energy",&out_neutron_start_energy);
	outtree->Branch("neutron_end_energy",&out_neutron_end_energy);
	outtree->Branch("neutron_end_process",&out_neutron_end_process);
	outtree->Branch("neutron_n_daughters",&out_neutron_ndaughters);
	
	// gamma
	outtree->Branch("gamma_energy",&out_gamma_energy);
	outtree->Branch("gamma_time",&out_gamma_time);
	
	return 1;
}

void TruthNeutronCaptures::ClearOutputTreeBranches(){
	// clear any vector branches
	
	out_primary_pdg.clear();
	out_primary_energy.clear();
	out_primary_start_pos.clear();
	out_primary_end_pos.clear();
	
	// parent nuclide
	out_nuclide_pdg.clear();
	out_nuclide_creation_pos.clear();
	out_nuclide_decay_pos.clear();
	out_nuclide_daughter_pdg.clear();
	
	// neutrons
	out_neutron_start_pos.clear();
	out_neutron_end_pos.clear();
	out_neutron_start_energy.clear();
	out_neutron_end_energy.clear();
	out_neutron_end_process.clear();
	out_neutron_ndaughters.clear();
	
	// gammas
	out_gamma_energy.clear();
	out_gamma_time.clear();
	
	return;
}

void TruthNeutronCaptures::PrintBranches(){
	std::cout<<"==========================================================="<<std::endl;
	std::cout<<"PRINTING EVENT"<<std::endl;
	std::cout<<"filename: "<<out_filename<<std::endl;
	std::cout<<"water transparency: "<<out_water_transparency<<std::endl;
	
//	std::cout<<"skdetsim_version: "<<out_skdetsim_version<<std::endl
//			 <<"tba_table_version: "<<out_tba_table_version<<std::endl
//			 <<"neutron_process_map: {";
//	for(std::map<int,std::string>::const_iterator aprocess=neutron_process_map.begin();
//	    aprocess!=neutron_process_map.end(); ++aprocess){
//			std::cout<<"["<<aprocess->first<<"="<<aprocess->second<<"], ";
//	}
//	std::cout<<"\b\b}"<<std::endl;
	std::cout<<"entry_number: "<<out_entry_number<<std::endl
			 <<"subevent_number:" <<out_subevent_number<<std::endl
			 <<"num primaries: "<<out_primary_pdg.size()<<std::endl;
	if(out_primary_pdg.size()){
		std::cout<<"primary vertex:"
				 <<" ("<<out_primary_start_pos.at(0).X()
				 <<", "<<out_primary_start_pos.at(0).Y()
				 <<", "<<out_primary_start_pos.at(0).Z()<<")"<<std::endl;
	}
	
	// print primaries
	for(int primary_i=0; primary_i<out_primary_pdg.size(); ++primary_i){
		std::cout<<"\tprimary ("<<primary_i<<"): "<<std::endl
				 <<"\t\tprimary pdg: "<<out_primary_pdg.at(primary_i)<<std::endl
				 <<"\t\tprimary energy: "<<out_primary_energy.at(primary_i)<<std::endl;
//		std::cout<<"\t\tprimary start pos:"
//				 <<" ("<<out_primary_start_pos.at(primary_i).X()
//				 <<", "<<out_primary_start_pos.at(primary_i).Y()
//				 <<", "<<out_primary_start_pos.at(primary_i).Z()<<")"<<std::endl
//				 <<"\t\tprimary end pos:"
//				 <<" ("<<out_primary_end_pos.at(primary_i).X()
//				 <<", "<<out_primary_end_pos.at(primary_i).Y()
//				 <<", "<<out_primary_end_pos.at(primary_i).Z()<<")"<<std::endl;
	}
	
	int total_neutrons=0;
	int total_gammas=0;
	// print nuclides
	std::cout<<"num nuclides: "<<out_nuclide_pdg.size()<<std::endl;
	for(int nuclide_i=0; nuclide_i<out_nuclide_pdg.size(); ++nuclide_i){
		std::cout<<"\tnuclide "<<nuclide_i<<": "<<std::endl
				 <<"\t\tnuclide pdg: "<<out_nuclide_pdg.at(nuclide_i)<<std::endl
				 <<"\t\tdaughter pdg: "<<out_nuclide_daughter_pdg.at(nuclide_i)<<std::endl
				 <<"\t\tnuclide start pos:"
				 <<" ("<<out_nuclide_creation_pos.at(nuclide_i).X()
				 <<", "<<out_nuclide_creation_pos.at(nuclide_i).Y()
				 <<", "<<out_nuclide_creation_pos.at(nuclide_i).Z()<<")"<<std::endl
				 <<"\t\tnuclide end pos:"
				 <<" ("<<out_nuclide_decay_pos.at(nuclide_i).X()
				 <<", "<<out_nuclide_decay_pos.at(nuclide_i).Y()
				 <<", "<<out_nuclide_decay_pos.at(nuclide_i).Z()<<")"<<std::endl;
		
		// each nuclide may have multiple neutrons
		// (every neutron should have a nuclide ... backgrounds?)
		std::cout<<"\t\tnum neutrons from this nuclide: "<<out_neutron_start_energy.at(nuclide_i).size()<<std::endl;
		for(int neutron_i=0; neutron_i<out_neutron_start_energy.at(nuclide_i).size(); ++neutron_i){
			std::cout<<"\t\tneutron "<<neutron_i<<": "<<std::endl
					 <<"\t\t\tneutron start energy: "
					 <<out_neutron_start_energy.at(nuclide_i).at(neutron_i)<<std::endl
					 <<"\t\t\tneutron start pos:"
					 <<" ("<<out_neutron_start_pos.at(nuclide_i).at(neutron_i).X()
					 <<", "<<out_neutron_start_pos.at(nuclide_i).at(neutron_i).Y()
					 <<", "<<out_neutron_start_pos.at(nuclide_i).at(neutron_i).Z()<<")"<<std::endl
					 <<"\t\t\tneutron end pos:"
					 <<" ("<<out_neutron_end_pos.at(nuclide_i).at(neutron_i).X()
					 <<", "<<out_neutron_end_pos.at(nuclide_i).at(neutron_i).Y()
					 <<", "<<out_neutron_end_pos.at(nuclide_i).at(neutron_i).Z()<<")"<<std::endl
					 <<"\t\t\tneutron end energy:"
					 <<out_neutron_end_energy.at(nuclide_i).at(neutron_i)<<std::endl
					 <<"\t\t\tneutron end process: "<<out_neutron_end_process.at(nuclide_i).at(neutron_i)
					 <<std::endl;
			total_neutrons++;
			
			// each neutron may have multiple decay gammas
			std::cout<<"\t\t\tnum gammas from this neutron: "
					 <<out_gamma_energy.at(nuclide_i).at(neutron_i).size()<<std::endl;
			for(int gamma_i=0; gamma_i<out_gamma_energy.at(nuclide_i).at(neutron_i).size(); ++gamma_i){
				std::cout<<"\t\t\tgamma "<<gamma_i<<": "<<std::endl
						 <<"\t\t\t\tgamma energy: "<<out_gamma_energy.at(nuclide_i).at(neutron_i).at(gamma_i)
						 <<std::endl
						 <<"\t\t\t\tgamma time: "<<out_gamma_time.at(nuclide_i).at(neutron_i).at(gamma_i)
						 <<std::endl;
				total_gammas++;
			}
		}
	}
	std::cout<<"total neutrons in event: "<<total_neutrons<<std::endl
			 <<"total gammas in event: "<<total_gammas<<std::endl;
	// shallow sanity checks
	std::cout<<"consistency check: num neutrons outer vs num nuclides: "
			 <<(out_neutron_start_energy.size()==out_nuclide_pdg.size())<<std::endl
			 <<"consistency check: num gammas outer vs num nuclides: "
			 <<(out_gamma_energy.size()==out_nuclide_pdg.size())<<std::endl;
	std::cout<<"==========================================================="<<std::endl;
}

int TruthNeutronCaptures::WriteTree(){
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

void TruthNeutronCaptures::CloseFile(){
	outtree->ResetBranchAddresses();
	outfile->Write("*",TObject::kOverwrite);
	outfile->Close();
	delete outfile;
	outfile = nullptr;
};


