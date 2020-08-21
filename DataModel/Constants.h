/* vim:set noexpandtab tabstop=4 wrap */
#ifndef Constants_h
#define Constants_h

#include <string>
#include <map>

#include "TDatabasePDG.h"

class TVector3;
class TLorentzVector;

std::string G3_process_code_to_string(int process_code);
std::string numnu_code_to_string(int numnu_code);
std::string neut_mode_to_string(int neut_code);
std::string PdgToString(int code);
int StringToPdg(std::string name);
int PdgToG3ParticleCode(int code);
int G3ParticleCodeToPdg(int code);
int StringToG3ParticleCode(std::string name);
std::string G3ParticleCodeToString(int code);
double PdgToMass(int code);

void PrintVector(TVector3& avec, bool newline=false);
void PrintVector(TLorentzVector& avec, bool newline=false);

namespace constants{
	static const TDatabasePDG particleDb;
	// https://root.cern.ch/root/html532/src/TDatabasePDG.cxx.html
	// https://root.cern/doc/v608/classTDatabasePDG.html
	
	static const std::map<int,std::string> G3_process_code_to_string{
		// TODO fill in the missing gaps?
		{5, "Decay"},
		{6, "Pair production"},
		{7, "Compton Scatter"},
		{8, "Photo-electric"},
		{9, "Bremsstrahlung"},
		{12, "Hadronic Interaction"},
		{13, "Hadronic Elastic Coherent Scattering"},
		{18, "Neutron Capture"},
		{20, "Hadronic Inelastic"},
		{21, "Muon-Nuclear Interaction"},
		{23, "Photonuclear"},
		{30, "Below tracking threshold"}
	};
	
	static const std::map<int,std::string> g3_particle_code_to_string{
		{1,"Gamma"},
		{2,"Positron"},
		{3,"Electron"},
		{4,"Neutrino"},
		{5,"Muon +"},
		{6,"Muon -"},
		{7,"Pion 0"},
		{8,"Pion +"},
		{9,"Pion -"},
		{10,"Kaon 0 Long"},
		{11,"Kaon +"},
		{12,"Kaon -"},
		{13,"Neutron"},
		{14,"Proton"},
		{15,"Antiproton"},
		{16,"Kaon 0 Short"},
		{17,"Eta"},
		{18,"Lambda"},
		{19,"Sigma +"},
		{20,"Sigma 0"},
		{21,"Sigma -"},
		{22,"Xi 0"},
		{23,"Xi -"},
		{24,"Omega -"},
		{25,"Antineutron"},
		{26,"Antilambda"},
		{27,"Antisigma -"},
		{28,"Antisigma 0"},
		{29,"Antisigma +"},
		{30,"Antixi 0"},
		{31,"Antixi +"},
		{32,"Antiomega +"},
		{45,"Deuteron"},
		{46,"Triton"},
		{47,"Alpha"},
		{48,"Geantino"},
		{49,"He3"},
		{50,"Cerenkov"}
	};
	
	static const std::map<std::string,int> string_to_g3_particle_code{
		{"Gamma",1},
		{"Positron",2},
		{"Electron",3},
		{"Neutrino",4},
		{"Muon +",5},
		{"Muon -",6},
		{"Pion 0",7},
		{"Pion +",8},
		{"Pion -",9},
		{"Kaon 0 Long",10},
		{"Kaon +",11},
		{"Kaon -",12},
		{"Neutron",13},
		{"Proton",14},
		{"Antiproton",15},
		{"Kaon 0 Short",16},
		{"Eta",17},
		{"Lambda",18},
		{"Sigma +",19},
		{"Sigma 0",20},
		{"Sigma -",21},
		{"Xi 0",22},
		{"Xi -",23},
		{"Omega -",24},
		{"Antineutron",25},
		{"Antilambda",26},
		{"Antisigma -",27},
		{"Antisigma 0",28},
		{"Antisigma +",29},
		{"Antixi 0",30},
		{"Antixi +",31},
		{"Antiomega +",32},
		{"Deuteron",45},
		{"Triton",46},
		{"Alpha",47},
		{"Geantino",48},
		{"He3",49},
		{"Cerenkov",50}
	};
	
	static const std::map<int,int> g3_particle_code_to_pdg{
		// from https://root.cern.ch/root/html532/src/TDatabasePDG.cxx.html#228
		// or use Int_t TDatabasePDG::ConvertGeant3ToPdg(Int_t Geant3number)
		{1,22},        // photon
		{25,-2112},    // anti-neutron
		{2,-11},       // e+
		{26,-3122},    // anti-Lambda
		{3,11},        // e-
		{27,-3222},    // Sigma-
		{4,12},        // e-neutrino : Geant3 just has "neutrino"... which to map it to?
		{28,-3212},    // Sigma0
		{5,-13},       // mu+
		{29,-3112},    // Sigma+ (PB)*/
		{6,13},        // mu-
		{30,-3322},    // Xi0
		{7,111},       // pi0
		{31,-3312},    // Xi+
		{8,211},       // pi+
		{32,-3334},    // Omega+ (PB)
		{9,-211},      // pi-
		{33,-15},      // tau+
		{10,130},      // K long
		{34,15},       // tau-
		{11,321},      // K+
		{35,411},      // D+
		{12,-321},     // K-
		{36,-411},     // D-
		{13,2112},     // n
		{37,421},      // D0
		{14,2212},     // p
		{38,-421},     // D0
		{15,-2212},    // anti-proton
		{39,431},      // Ds+
		{16,310},      // K short
		{40,-431},     // anti Ds-
		{17,221},      // eta
		{41,4122},     // Lamba_c+
		{18,3122},     // Lambda
		{42,24},       // W+
		{19,3222},     // Sigma+
		{43,-24},      // W-
		{20,3212},     // Sigma0
		{44,23},       // Z
		{21,3112},     // Sigma-
		{45,3329},        // deuteron
		{22,3322},     // Xi0
		{46,0},        // triton
		{23,3312},     // Xi-
		{47,0},        // alpha
		{24,3334},     // Omega- (PB)
		{48,0}         // Geantino
	};
	
	static const std::map<int,int> pdg_to_g3_particle_code{
		// from https://root.cern.ch/root/html532/src/TDatabasePDG.cxx.html#228
		// or use Int_t TDatabasePDG::ConvertPdgToGeant3(Int_t pdgNumber)
		{22,1},        // photon
		{-2112,25},    // anti-neutron
		{-11,2},       // e+
		{-3122,26},    // anti-Lambda
		{11,3},        // e-
		{-3222,27},    // Sigma-
		{12,4},        // e-neutrino : Geant3 just has "neutrino"... which to map it to?
		{-3212,28},    // Sigma0
		{-13,5},       // mu+
		{-3112,29},    // Sigma+ (PB)*/
		{13,6},        // mu-
		{-3322,30},    // Xi0
		{111,7},       // pi0
		{-3312,31},    // Xi+
		{211,8},       // pi+
		{-3334,32},    // Omega+ (PB)
		{-211,9},      // pi-
		{-15,33},      // tau+
		{130,10},      // K long
		{15,34},       // tau-
		{321,11},      // K+
		{411,35},      // D+
		{-321,12},     // K-
		{-411,36},     // D-
		{2112,13},     // n
		{421,37},      // D0
		{2212,14},     // p
		{-421,38},     // D0
		{-2212,15},    // anti-proton
		{431,39},      // Ds+
		{310,16},      // K short
		{-431,40},     // anti Ds-
		{221,17},      // eta
		{4122,41},     // Lamba_c+
		{3122,18},     // Lambda
		{24,42},       // W+
		{3222,19},     // Sigma+
		{-24,43},      // W-
		{3212,20},     // Sigma0
		{23,44},       // Z
		{3112,21},     // Sigma-
		{3329,45},        // deuteron
		{3322,22},     // Xi0
		{0,46},        // triton
		{3312,23},     // Xi-
		{0,47},        // alpha
		{3334,24},     // Omega- (PB)
		{0,48}         // Geantino
	};
	
	static const std::map<int,std::string> numnu_code_to_string{
		{1,"is_neutrino"},   // initial state (incident) neutrino
		{2,"is_target"},     // initial state (struck) target
		{3,"fs_lepton"},     // final state (outgoing) target
		{4,"fs_target"},     // final state (outgoing) lepton
		{5,"fs_other"}       // final state (outgoing) other particle (codes 5 and up)
	};
	
	static const std::map<int,std::string> neut_mode_to_string{
		{1, "CC quasi-elastic"},
		{11, "CC single pi from delta resonance"},
		{12, "CC single pi from delta resonance"},
		{13, "CC single pi from delta resonance"},
		{16, "CC coherent pi production"},
		{21, "CC multi pi production"},
		{27, "CC diffractive pion production"},
		{31, "NC single pi from delta resonance"},
		{32, "NC single pi from delta resonance"},
		{33, "NC single pi from delta resonance"},
		{34, "NC single pi from delta resonance"},
		{36, "NC coherent pi"},
		{41, "NC multi pi production"},
		{47, "NC diffractive pion production"},
		{51, "NC elastic"},
		{52, "NC elastic"}
	};
	
	static const std::map<int,std::string> pdg_to_string{
		// FIXME many of these are probably WRONG, update from TParticlePDG
		{2212,"Proton"},
		{-2212,"Anti Proton"},
		{11,"Electron"},
		{-11,"Positron"},
		{12,"Electron Neutrino"},
		{-12,"Anti Electron Neutrino"},
		{22,"Gamma"},
		{2112,"Neutron"},
		{-2112,"Anti Neutron"},
		{-13,"Muon+"},
		{13,"Muon-"},
		{130,"Kaonlong"},
		{211,"Pion+"},
		{-211,"Pion-"},
		{321,"Kaon+"},
		{-321,"Kaon-"},
		{3122,"Lambda"},
		{-3122,"Antilambda"},
		{310,"Kaonshort"},
		{3112,"Sigma-"},
		{3222,"Sigma+"},
		{3212,"Sigma0"},
		{111,"Pion0"},
		{311,"Kaon0"},
		{-311,"Antikaon0"},
		{14,"Muon Neutrino"},
		{-14,"Anti Muon Neutrino"},
		{-3222,"Anti Sigma-"},
		{-3212,"Anti Sigma0"},
		{-3112,"Anti Sigma+"},
		{3322,"Xsi0"},
		{-3322,"Anti Xsi0"},
		{3312,"Xsi-"},
		{-3312,"Xsi+"},
		{3334,"Omega-"},
		{-3334,"Omega+"},
		{-15,"Tau+"},
		{15,"Tau-"},
		{100,"OpticalPhoton"},
		{3328,"Alpha"},
		{3329,"Deuteron"},
		{3330,"Triton"},
		{3351,"Li7"},
		{3331,"C10"},
		{3345,"B11"},
		{3332,"C12"},
		{3350,"C13"},
		{3349,"N13"},
		{3340,"N14"},
		{3333,"N15"},
		{3334,"N16"},
		{3335,"O16"},
		{3346,"Al27"},
		{3341,"Fe54"},
		{3348,"Mn54"},
		{3342,"Mn55"},
		{3352,"Mn56"},
		{3343,"Fe56"},
		{3344,"Fe57"},
		{3347,"Fe58"},
		{3353,"Eu154"},
		{3336,"Gd158"},
		{3337,"Gd156"},
		{3338,"Gd157"},
		{3339,"Gd155"}
	};
	
	static const std::map<std::string,int> string_to_pdg{
	// FIXME many of these are probably WRONG, update from TParticlePDG
		{"Proton",2212},
		{"Anti Proton",-2212},
		{"Electron",11},
		{"Positron",-11},
		{"Electron Neutrino",12},
		{"Anti Electron Neutrino",-12},
		{"Gamma",22},
		{"Neutron",2112},
		{"Anti Neutron",-2112},
		{"Muon+",-13},
		{"Muon-",13},
		{"Kaonlong",130},
		{"Pion+",211},
		{"Pion-",-211},
		{"Kaon+",321},
		{"Kaon-",-321},
		{"Lambda",3122},
		{"Antilambda",-3122},
		{"Kaonshort",310},
		{"Sigma-",3112},
		{"Sigma+",3222},
		{"Sigma0",3212},
		{"Pion0",111},
		{"Kaon0",311},
		{"Antikaon0",-311},
		{"Muon Neutrino",14},
		{"Anti Muon Neutrino",-14},
		{"Anti Sigma-",-3222},
		{"Anti Sigma0",-3212},
		{"Anti Sigma+",-3112},
		{"Xsi0",3322},
		{"Anti Xsi0",-3322},
		{"Xsi-",3312},
		{"Xsi+",-3312},
		{"Omega-",3334},
		{"Omega+",-3334},
		{"Tau+",-15},
		{"Tau-",15},
		{"OpticalPhoton",100},
		{"Alpha",3328},
		{"Deuteron",3329},    // from WChSandbox ... is this custom?
		{"Deuterium",100045}, // according to sonia...
		{"Triton",3330},
		{"Li7",3351},
		{"C10",3331},
		{"B11",3345},
		{"C12",3332},
		{"C13",3350},
		{"N13",3349},
		{"N14",3340},
		{"N15",3333},
		{"N16",3334},
		{"O16",3335},
		{"Al27",3346},
		{"Fe54",3341},
		{"Mn54",3348},
		{"Mn55",3342},
		{"Mn56",3352},
		{"Fe56",3343},
		{"Fe57",3344},
		{"Fe58",3347},
		{"Eu154",3353},
		{"Gd158",3336},
		{"Gd156",3337},
		{"Gd157",3338},
		{"Gd155",3339}
	};
	
}

#endif // define Constants_h
