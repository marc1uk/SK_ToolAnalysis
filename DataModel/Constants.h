/* vim:set noexpandtab tabstop=4 wrap */

#include <string>
#include <map>

std::string G3_process_code_to_string(int process_code);
std::string numnu_code_to_string(int numnu_code);
std::string neut_mode_to_string(int neut_code);

namespace constants{
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
	
}
