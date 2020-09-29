/* vim:set noexpandtab tabstop=4 wrap */
#include "MTreeReader.h"

#include "TROOT.h"
#include "TInterpreter.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TLeafElement.h"
//#include "TParameter.h"

#include "type_name_as_string.h"

#include <iostream>
#include <vector>
#include <sstream>
#include <algorithm> // find

// TODO constructor/loader for tchains or tree pointers
MTreeReader::MTreeReader(std::string filename, std::string treename){
	Load(filename, treename);
}

int MTreeReader::Load(std::string filename, std::string treename){
	int ok = LoadFile(filename);
	if(not ok) return ok;
	LoadTree(treename);
	if(not ok) return ok;
	ok = ParseBranches();
	return ok;
}

int MTreeReader::LoadFile(std::string filename){
	if(verbosity) std::cout<<"getting file"<<std::endl;
	thefile = TFile::Open(filename.c_str());
	if(not thefile){
		std::cerr<<"MTreeReader failed to load file "<<filename<<std::endl;
		return 0;
	}
	return 1;
}

int MTreeReader::LoadTree(std::string treename){
	if(verbosity) std::cout<<"getting tree "<<treename<<std::endl;
	thetree = (TTree*)thefile->Get(treename.c_str());
	if(not thetree){
		std::cerr<<"MTreeReader could not find tree "<<treename<<" in file "
				 <<thetree->GetCurrentFile()->GetName()<<std::endl;
		return 0;
	}
	// check the tree has entries
	if(verbosity) std::cout<<"getting num entries"<<std::endl;
	auto numentries = thetree->GetEntriesFast();
	if(numentries==0){
		std::cerr<<"MTreeReader tree "<<thetree->GetName()<<" in file "<<thetree->GetCurrentFile()->GetName()
				 <<" has no entries"<<std::endl;
		return 0;
	}
	if(verbosity) std::cout<<"tree has "<<numentries<<" entries"<<std::endl;
	return 1;
}

int MTreeReader::ParseBranches(){
	
	// preload first entry to get branches
	if(verbosity) std::cout<<"getting entry 0"<<std::endl;
	thetree->GetEntry(0);
	
	if(verbosity) std::cout<<"getting leaves"<<std::endl;
	for(int i=0; i<thetree->GetListOfLeaves()->GetEntriesFast(); ++i){
		TLeaf* lf=(TLeaf*)thetree->GetListOfLeaves()->At(i);
		std::string branchname = lf->GetName();
		leaf_pointers.emplace(branchname,lf);
		branch_pointers.emplace(branchname,lf->GetBranch());
		// the following is fine for objects, primitives or containers
		// but only returns the primitive type for c-style arrays
		branch_types.emplace(branchname,lf->GetTypeName());
		// the branch title includes dimensionality for c-style arrays
		// e.g. "mybranchname     int[nparts][3]/F"
		std::string branchtitle=lf->GetBranch()->GetTitle();
		branch_titles.emplace(branchname,branchtitle);
		
		// handle object pointers
		if (lf->IsA() == TLeafElement::Class()) {
			// could be TObjects, or could be stl containers
			// is intptr_t any better than (void*)? probably not.
			if(verbosity>1) std::cout<<"branch "<<branchname
				<<" stores object at "<<lf->GetValuePointer()<<std::endl;
			TBranchElement* bev = (TBranchElement*)lf->GetBranch();
			intptr_t objpp=reinterpret_cast<intptr_t>(bev->GetObject());
			branch_value_pointers.emplace(branchname,objpp);
			branch_isobject.emplace(branchname,true);
			branch_isarray.emplace(branchname,false);
			
			// both classes inheriting from TObject and STL containers
			// are flagged as TLeafElements, so need further check
			TClass* ac = TClass::GetClass(lf->GetTypeName());
			branch_istobject.emplace(branchname,ac->InheritsFrom("TObject"));
		}
		// handle arrays
		//else if(lf->GetLen()>1){  // flattened length. Unsuitable when dynamic size happens to be 1!
		else if(branchtitle.find_first_of("[",0)!=std::string::npos){  // hope for no '[' in branch names
			// we'll need to parse the title to retrieve the actual dimensions
			intptr_t objpp=reinterpret_cast<intptr_t>(lf->GetValuePointer());
			branch_value_pointers.emplace(branchname,objpp);
			branch_isobject.emplace(branchname,false);
			branch_isarray.emplace(branchname,true);
			branch_istobject.emplace(branchname,false);
		}
		// handle basic types
		else {
			intptr_t objpp=reinterpret_cast<intptr_t>(lf->GetValuePointer());
			branch_value_pointers.emplace(branchname,objpp);
			branch_isobject.emplace(branchname,false);
			branch_isarray.emplace(branchname,false);
			branch_istobject.emplace(branchname,false);
		}
	}
	
	// for all array branches, parse their titles to extract information about dimensions
	for(auto&& abranch : branch_isarray){
		if(abranch.second) ParseBranchDims(abranch.first);
	}
	return 1;
}

int MTreeReader::UpdateBranchPointer(std::string branchname){
	if(leaf_pointers.count(branchname)==0) return 0;
	TLeaf* lf = leaf_pointers.at(branchname);
	intptr_t objpp;
	if(branch_isobject.at(branchname)){
		// objects need another level of indirection
		TBranchElement* bev = (TBranchElement*)lf->GetBranch();
		objpp=reinterpret_cast<intptr_t>(bev->GetObject());
	} else {
		objpp=reinterpret_cast<intptr_t>(lf->GetValuePointer());
	}
	branch_value_pointers.at(branchname)=objpp;
	return 1;
}

int MTreeReader::UpdateBranchPointers(){
	// assume we only need to re-check dynamic arrays? Objects won't move, right?
	// dynamically sized arrays may do, if more space is required...
	for(auto&& abranch : branch_isarray){
		if(abranch.second){
			// fixed sized arrays won't need to be reallocated
			// so only update pointers if we don't have a cached (constant) size
			if(branch_dims_cache.count(abranch.first)==0){
				int ok = UpdateBranchPointer(abranch.first);
				if(not ok) return ok;
			}
		}
	}
	return 1;
}

int MTreeReader::ParseBranchDims(std::string branchname){
	// parse branch title for sequences of type '[X]' suggesting an array
	// extract 'X'. Scan the list of branch names for 'X', in which case
	// this is a variable length array, otherwise it's a fixed size so use stoi.
	std::string branchtitle = branch_titles.at(branchname);
	size_t startpos = 0;
	size_t endpos = 0;
	// each subsequent entry of the array represents a dimension
	// the pair will either hold a branch name (in the first element)
	// or a numeric size (in the second element)
	std::vector<std::pair<std::string,int>> this_branch_dimensions;
	while(true){
		startpos = branchtitle.find_first_of("[", endpos+1);
		if(startpos==std::string::npos) break;
		endpos = branchtitle.find_first_of("]",startpos+1);
		if(endpos==std::string::npos){
			std::cerr<<"Trailing '[' character while parsing branch titles for branch "
					 <<branchname<<", title "<<branchtitle<<"!"<<std::endl;
			return 0;
		}
		std::string sizestring = branchtitle.substr(startpos+1,endpos-startpos-1);
		//std::cout<<"extracted array size label "<<sizestring<<" for branch "<<branchname
		//		 <<" dimension "<<this_branch_dimensions.size()<<std::endl;
		// check if this string is the name of another branch - i.e. variable size array
		if(branch_pointers.count(sizestring)){
			// it is - we'll need to retrieve the corresponding branch entry to get the size on each entry
			this_branch_dimensions.emplace_back(std::pair<std::string, int>{sizestring,0});
		} else {
			// it ought to be a static size, try to convert from string to int
			// catch the exception thrown in the event that it can't be converted
			try{
				this_branch_dimensions.emplace_back(std::pair<std::string,int>{"",std::stoi(sizestring)});
			}
			catch(const std::invalid_argument& ia) {
				std::cerr<<"Failed to extract array dimension from branch title "<<branchtitle
						 <<" - not a recognised branchname and string to int conversion failed with "
						 <<ia.what()<<std::endl;
				return 0;
			}
			catch(const std::out_of_range& oor){
				std::cerr<<"Failed to extract array dimension from branch title "<<branchtitle
						 <<" - not a recognised branchname and string to int conversion failed with "
						 <<oor.what()<<std::endl;  // occurs if string exceeds range representable by int
				return 0;
			}
		}
		// loop back round for any further dimensions
	}
	branch_dimensions.emplace(branchname,this_branch_dimensions);
	if(branch_dimensions.at(branchname).size()==0){
		std::cerr<<"Failed to identify any dimensions for branch "<<branchname
				 <<" despite TLeaf::GetLength() returning >1"<<std::endl;
		return 0;
	}
	
	int vlevel=2;
	if(verbosity>vlevel) std::cout<<"end of branch title parsing, found "
								  <<branch_dimensions.at(branchname).size()<<" dimensions, [";
	// loop over the vector of dimensions
	std::string dims_string = "";
	bool allstatics=true;
	std::vector<size_t> cached_dims;
	for(auto&& dims_pair : branch_dimensions.at(branchname)){
		if(dims_pair.first==""){
			if(verbosity>vlevel) std::cout<<"(N)"<<dims_pair.second;  // static numeric size
			cached_dims.push_back(dims_pair.second);
			dims_string.append(std::string("[") + std::to_string(dims_pair.second) + std::string("]"));
		} else {
			if(verbosity>vlevel) std::cout<<"(B)"<<dims_pair.first; // name of branch that stores variable size
			dims_string.append(std::string("[") + dims_pair.first + std::string("]"));
			allstatics=false;
		}
		if((verbosity>vlevel)&&(dims_pair!=(branch_dimensions.at(branchname).back()))) std::cout<<"], [";
	}
	// if all dimensions are static we can cache the results for quicker lookup
	if(allstatics){ branch_dims_cache.emplace(branchname,cached_dims); }
	// append dimensions to the type string, since they aren't properly indicated by TLeaf::GetTypeName
	if(branch_types.count(branchname)){
		branch_types.at(branchname).append(dims_string);
	}
	
	return 1;
}

std::vector<size_t> MTreeReader::GetBranchDims(std::string branchname){
	// get the dimensions of the array for this entry
	// if all dimensions are constant we should have them cached
	if(branch_dims_cache.count(branchname)) return branch_dims_cache.at(branchname);
	
	// otherwise we must retrieve at least one size from another branch
	std::vector<size_t> dimstemp;  // dimensions for this entry
	if(branch_dimensions.count(branchname)==0){
		std::cerr<<"GetBranchDims called but no dimensions for this branch!"<<std::endl;
		return dimstemp;
	}
	// loop over dimensions
	for(auto&& adim : branch_dimensions.at(branchname)){
		if(adim.first==""){
			// this dimension is constant
			dimstemp.push_back(adim.second);
		} else {
			// this dimensions is a branch name - get the entry value of that branch
			std::string sizebranchname = adim.first;
			int lengththisentry;
			int get_ok = GetBranchValue(sizebranchname, lengththisentry);
			if(not get_ok){
				std::cerr<<"Failed to retrieve value for branch "<<sizebranchname
						 <<" required while obtaining this entry's dimensions for array in branch "
						 <<branchname<<std::endl;
				return dimstemp; // TODO throw an exception? We should do more than just print an error...
			}
			dimstemp.push_back(lengththisentry);
		}
	}
	return dimstemp;
}

int MTreeReader::Clear(){
	// loop over all branches
	for(auto&& isobject : branch_istobject){
		// skip if doesn't inherit from TObject so may not have Clear() method
		// XXX note, maybe we should check if it has a 'clear' method (stl container)
		// and invoke that if not? Should be safe even without doing that though.
		if(not isobject.second) continue;
		// get pointer to the object otherwise
		TObject* theobject = reinterpret_cast<TObject*>(branch_value_pointers.at(isobject.first));
		if(not theobject){
			// no object... is this an error?
			std::cerr<<"MTreeReader AutoClear error: failure to get pointer to TObject "
					 <<"for branch "<<isobject.second<<std::endl;
			continue;  // TODO throw suitable exception
		}
		theobject->Clear();
	}
	return 1; // TODO check for errs
}

int MTreeReader::GetEntry(long entry_number){
	// if we've been requested to invoke Clear() on all objects before each Get, do so
	if(verbosity>3) std::cout<<"MTreeReader GetEntry "<<entry_number<<std::endl;
	if(autoclear){
		int clear_ok = Clear();
		if(not clear_ok){ return -2; }
	}
	
	// load data from tree
	// The function returns the number of bytes read from the input buffer.
	// If entry does not exist the function returns 0. If an I/O error occurs, the function returns -1.
	auto bytesread = thetree->GetEntry(entry_number);
	return bytesread;
}

long MTreeReader::GetEntriesFast(){
	return thetree->GetEntriesFast();  // FIXME for TChain
}

long MTreeReader::GetEntries(){
	return thetree->GetEntries();
}

MTreeReader::~MTreeReader(){
	if(thechain) thechain->ResetBranchAddresses();  // are these mutually exclusive?
	if(thetree) thetree->ResetBranchAddresses();    // 
	if(thefile) thefile->Close();
	delete thefile;
}

// misc operations
void MTreeReader::SetVerbosity(int verbin){
	verbosity=verbin;
}

void MTreeReader::SetAutoClear(bool autoclearin){
	autoclear=autoclearin;
}

// file/tree level getters
TFile* MTreeReader::GetFile(){
	return thefile;
}

TTree* MTreeReader::GetTree(){
	return thetree;
}

// branch map getters
std::map<std::string,std::string> MTreeReader::GetBranchTypes(){
	return branch_types;
}

// specific branch getters
TBranch* MTreeReader::GetBranch(std::string branchname){
	if(branch_pointers.count(branchname)){
		return branch_pointers.at(branchname);
	} else {
		std::cerr<<"No such branch "<<branchname<<std::endl;
		return nullptr;
	}
}

std::string MTreeReader::GetBranchType(std::string branchname){
	if(branch_types.count(branchname)){
		return branch_types.at(branchname);
	} else {
		std::cerr<<"No such branch "<<branchname<<std::endl;
		return "";
	}
	return ""; // dummy to silence warning
}

int MTreeReader::DisableBranches(std::vector<std::string> branchnames){
	int success=1;
	// disable branches by name
	for(auto&& branchname : branchnames){
		if(branch_pointers.count(branchname)){
			branch_pointers.at(branchname)->SetStatus(0);
		} else {
			std::cerr<<"No such branch "<<branchname<<std::endl;
			success=0;
		}
	}
	return success;
}

int MTreeReader::EnableBranches(std::vector<std::string> branchnames){
	int success=1;
	// disable branches by name
	for(auto&& branchname : branchnames){
		if(branch_pointers.count(branchname)){
			branch_pointers.at(branchname)->SetStatus(1);
		} else {
			std::cerr<<"No such branch "<<branchname<<std::endl;
			success=0;
		}
	}
	return success;
}

int MTreeReader::OnlyDisableBranches(std::vector<std::string> branchnames){
	// enable all branches except those named
	int num_named_branches=branchnames.size();
	for(auto&& abranch : branch_pointers){
		if(std::find(branchnames.begin(),branchnames.end(),abranch.first)!=branchnames.end()){
			abranch.second->SetStatus(0);
			--num_named_branches;
		} else {
			abranch.second->SetStatus(1);
		}
	}
	// return whether we found all branches in the list given
	return (num_named_branches==0);
}

int MTreeReader::OnlyEnableBranches(std::vector<std::string> branchnames){
	// disable all branches except those named
	int num_named_branches=branchnames.size();
	for(auto&& abranch : branch_pointers){
		if(std::find(branchnames.begin(),branchnames.end(),abranch.first)!=branchnames.end()){
			abranch.second->SetStatus(1);
			--num_named_branches;
		} else {
			abranch.second->SetStatus(0);
		}
	}
	// return whether we found all branches in the list given
	return (num_named_branches==0);
}

