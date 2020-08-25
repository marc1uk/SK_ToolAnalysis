/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef MTreeReader_H
#define MTreeReader_H

#include <string>
#include <map>
#include <utility> // pair

#include "basic_array.h"

class TFile;
class TChain;
class TTree;
class TBranch;
class TLeaf;

class MTreeReader {
	public:
	
	MTreeReader(std::string filename, std::string treename);
	MTreeReader(){};
	~MTreeReader();
	int Load(std::string filename, std::string treename);
	int LoadFile(std::string filename);
	int LoadTree(std::string treename);
	// TODO constructor/loader for tchains or tree pointers
	// TODO accept vector of files from LoadFileList Algorithms function
	
	// get a pointer to an object
	template<typename T>
	int GetBranchValue(std::string branchname, const T* &pointer_in){
		if(branch_value_pointers.count(branchname)==0){
			std::cerr<<"No such branch "<<branchname<<std::endl;
			std::cerr<<"known branches: {";
			for(auto&& abranch : branch_value_pointers) std::cout<<abranch.first<<", ";
			std::cerr<<"\b\b}"<<std::endl;
			return 0;
		}
		pointer_in = reinterpret_cast<const T*>(branch_value_pointers.at(branchname));
		if(verbosity>3) std::cout<<"retrieved pointer to "<<type_name<T>()<<" at "<<pointer_in<<std::endl;
		return 1;
	}
	
	// get the object itself - valid for primitives only (could add objects if we made copies)
	template<typename T>
	int GetBranchValue(std::string branchname, T& ref_in){
		// check we know this branch
		if(branch_value_pointers.count(branchname)==0){
			std::cerr<<"No such branch "<<branchname<<std::endl;
			return 0;
		}
		// check if the branch is a primitive
		if(branch_isobject.at(branchname)||branch_isarray.at(branchname)){
			std::cerr<<"Branch "<<branchname
				 <<" is not a primitive; please pass a suitable pointer"
				 <<" or basic_array to GetBranchValue()"<<std::endl;
			// TODO copy-construct an object, if they really want
			// requires a suitable copy constructor (or operator=) exists for the class
			return 0;
		}
		// else for primitives, de-reference the pointer to allow the user a copy
		T* objp = reinterpret_cast<T*>(branch_value_pointers.at(branchname));
		ref_in = *objp;
		return 1;
	}
	
	// specialization for arrays
	template<typename T>
	int GetBranchValue(std::string branchname, basic_array<T>& ref_in){
		// check we know this branch
		if(branch_value_pointers.count(branchname)==0){
			std::cerr<<"No such branch "<<branchname<<std::endl;
			return 0;
		}
		// check if the branch is an array - this template specialization is only for arrays
		if(not branch_isarray.at(branchname)){
			std::cerr<<"Branch "<<branchname
				 <<" is not an array; please check your datatype to GetBranchValue()"<<std::endl;
			return 0;
		}
		// for dynamic arrays we may need to update our pointer to the stored array
		UpdateBranchPointer(branchname);
		// next we need to know the array dimensions, which may vary by entry
		std::vector<size_t> branchdims = GetBranchDims(branchname);
		// finally construct and return the wrapper
		ref_in = basic_array<T>(branch_value_pointers.at(branchname),branchdims);
		return 1;
	}
	
	// misc operations
	void SetVerbosity(int verbin);
	
	// file/tree level getters
	TFile* GetFile();
	TTree* GetTree();
	
	// tree operations
	int Clear();
	void SetAutoClear(bool autoclearin);
	int GetEntry(long entry_number);
	long GetEntriesFast();
	long GetEntries();
	
	// maps of branch properties
	std::map<std::string,std::string> GetBranchTypes();
	
	// specific branch properties
	TBranch* GetBranch(std::string branchname);
	std::string GetBranchType(std::string branchname);
	std::vector<size_t> GetBranchDims(std::string branchname);
	
	private:
	// functions
	int ParseBranches();
	int ParseBranchDims(std::string branchname);
	int UpdateBranchPointer(std::string branchname);
	int UpdateBranchPointers();
	
	// variables
	std::map<std::string,TBranch*> branch_pointers;  // branch name to TBranch*
	std::map<std::string,bool> branch_istobject;     // branch inherits from TObject so has Clear method
	std::map<std::string,TLeaf*> leaf_pointers;      // branch name to TLeaf*
	std::map<std::string,std::string> branch_titles; // branch name to title (including type string)
	std::map<std::string,intptr_t> branch_value_pointers; // branch name to pointer to value, cast to intptr_t
	std::map<std::string,std::string> branch_types;  // branch name to string describing type - not good for arrays
	std::map<std::string,bool> branch_isobject;      // does branch hold an object
	std::map<std::string,bool> branch_isarray;       // does branch hold a (c-style) array
	std::map<std::string,std::vector<std::pair<std::string,int>>> branch_dimensions; // dims of variable size arrays
	std::map<std::string,std::vector<size_t>> branch_dims_cache; // dims of constant sized arrays
	
	// XXX TODO FIXME XXX we should hold unique_ptr to these, and then implement std::move
	// on copy construction and assignment, otherwise i think MTreeReader mytreeReader = MTreeReader(file, tree)
	// fails because it constructs a temporary, uses operator= to copy over the members (the pointers),
	// then call the destructor on the temporary which closes the file and deletes the trees!
	// for now workaround is to use myTreeReader.Load(file, tree) not operator=
	TFile* thefile=nullptr;
	TTree* thetree=nullptr;
	TChain* thechain=nullptr;
	bool autoclear=true;  // call 'Clear' method on all object branches before GetEntry
	int verbosity=1; // TODO add to constructor
};

/*
// mechanism to check if a given type has a Clear() method.
// ROOT seems to flag both objects and stl containers as inheriting from TObject
// with InheritsFrom("TLeafElement")
// from https://stackoverflow.com/a/29772824
template<typename T>
struct has_clear {
	// NOTE: sig_matches() must come before fn_exists() as it is used for its type.
	// Also, no function bodies are needed as they are never called.
	
	// This matching sig results in a return type of true_type
	template<typename Q>
	static auto sig_matches(void(Q::*)()) -> std::true_type;
	
	// If the member function Q::Clear exists and a sig_matches() function
	// exists with the required sig, then the return type is the return type of
	// sig_matches(), otherwise this function can't exist because at least one
	// the types don't exist so match against fn_exists(...).
	template <typename Q>
	static auto fn_exists(std::nullptr_t) -> decltype(sig_matches<Q>(&Q::Clear));
	
	// Member function either doesn't exist or doesn't match against a 
	// sig_matches() function.
	template<typename Q>
	static auto fn_exists(...) -> std::false_type;
	
	// Intermediate storage of type for clarity
	typedef decltype(fn_exists<T>(nullptr)) type;
	
	// Storing the resulting value
	static int const value = type::value;
};
*/


#endif // defined MTreeReader_H
