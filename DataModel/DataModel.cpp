#include "DataModel.h"

DataModel::DataModel(){ rootTApp = new TApplication("rootTApp",0,0); }
DataModel::~DataModel(){ if(rootTApp) delete rootTApp; }

/*
TTree* DataModel::GetTTree(std::string name){

  return m_trees[name];

}


void DataModel::AddTTree(std::string name,TTree *tree){

  m_trees[name]=tree;

}


void DataModel::DeleteTTree(std::string name){

  m_trees.erase(name);

}

*/

TApplication* DataModel::GetTApp(){
	if(rootTApp==nullptr){
		rootTApp = new TApplication("rootTApp",0,0);
	}
	return rootTApp;
}

bool DataModel::RegisterReader(std::string readerName, std::function<bool()> hasAFT, std::function<bool()> loadSHE, std::function<bool()> loadAFT, std::function<bool(int)> loadCommon){
	hasAFTs.emplace(readerName, hasAFT);
	loadSHEs.emplace(readerName, loadSHE);
	loadAFTs.emplace(readerName, loadAFT);
	loadCommons.emplace(readerName, loadCommon);
	return true;
}

bool DataModel::HasAFT(std::string ReaderName){
	//if(TreeReaders.count(ReaderName) && TreeReaders[ReaderName]){
	//	return TreeReaders[ReaderName]->HasAFT();
	if(hasAFTs.count(ReaderName)){
		return hasAFTs.at(ReaderName)();
	} else {
		std::cerr<<"HasAFT requested for Unknown reader "<<ReaderName<<std::endl;
	}
	return false;
}

bool DataModel::LoadSHE(std::string ReaderName){
	//if(TreeReaders.count(ReaderName) && TreeReaders[ReaderName]){
	//	return TreeReaders[ReaderName]->LoadSHE();
	if(loadSHEs.count(ReaderName)){
		return loadSHEs.at(ReaderName)();
	} else {
		std::cerr<<"LoadSHE requested for Unknown reader "<<ReaderName<<std::endl;
	}
	return false;
}

bool DataModel::LoadAFT(std::string ReaderName){
	//if(TreeReaders.count(ReaderName) && TreeReaders[ReaderName]){
	//	return TreeReaders[ReaderName]->LoadAFT();
	if(loadAFTs.count(ReaderName)){
		return loadAFTs.at(ReaderName)();
	} else {
		std::cerr<<"LoadAFT requested for Unknown reader "<<ReaderName<<std::endl;
	}
	return false;
}

bool DataModel::LoadEntry(std::string ReaderName, int entry_i){
	if(loadCommons.count(ReaderName)){
		return loadCommons.at(ReaderName)(entry_i);
	} else {
		std::cerr<<"LoadEntry requested for Unknown reader "<<ReaderName<<std::endl;
	}
	return false;
}
