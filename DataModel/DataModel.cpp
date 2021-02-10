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

