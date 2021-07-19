#ifndef DATAMODEL_H
#define DATAMODEL_H

#include <map>
#include <unordered_map>
#include <string>
#include <vector>
#include <functional>

//#include "TTree.h"
#include "TApplication.h"

#include "Store.h"
#include "BoostStore.h"
#include "Logging.h"
#include "Utilities.h"

class MTreeReader;
class MTreeSelection;
class TreeReader;

#include <zmq.hpp>

/**
* \class DataModel
 *
 * This class Is a transient data model class for your Tools within the ToolChain. If Tools need to comunicate they pass all data objects through the data model. There fore inter tool data objects should be deffined in this class.
 *
 *
 * $Author: B.Richards $
 * $Date: 2019/05/26 18:34:00 $
 * Contact: b.richards@qmul.ac.uk
 *          
 */

class DataModel {


 public:
  
  DataModel(); ///< Simple constructor
  ~DataModel(); ///< Simple destructor

  //TTree* GetTTree(std::string name);
  //void AddTTree(std::string name,TTree *tree);
  //void DeleteTTree(std::string name);
  TApplication* GetTApp();

  Store vars; ///< This Store can be used for any variables. It is an inefficent ascii based storage    
  BoostStore CStore; ///< This is a more efficent binary BoostStore that can be used to store a dynamic set of inter Tool variables.
  std::map<std::string,BoostStore*> Stores; ///< This is a map of named BooStore pointers which can be deffined to hold a nammed collection of any tipe of BoostStore. It is usefull to store data that needs subdividing into differnt stores.
  std::map<std::string,MTreeReader*> Trees; ///< A map of MTreeReader pointers, used to read ROOT trees
  std::map<std::string,MTreeSelection*> Selectors; ///< A map of MTreeSelection pointers used to read event selections
//  std::map<std::string,TreeReader*> TreeReaders; ///< A map of TreeReader tool pointers, used to invoke LoadSHE/AFT
  std::unordered_map<std::string, std::function<bool()>> hasAFTs;
  std::unordered_map<std::string, std::function<bool()>> loadSHEs;
  std::unordered_map<std::string, std::function<bool()>> loadAFTs;
  std::unordered_map<std::string, std::function<bool(int)>> loadCommons;
  
  Logging *Log; ///< Log class pointer for use in Tools, it can be used to send messages which can have multiple error levels and destination end points  

  zmq::context_t* context; ///< ZMQ contex used for producing zmq sockets for inter thread,  process, or computer communication

  // These call the corresponding TreeReader functions.
  // The TreeReader instance is obtained from the name specified in their config file.
  bool HasAFT(std::string ReaderName);
  bool LoadSHE(std::string ReaderName);
  bool LoadAFT(std::string ReaderName);
  bool LoadEntry(std::string ReaderName, int entry_i);
  
  bool RegisterReader(std::string readerName, std::function<bool()> hasAFT, std::function<bool()> loadSHE, std::function<bool()> loadAFT, std::function<bool(int)> loadCommon);


 private:


  
  //std::map<std::string,TTree*> m_trees; 
  TApplication* rootTApp=nullptr;
  
  
  
};



#endif
