#include "Factory.h"

Tool* Factory(std::string tool){
Tool* ret=0;

// if (tool=="Type") tool=new Type;
if (tool=="DummyTool") ret=new DummyTool;

if (tool=="TruthNeutronCaptures") ret=new TruthNeutronCaptures;
if (tool=="LoadFileList") ret=new LoadFileList;
if (tool=="RootReadTest") ret=new RootReadTest;
if (tool=="PlotNeutronCaptures") ret=new PlotNeutronCaptures;
return ret;
}

