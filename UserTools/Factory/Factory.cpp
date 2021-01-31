#include "Factory.h"

Tool* Factory(std::string tool){
Tool* ret=0;

// if (tool=="Type") tool=new Type;
if (tool=="DummyTool") ret=new DummyTool;

if (tool=="TruthNeutronCaptures") ret=new TruthNeutronCaptures;
if (tool=="TruthNeutronCaptures_v2") ret=new TruthNeutronCaptures_v2;
if (tool=="TruthNeutronCaptures_v3") ret=new TruthNeutronCaptures_v3;
if (tool=="LoadFileList") ret=new LoadFileList;
if (tool=="RootReadTest") ret=new RootReadTest;
if (tool=="PlotNeutronCaptures") ret=new PlotNeutronCaptures;
if (tool=="PurewaterLi9Rate") ret=new PurewaterLi9Rate;
if (tool=="GracefulStop") ret=new GracefulStop;
if (tool=="PurewaterLi9Plots") ret=new PurewaterLi9Plots;
if (tool=="LoadBetaSpectraFluka") ret=new LoadBetaSpectraFluka;
return ret;
}

