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
if (tool=="GracefulStop") ret=new GracefulStop;
if (tool=="LoadBetaSpectraFluka") ret=new LoadBetaSpectraFluka;
if (tool=="PlotMuonDtDlt") ret=new PlotMuonDtDlt;
if (tool=="FitLi9Lifetime") ret=new FitLi9Lifetime;
if (tool=="FitPurewaterLi9NcaptureDt") ret=new FitPurewaterLi9NcaptureDt;
if (tool=="FitSpallationDt") ret=new FitSpallationDt;
if (tool=="PurewaterSpallAbundanceCuts") ret=new PurewaterSpallAbundanceCuts;
if (tool=="TreeReader") ret=new TreeReader;
if (tool=="lf_allfit") ret=new lf_allfit;
return ret;
}

