#include "evDisp.h"

#include "TH2.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TView.h"

#include "skroot.h"
#include "fortran_routines.h"

#include "ConnectionTable.h"
#include "TableReader.h"
#include "TableEntry.h"

evDisp::evDisp():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

bool evDisp::Initialise(std::string configfile, DataModel &data){

	if(configfile!="") m_variables.Initialise(configfile);
	//m_variables.Print();
	m_data= &data;
	
	m_variables.Get("verbosity",verbosity);
	m_variables.Get("treeReaderName",treeReaderName);
	m_variables.Get("plotVar",plotVar);
	m_variables.Get("dataSrc",dataSrc);
	m_variables.Get("evtSrc",evtSrc);
	m_variables.Get("plotStyle",plotStyle);
	
	if(m_data->Trees.count(treeReaderName)==0){
		Log("Failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,verbosity);
		return false;
	}
	myTreeReader = m_data->Trees.at(treeReaderName);
	
	// class to get PMT positions from cable number
	myConnectionTable = new ConnectionTable();
	
	if(plotStyle==0){
		// polymarker versions
		topCapHitMap = new TGraph2D();
		bottomCapHitMap = new TGraph2D();
		barrelHitMap = new TGraph2D();
	} else if(plotStyle==1){
		// histgrams for plotting
		topCapHeatMap = new TH2D("topCapHeatMap", "topCapHeatMap", 50, -2000, 2000, 50, -2000, 2000);
		bottomCapHeatMap = new TH2D("bottomCapHeatMap", "bottomCapHeatMap", 50, -2000, 2000, 50, -2000, 2000);
		barrelHeatMap = new TH2D("barrelSideHeatMap", "barrelSideHeatMap", 150, 6000, 6000, 50, -2000, 2000);
	}
	
	// canvas for plotting
	displayCanvas = new TCanvas("displayCanvas", "displayCanvas", 2000, 2000);
	/* divide up the canvas:
	 __________
	| ———  ——— |
	||TOP||BOT||
	| ———  ——— |
	| ———————— |
	|| BARREL ||
	| ———————— |
	 ‾‾‾‾‾‾‾‾‾‾
	*/
	gStyle->SetPalette(77);
	displayCanvas->Divide(1,2);
	displayPad = displayCanvas->cd(1);
	displayPad->Divide(2,1);
	
	return true;
}


bool evDisp::Execute(){
	
	// not sure where we want to get our hits:
	// sktqz_, rawtqinfo_, TQREAL, or somewhere else?
	
	/* rfm data files don't have TQREAL populated, even for events with a valid HEADER.
	   Instead they have a valid TQLIST branch, which is a TClonesArray of TQ class objects.
	   This gets read out by skrawread which calls TQRAWSK, which calls TQSKZ, which calls skroot_get_idtq
	   to read the TQLIST branch. The returned data gets passed into to TQSKZ, which puts it into the 
	   common block array IQISKZ. After returning to TQRAWSK, this calls SKTQCONV_QB to do charge conversion
	   from TDC/ADC counts to ns/p.e., and puts the result into TBUF_RAW (or QBUF_RAW, there are basically
	   arrays of T, Q, and cable number here). Finally, skroot_set_tree called by lowfit calls 
	   skroot_set_tqreal to transfer the converted hits from QBUF_RAW into the TQREAL branch.
	   The upshot is TQLIST is used to populate TQREAL after conversion from TDC/ADC counts to ns/p.e.
	   
	   Once populated TQREAL has members `int nhits`, and a `std::vector<float> T, Q` for hit times and charges
	   and a `std::vector<int> cables` for cable numbers of the hits.
	   
	   sktqz_ contains members `int nqiskz`, `int icabiz[]`, `float tiskz[]` and `float qiskz[]` among others
	*/
	
	if(dataSrc==0){
		// if we have both SHE and AFT available, load the appropriate common block data
		if(evtSrc==0) m_data->LoadSHE(treeReaderName);
		else          m_data->LoadAFT(treeReaderName);
	}
	if(dataSrc==1) GetData();  // get data from TreeReader if not using SK common blocks
	
	switch (dataSrc){
		case 0: {
			// sktqz_ common block
			totalPMTsActivated = sktqz_.nqiskz;
			break;
		}
		case 1: {
			// TQReal branch
			totalPMTsActivated = myTQReal->cables.size();
			break;
		}
		default: {
			// unknown
			Log(toolName+" unknown dataSrc: "+std::to_string(dataSrc),v_error,verbosity);
			totalPMTsActivated = 0;
			break;
		}
	}
	
	// resize internal arrays of TGraph2Ds
	if(topCapHitMap)    topCapHitMap->Set(totalPMTsActivated);
	if(bottomCapHitMap) bottomCapHitMap->Set(totalPMTsActivated);
	if(barrelHitMap)    barrelHitMap->Set(totalPMTsActivated);
	
	for (int pmtNumber = 0; pmtNumber < totalPMTsActivated; ++pmtNumber){
		switch (dataSrc){
			case 0: {
				// sktqz_ common block
				cableNumber = sktqz_.icabiz[pmtNumber];
				charge = sktqz_.qiskz[pmtNumber];
				time = sktqz_.tiskz[pmtNumber];
				break;
			}
			case 1: {
				// TQReal branch
				cableNumber = myTQReal->cables.at(pmtNumber);
				charge = myTQReal->Q.at(pmtNumber);
				time = myTQReal->T.at(pmtNumber);
				break;
			}
			default: {
				// unknown
				Log(toolName+" unknown dataSrc: "+std::to_string(dataSrc),v_error,verbosity);
				totalPMTsActivated = 0;
				break;
			}
		}
		
		
		//std::cout << "cable number is: " << cableNumber << std::endl;
		/*
		 41 # -- for inner-PMTs (1-11146)
		 42 # -- for muon VETO(11151,11152,11153,11154)
		 43 # -- for calibration ID (11155-?, see skveto.h for details)
		 44 # -- for muon chamber(only hut3 and hut4)
		 45 # -- for trigger ID QB (15001-15240, see skhead.h for details)
		 46 # -- for anti-PMT(20001-21885)
		*/
		
		// get tube position
		myConnectionTable->GetTubePosition(cableNumber, tubePosition);
		tubeRadialCoordinate = sqrt(pow(tubePosition[0], 2.f) + pow(tubePosition[1],2.f));
		tubeAngularCoordinate = acos(tubePosition[0] / tubeRadialCoordinate);
		if(tubePosition[1] > 0) tubeAngularCoordinate *= -1.;
		
		// for now we only plot one thing
		var = (plotVar) ? time : charge;
		
		// identify top/bottom cap PMTs by their z position
		if (tubePosition[2] == zMax){
			if(topCapHeatMap) topCapHeatMap->Fill(tubePosition[0], tubePosition[1], var);
			if(topCapHitMap)  topCapHitMap->SetPoint(pmtNumber, tubePosition[0], tubePosition[1], var);
		} else if (tubePosition[2] == zMin){
			if(bottomCapHeatMap) bottomCapHeatMap->Fill(tubePosition[0], tubePosition[1], var);
			if(bottomCapHitMap)  bottomCapHitMap->SetPoint(pmtNumber, tubePosition[0], tubePosition[1], var);
		} else {
			if(barrelHeatMap) barrelHeatMap->Fill( tubeAngularCoordinate, tubePosition[2], var);
			if(barrelHitMap)  barrelHitMap->SetPoint(pmtNumber, tubeAngularCoordinate, tubePosition[1], var);
		}
	}
	
	// we only want the 'z' (charge/time) to be represented by colour, we don't actually want
	// it in 3D space, but apparently this isn't possible so just orient the canvas appropriately.
	// This is the officially suggested solution but damn is it hacky :/
	/* Reference:
	// x-z plane
	canvas->SetTheta(0);
	canvas->SetPhi(0);
	// y-z plane
	canvas->SetTheta(0);
	canvas->SetPhi(-90);
	// x-y plane
	canvas->SetTheta(90);
	canvas->SetPhi(0);
	*/
	
	// select caps (top half)
	displayCanvas->cd(1);
	// draw top cap (LHS)
	displayPad->cd(1);
	if(plotStyle==0) topCapHitMap->Draw("PCOL");
	if(plotStyle==1) topCapHeatMap->Draw("COL");
	gPad->GetView()->TopView();
	// draw bottom cap (RHS)
	displayPad->cd(2);
	if(plotStyle==0) bottomCapHitMap->Draw("PCOL");
	if(plotStyle==1) bottomCapHeatMap->Draw("COL");
	gPad->GetView()->TopView();
	// select barrel (bottom half)
	displayCanvas->cd(2);
	// draw barrel
	if(plotStyle==0) barrelHitMap->Draw("PCOL");
	if(plotStyle==1) barrelHeatMap->Draw("COL");
	gPad->GetView()->TopView();
	
	gPad->WaitPrimitive();
	
	//displayCanvas->SaveAs("HeatMap.png");
	
	return true;
}


bool evDisp::Finalise(){
	
	Log(toolName+" performing cleanup",v_debug,verbosity);
	// graphs
	if(topCapHitMap) delete topCapHitMap;
	if(bottomCapHitMap) delete bottomCapHitMap;
	if(barrelHitMap) delete barrelHitMap;
	// histograms
	if(topCapHeatMap) delete topCapHeatMap;
	if(bottomCapHeatMap) delete bottomCapHeatMap;
	if(barrelHeatMap) delete barrelHeatMap;
	// canvas
	delete displayCanvas;
	// PMT tube properties table
	delete myConnectionTable;
	
	return true;
}

bool evDisp::GetData(){
	myTreeReader->Get("TQREAL", myTQReal);
	return true;
}
