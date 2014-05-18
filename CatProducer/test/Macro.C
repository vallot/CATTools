#include <iomanip>
#include "../interface/CatMuon.h"
#include "../interface/CatElectron.h"
#include "../interface/CatJe t.h"
#include "../interface/CatCaloJet.h"
#include "../interface/CatPFJet.h"
#include "../interface/CatMET.h"
#include "../interface/CatGenEvent.h"
#include "../interface/CatEvent.h"
#include "../interface/CatRun.h"
#include "../interface/CatParticle.h"
#include "../interface/CatMCParticle.h"
#include "../interface/CatVertex.h"
#include "../interface/CatHLTInfo.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TBranch.h>
#include <TTree.h>

#include <fstream>
#include <sstream>

using namespace cat;

int main(){
        int verbosity                 = 1;
	//0 muet
	//1 Main Info
	//2
	//3 
	//4 Info for each event
	//5 Debug

	bool realData						= false;
	bool doHLT                    = true;
	bool doMC                     = false;
	bool doCaloJet                = true;
	bool doMuon                   = false;
	bool doElectron               = false;
	bool doMET                    = false;
	bool doGenEvent               = false;

	cout << "Opening file" << endl;
	TFile* f = new TFile("21052010_123157_TOPTREE.root");
	//TFile* f=new TFile("cat_pythia.root");
	//TFile* f=new TFile("/user/jmmaes/CMSSW/CMSSW_2_2_4_Sanity/CMSSW_2_2_4/src/CATTools/CatProducer/test/cat.root");
	TTree* runTree = (TTree*) f->Get("runTree");
	TTree* eventTree = (TTree*) f->Get("eventTree");

	TBranch* run_br = (TBranch *) runTree->GetBranch("runInfos");
	CatRun* runInfos = 0;
	run_br->SetAddress(&runInfos);
	
	TBranch* event_br = (TBranch *) eventTree->GetBranch("Event");
	CatEvent* event = 0;
	event_br->SetAddress(&event);
	
	if(verbosity > 0) cout<<"Declaring Branches and TClonesArray"<<endl;
	//Declartion of Branches and TClonesArray
	TBranch* mcParticles_br;
	TClonesArray* mcParticles;
	TBranch* caloJets_br;
	TClonesArray* caloJets;
	TBranch* muons_br;
	TClonesArray* muons;
	TBranch* electrons_br;
	TClonesArray* electrons;
	TBranch* genEvents_br;
	TClonesArray* genEvents;
	TBranch* mets_br;
	TClonesArray* mets;

	if(doMC)
	{
		mcParticles_br = (TBranch *) eventTree->GetBranch("MCParticles");
		mcParticles = new TClonesArray("cat::CatParticle", 0);
		mcParticles_br->SetAddress(&mcParticles);
	}

	if(doCaloJet)
	{
		caloJets_br = (TBranch *) eventTree->GetBranch("CaloJets_selectedPatJetsAK5Calo");
//		caloJets_br = (TBranch *) eventTree->GetBranch("Jets_iterativeCone5CaloJets");
		caloJets = new TClonesArray("cat::CatCaloJet", 0);
		caloJets_br->SetAddress(&caloJets);
	}

	if(doMuon)
	{
		muons_br = (TBranch *) eventTree->GetBranch("Muons");
		muons = new TClonesArray("cat::CatMuon", 0);
		muons_br->SetAddress(&muons);
	}
		
	if(doElectron)
	{
		electrons_br = (TBranch *) eventTree->GetBranch("Electrons");
		electrons = new TClonesArray("cat::CatElectron", 0);
		electrons_br->SetAddress(&electrons);
	}
		
        if(doMET)
	{
		mets_br = (TBranch *) eventTree->GetBranch("MET");
		mets = new TClonesArray("cat::CatMET", 0);
		mets_br->SetAddress(&mets);
	}
        
	if(doGenEvent)
	{
		genEvents_br = (TBranch *) eventTree->GetBranch("GenEvent");
		genEvents = new TClonesArray("cat::CatGenEvent", 0);
		genEvents_br->SetAddress(&genEvents);
	}
		
		
	// HLT Run Summary
	if (doHLT)
	{
	  cout << "runTree->GetEntries()="<<runTree->GetEntries()<<endl;
	  runTree->GetEvent(0);
	  cout << dec << endl;
		cout << "HLT-Report " << "---------- Event  Summary ------------\n";
		cout << "HLT-Report"
				<< " Events total = " << runInfos->nHLTEvents()
				<< "  wasrun = " << runInfos->nHLTWasRun()
				<< "  passed = " << runInfos->nHLTAccept()
				<< "  errors = " << runInfos->nHLTErrors()
				<< "\n";

		cout << endl;
	
	}

   //Declaration of histograms
   TH1F h_PtJets("PtCaloJets","Pt of caloJets",50,0,500);
	TH1F h_distrib("distrib","",100,-1,1);
	//
		vector< vector<int> > runLumiInfo;

	if(realData)
	{
		if(verbosity > 3) cout << "Reading in JSON file " << endl;

		string inputJSON;

		ifstream myfile ("JSON_test.txt");
		if (myfile.is_open())
		{
			getline (myfile,inputJSON); // Only the first line is needed
			myfile.close();
		}
	
		vector<string> splittedInputJSON;
		size_t begin = 2, end = 2;

		while(end < inputJSON.size())
		{
			end = inputJSON.find("]], \"",begin);
			string splitted = inputJSON.substr(begin, end - begin + 1);
			begin = end + 5;
		
			size_t tempEnd = splitted.find("\": [[", 0);
			string runNr = splitted.substr(0, tempEnd);
			stringstream ss(runNr);
			int runNumber = 0;
			ss >> runNumber;
		
			string remain = splitted.substr(tempEnd + 4, splitted.size() - ( tempEnd + 3 ) );
			size_t tempEnd2 = remain.find("]", 0);
			size_t tempBegin2 = 0;

			while(tempEnd2 < remain.size())
			{
				string lumiInfo = remain.substr(tempBegin2 + 1, tempEnd2 - tempBegin2 - 1);
				tempBegin2 = tempEnd2 + 3;
				tempEnd2 = remain.find("]", tempBegin2);
			
				// parse lumiInfo string
				size_t tempBegin3 = lumiInfo.find(", ",0);
				string minLS = lumiInfo.substr(0,tempBegin3);
				string maxLS = lumiInfo.substr(tempBegin3 + 2, lumiInfo.size());
				int minLumiSection = 0;
				int maxLumiSection = 0;
				stringstream ssMin(minLS);		
				stringstream ssMax(maxLS);
				ssMin >> minLumiSection;
				ssMax >> maxLumiSection;
		
				vector<int> tempInfo;
				tempInfo.push_back(runNumber);
				tempInfo.push_back(minLumiSection);
				tempInfo.push_back(maxLumiSection);
				runLumiInfo.push_back(tempInfo);
			}
		}
	}

	unsigned int nEvents = (int)eventTree->GetEntries();

	//nEvents = 5000;

	int previous_run = 0;

	for(unsigned int ievt=0; ievt<nEvents; ievt++)
	{
		eventTree->GetEvent(ievt);
		runTree->GetEvent(0);
		if(verbosity>3) cout <<"event "<< ievt <<endl;
		if(verbosity>3) cout<<"event->nb()="<<event->nb()<<endl;

		bool goodEvent = true;
		if(realData)
		{
			goodEvent = false;
			for(unsigned int k=0; k<runLumiInfo.size(); k++)
			{
				if(event->runId() == runLumiInfo[k][0] && event->lumiBlockId() >= runLumiInfo[k][1] && event->lumiBlockId() <= runLumiInfo[k][2])
					goodEvent = true;
			}
		}
		
		if (doHLT)
		{

		  cat::CatHLTInfo hltInfo = runInfos->getHLTinfo(event->runId()); // get the triggerconfiguration for this run
		  
		  // get the trigger bits
		  
		  int HLT_Muon_bit = hltInfo.hltPath("HLT_Mu9");
		  
		  if(verbosity>5) cout << "HLT_Muon_bit: " <<  HLT_Muon_bit << endl;

		  if( !event->trigHLT(HLT_Muon_bit) ) // skip event if it doesn't pass our selected muon trigger

		    goodEvent=false;

		  // print hlt summary for every new run 

		  if ( event->runId() != previous_run ) {

		    previous_run = event->runId();

		    stringstream runID; runID << event->runId();
		  
		    cout << "HLT-Report " << "---------- HLTrig Summary for RUN "+runID.str()+" ------------\n";
		    cout << "HLT-Report   HLT Bit#     WasRun     Passed     Errors  Name\n";

		    for (unsigned int i=0; i!=hltInfo.nHLTPaths(); ++i)
		      {

			printf("HLT-Report %10u %10u %10u %10u  ", i, hltInfo.hltWasRun(i), hltInfo.hltAccept(i), hltInfo.hltErrors(i));
			cout << hltInfo.hltNames(i) << endl;
		      }
		    
		    cout << endl;
		    //for(unsigned int ipath=36; ipath<40; ipath++) cout << "   " << runInfos->hltNames(ipath) << " decision=" << event->trigHLT(ipath) <<endl;
		  }

		}

		if(goodEvent == false) continue; // skip events that don't pass JSON or HLT

		//access to muons
		if(doMuon){
	  	  if(verbosity>4) cout<<"Access to muons"<<endl;
                  if(verbosity>3) cout<<"Nof muons: "<<muons->GetEntriesFast()<<endl;
		  CatMuon* muon;
		  for(int i=0;i<muons->GetEntriesFast();i++){
		    muon = (CatMuon*) muons->At(i);
		  }
		}
		
		//access to electrons
		if(doElectron){
		  if(verbosity>4) cout<<"Access to electrons"<<endl;
                  if(verbosity>3) cout<<"Nof electrons: "<<electrons->GetEntriesFast()<<endl;
		  CatElectron* electron;
		  for(int i=0;i<electrons->GetEntriesFast();i++){
		    electron = (CatElectron*) electrons->At(i);
		  }
	        } 
		
		//access to caloJets
		if(doCaloJet){
		  if(verbosity>4) cout<<"Access to caloJets"<<endl;
        if(verbosity>3) cout<<"Nof caloJets: "<<caloJets->GetEntriesFast()<<endl;
		  CatCaloJet* caloJet;
		  for(int i=0;i<caloJets->GetEntriesFast();i++){
		    caloJet = (CatCaloJet*) caloJets->At(i);
		    h_PtJets.Fill(caloJet->Pt());
		    //h_distrib.Fill(caloJet->btag_combinedSecondaryVertexBJetTags());
		    h_distrib.Fill(caloJet->btag_trackCountingHighEffBJetTags());
		  }
		}
		
		//access to GenEvent
		if(doGenEvent){
		  if(verbosity>4) cout<<"Access to GenEvent"<<endl;
		  if(genEvents->GetEntriesFast()>0){
		    CatGenEvent* genEvent = (CatGenEvent*) genEvents->At(0);
		    if(genEvent->isSemiLeptonic())
		    h_PtJets.Fill(genEvent->hadronicDecayTop().Pt());
		    
		  }
		  else if(verbosity>0) cout<<" No access to GenEvent in this entry"<<endl;
		}
		
		//access to MET
		if(doMET){
		  if(verbosity>4) cout<<"Access to MET"<<endl;
		  if(mets->GetEntriesFast()>0){
		    CatMET* met = (CatMET*) mets->At(0);
		  }
		  else if(verbosity>0) cout<<" No access to MET in this entry"<<endl;
		}



	} // end of loop over evts

	

      if(verbosity>1) cout<<"Writting histograms in the root-file ... "<<endl;
      
      TFile* fout = new TFile("MacroTreeOutput.root","RECREATE");
      fout->cd();
      //Write histograms
      h_PtJets.Write();
      h_distrib.Write();
      fout->Close();

      if(verbosity>1) cout<<"End of the Macro"<<endl;

}
