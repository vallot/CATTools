// HOW TO COMPILE
//
// g++  -L../src -l Toto -I `root-config --incdir` `root-config --libs` SanityCheck.cc -o SanityCheck
//

#include <iomanip>

#include "../interface/CatGenEvent.h"
#include "../interface/CatEvent.h"
#include "../interface/CatRun.h"

#include "../interface/CatLepton.h"
#include "../interface/CatMuon.h"
#include "../interface/CatElectron.h"

#include "../interface/CatPhoton.h"
#include "../interface/CatJet.h"
#include "../interface/CatCaloJet.h"
#include "../interface/CatJPTJet.h"
#include "../interface/CatPFJet.h"

#include "../interface/CatMET.h"
#include "../interface/CatCaloMET.h"
#include "../interface/CatPFMET.h"

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
#include <TChain.h>
#include <TClonesArray.h>

#include <fstream>
#include <sstream>

#include <sys/stat.h>

using namespace cat;
using namespace std;

void Tokenize(const string& str,
                      vector<string>& tokens,
                      const string& delimiters = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}


int main(int argc, char *argv[]){

  std::map<std::string, TH1F* > histos;

  std::map<std::string, std::pair<std::string,TClonesArray*> > arrays;

  std::vector<TClonesArray*> objects;
    
  int verbosity                 = 5;
  //0 muet
  //1 Main Info
  //2
  //3 
  //4 Info for each event
  //5 Debug

  bool onlyPU = false;

   // RETRIEVING LIST OF FILENAMES TO CHECK

  if (argc < 3) {

    cout << "Usage: ./SanityCheck --inputfiles file1;file2;fileN (--onlyPU)\n\n" << endl;

    exit(0);

  } else if (argc == 3 && !strstr(argv[1],"--inputfiles")) {

    cout << "Usage: ./SanityCheck --inputfiles file1;file2;fileN (--onlyPU)\n\n" << endl;

    exit(0);

  }

  for (unsigned int i=0;i<argc;i++)
    if (strstr(argv[i],"--onlyPU"))
      onlyPU=true;

  vector<string> fileNames;
  
  Tokenize(argv[2], fileNames, ";");

  if (fileNames.size() == 1) {

    if (!strstr(fileNames[0].c_str(),".root")) {

      cout << "A file containing a file-list was provided" << endl;

      fileNames.clear();

      string line;
      ifstream myfile (argv[2]);
      if (myfile.is_open())
	{
	  while ( myfile.good() )
	    {
	      getline (myfile,line);
	      if (strstr(line.c_str(),".root")) 
		fileNames.push_back((string)line);
	    }
	  myfile.close();
	}
      
      else cout << "Unable to open file"; 


    }

  }

  // CHECKING THE FILECONTENT FOR FILE 0 AND COUNT EVENTS FOR ALL FILES

  unsigned int nEvents = 0; 

  if (verbosity > 0) cout << "Declaring TChains" << endl;
  TChain* eventTree = new TChain("eventTree");
  TChain* runTree = new TChain("runTree");


  for (int fileID=0; fileID < fileNames.size(); fileID++) {
  
    if (verbosity > 0) cout << "Adding file " << fileNames.at(fileID) <<  endl;

    /*TFile* f = TFile::Open(fileNames.at(fileID).c_str());
    
    TTree* runTree = (TTree*) f->Get("runTree");
    TTree* eventTree = (TTree*) f->Get("eventTree");
    
    TBranch* run_br = (TBranch *) runTree->GetBranch("runInfos");
    CatRun* runInfos = 0;
    run_br->SetAddress(&runInfos);
    
    TBranch* event_br = (TBranch *) eventTree->GetBranch("Event");
    CatEvent* event = 0;
    event_br->SetAddress(&event);*/

    eventTree->Add(fileNames.at(fileID).c_str());
    runTree->Add(fileNames.at(fileID).c_str());
    
    nEvents += eventTree->GetEntriesFast();
      
  }

  CatRun* runInfos = new CatRun();
  runTree->SetBranchStatus("runInfos*",1);
  runTree->SetBranchAddress("runInfos",&runInfos);

  CatEvent* event = new CatEvent();
  eventTree->SetBranchStatus("Event*",1);
  eventTree->SetBranchAddress("Event",&event);

  if (verbosity > 0) cout << "Searching the TBranches" << endl;

  int nArrays = 0;
  TClonesArray* myArrays[2000];
  string myArrayClass[2000];
  string myArrayName[2000];

  for (int i=1; i<eventTree->GetListOfBranches()->GetEntriesFast(); i++) {
    
    TBranch * branch = (TBranch *)eventTree->GetListOfBranches()->At(i);
    
    TObject* obj = branch->GetListOfLeaves()->At(0);
    
    std::string ObjName = obj->GetName();

    if (onlyPU && !strstr(ObjName.c_str(),"PrimaryVertex")) continue;
    
    string::size_type position = ObjName.find_last_of("_");
    
    std::string className = "";
    
    if (strstr(ObjName.c_str(),"CaloJet"))
      className="cat::CatCaloJet";
    else if (strstr(ObjName.c_str(),"PFJet"))
      className="cat::CatPFJet";
    else if (strstr(ObjName.c_str(),"JPTJet"))
      className="cat::CatJPTJet";
    else if (strstr(ObjName.c_str(),"GenJet"))
      className="cat::CatGenJet";
    else if (strstr(ObjName.c_str(),"MCParticles"))
      className="cat::CatMCParticle";
    else if (strstr(ObjName.c_str(),"NPGenEvent"))
      className="cat::CatNPGenEvent";
    else if (strstr(ObjName.c_str(),"GenEvent"))
      className="cat::CatGenEvent";
    else if (strstr(ObjName.c_str(),"Muon"))
      className="cat::CatMuon";
    else if (strstr(ObjName.c_str(),"Electron"))
      className="cat::CatElectron";
    else if (strstr(ObjName.c_str(),"Photon"))
      className="cat::CatPhoton";
    else if (strstr(ObjName.c_str(),"TCMET"))
      className="cat::CatMET";
    else if (strstr(ObjName.c_str(),"CaloMET"))
	className="cat::CatCaloMET";
    else if (strstr(ObjName.c_str(),"PFMET"))
      className="cat::CatPFMET";
    else if (strstr(ObjName.c_str(),"MET"))
      className="cat::CatMET";
    else if (strstr(ObjName.c_str(),"MHT"))
      className="cat::CatMHT";
    else if (strstr(ObjName.c_str(),"PrimaryVertex"))
      className="cat::CatVertex";
    
    if (verbosity > 1) cout << "  Found Branch " << className << " " << ObjName.substr(0,position) << endl;
      
    arrays[ObjName.substr(0,position)]=std::pair<std::string,TClonesArray*>(className,new TClonesArray());
    
    char branchStatus[100];
    myArrays[nArrays] = new TClonesArray(className.c_str(), 0);
    myArrayClass[nArrays]=className;
    myArrayName[nArrays]=ObjName.substr(0,position);
    
    string branchName = ObjName.substr(0,position);
    sprintf(branchStatus,"%s*",branchName.c_str());
    
    eventTree->SetBranchStatus(branchStatus,1);
    eventTree->SetBranchAddress(branchName.c_str(),&myArrays[nArrays]);
    
    nArrays++;
    
  } 

  /*char branchStatus[100];
  objects.push_back(new TClonesArray("cat::CatMuon", 0));
  sprintf(branchStatus,"%s*","Muons_selectedPatMuons");
  eventTree->SetBranchStatus(branchStatus,1);
  eventTree->SetBranchAddress("Muons_selectedPatMuons",&objects[0]);*/

  //arrays["Muons_selectedPatMuons"]=std::pair<std::string,TClonesArray*>("cat::CatMuon",objects[0]);

  if (verbosity > 1) cout << "eventTree->GetEntriesFast(): " << eventTree->GetEntries() << endl;
  if (verbosity > 1) cout << "runTree->GetEntriesFast(): " << runTree->GetEntries() << endl;

  if (verbosity > 0) cout << "Looping over " << eventTree->GetEntries() << " events " << endl;
  for(unsigned int ievt=0; ievt<eventTree->GetEntries(); ievt++) {
    //for(unsigned int ievt=0; ievt<1; ievt++) {

    eventTree->GetEntry(ievt);
    runTree->GetEntry(0);
    
    if (verbosity > 0 && ievt % 1000 == 0) std::cout<<"Processing the "<<ievt<<"th event." << flush<<"\r";

    // PileUp plot -> not yet in these toptrees

    string hist="pileup";
    if (histos.find(hist) == histos.end()) histos[hist]=new TH1F((hist).c_str(),(hist+";nPu").c_str(),76,-0.5,75.5);    
    histos[hist]->Fill(event->nPu(0));
    
    hist="kt6PFJets_rho";
    if (histos.find(hist) == histos.end()) histos[hist]=new TH1F((hist).c_str(),hist.c_str(),80,0,40);
    histos[hist]->Fill(event->kt6PFJets_rho());    
    
    
    for (unsigned int p=0; p<nArrays; p++) {
      
      //cout << myArrayClass[p] << " " << myArrayName[p] << " " << myArrays[p]->GetEntries() << endl;

      
      //
      // define global plots we want to see
      //
      
      string hist="";
      
      hist="_size";
      if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),100,0,100);
      
      for (unsigned int o=0; o<myArrays[p]->GetEntries(); o++) {
	
	// because c++ it not flexible we need to write the elaborate code below.... Sorry:'(
	
	//
	//
	// Dear object developer:) You can find some sample code here to add plots so all your changes will happen below this comment.....
	//
	//
	
	//
	// HLT
	//
	
	/* NO IMPLEMENTATION YET */
	
	//
	// PRIMARY VERTEX
	//
	
	if (strstr(myArrayClass[p].c_str(),"cat::CatVertex")) {
	  
	  if(o==0) histos[myArrayName[p]+"_size"]->Fill(myArrays[p]->GetEntries());
	  
	}
	
	//
	// GenEvent
	//
	
	else if (strstr(myArrayClass[p].c_str(),"cat::CatGenEvent") || strstr(myArrayClass[p].c_str(),"cat::CatNPGenEvent")) {
	  
	  if(o==0) histos[myArrayName[p]+"_size"]->Fill(myArrays[p]->GetEntries());
	  
	}
	
	//
	// MCParticles
	//
	
	else if (strstr(myArrayClass[p].c_str(),"cat::CatMCParticle")) {
	  
	  if(o==0) histos[myArrayName[p]+"_size"]->Fill(myArrays[p]->GetEntries());
	  
	}
	
	//
	// JETS
	//
	
	// now putting all jet types together, these can be addressed seperately in the future if needed
	else if (strstr(myArrayClass[p].c_str(),"cat::CatGenJet") || strstr(myArrayClass[p].c_str(),"cat::CatCaloJet") || strstr(myArrayClass[p].c_str(),"cat::CatPFJet") || strstr(myArrayClass[p].c_str(),"cat::CatJPTJet") ){
	  
	  if(o==0) histos[myArrayName[p]+"_size"]->Fill(myArrays[p]->GetEntries());
	  
	  // some local plots
	  hist="_pt";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),250,0,500);
	  hist="_eta";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),300,-6,6);	  
	  hist="_phi";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),250,-5,5);
	  
	  histos[myArrayName[p]+"_pt"]->Fill(((CatJet*)myArrays[p]->At(o))->Pt());
	  histos[myArrayName[p]+"_eta"]->Fill(((CatJet*)myArrays[p]->At(o))->Eta());
	  histos[myArrayName[p]+"_phi"]->Fill(((CatJet*)myArrays[p]->At(o))->Phi());
	 
	  // b-tagging

	  hist="_bTagDisc_TCHE";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),50,-10,30);
	  histos[myArrayName[p]+hist]->Fill(((CatJet*)myArrays[p]->At(o))->btag_trackCountingHighEffBJetTags());

	  hist="_bTagDisc_TCHP";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),50,-10,30);
	  histos[myArrayName[p]+hist]->Fill(((CatJet*)myArrays[p]->At(o))->btag_trackCountingHighPurBJetTags());
	  
	  hist="_bTagDisc_SSVHE";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),50,0,8);
	  histos[myArrayName[p]+hist]->Fill(((CatJet*)myArrays[p]->At(o))->btag_simpleSecondaryVertexHighEffBJetTags());

	  hist="_bTagDisc_SSVHP";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),50,0,8);
	  histos[myArrayName[p]+hist]->Fill(((CatJet*)myArrays[p]->At(o))->btag_simpleSecondaryVertexHighPurBJetTags());

	  hist="_bTagDisc_JP";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),50,0,3);
	  histos[myArrayName[p]+hist]->Fill(((CatJet*)myArrays[p]->At(o))->btag_jetProbabilityBJetTags());

	  hist="_bTagDisc_JBP";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),50,0,8);
	  histos[myArrayName[p]+hist]->Fill(((CatJet*)myArrays[p]->At(o))->btag_jetBProbabilityBJetTags());

	  hist="_bTagDisc_CSV";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),50,-1,2);
	  histos[myArrayName[p]+hist]->Fill(((CatJet*)myArrays[p]->At(o))->btag_combinedSecondaryVertexBJetTags());

	  hist="_bTagDisc_CSVRetrained";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),50,-1,2);
	  histos[myArrayName[p]+hist]->Fill(((CatJet*)myArrays[p]->At(o))->btag_combinedSecondaryVertexRetrainedBJetTags());

	  hist="_bTagDisc_softPFMuon";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),50,-1,2);
	  histos[myArrayName[p]+hist]->Fill(((CatJet*)myArrays[p]->At(o))->btag_softPFMuonRetrainedBJetsTags());

	  hist="_bTagDisc_softPFElectron";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),50,-1,2);
	  histos[myArrayName[p]+hist]->Fill(((CatJet*)myArrays[p]->At(o))->btag_softPFElectronRetrainedBJetsTags());

	  hist="_bTagDisc_Combined_CSV_SL";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),50,-1,2);
	  histos[myArrayName[p]+hist]->Fill(((CatJet*)myArrays[p]->At(o))->btag_combinedCSVSLBJetTags());
	  
	  hist="_bTagDisc_Combined_CSV_JP";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),50,-1,2);
	  histos[myArrayName[p]+hist]->Fill(((CatJet*)myArrays[p]->At(o))->btag_combinedCSVJPBJetTags());

	   hist="_bTagDisc_Combined_CSV_JP_SL";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),50,-1,2);
	  histos[myArrayName[p]+hist]->Fill(((CatJet*)myArrays[p]->At(o))->btag_combinedCSVJPSLBJetTags());

	}

	// 
	// MET
	//
	
	// now putting all met types together, these can be addressed seperately in the future if needed
	else if (strstr(myArrayClass[p].c_str(),"cat::CatCaloMET") || strstr(myArrayClass[p].c_str(),"cat::CatPFMET") || strstr(myArrayClass[p].c_str(),"cat::CatMET") ){
	  
	  if(o==0) histos[myArrayName[p]+"_size"]->Fill(myArrays[p]->GetEntries());
	  
	  // some local plots
	  hist="_Et";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),250,0,500);
	  
	  histos[myArrayName[p]+"_Et"]->Fill(((CatMET*)myArrays[p]->At(o))->Et());
	  
	  
	}
	
	//
	// MUON
	//
	
	else if (strstr(myArrayClass[p].c_str(),"cat::CatMuon")){
	  
	  if(o==0) histos[myArrayName[p]+"_size"]->Fill(myArrays[p]->GetEntries());
	  
	  // some local plots
	  hist="_pt";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),250,0,500);
	  hist="_eta";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),300,-6,6);	  
	  hist="_phi";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),250,-5,5);
	  hist="_chargedHadronIso";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),250,-5,25);
	  hist="_photonIso";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),250,-5,25);
	  hist="_neutralHadronIso";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),250,-5,25);
	  hist="_nofValidPixelHits";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),15,0,15);
	  hist="_nofMatchedStations";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),15,0,15);
	  
	  histos[myArrayName[p]+"_pt"]->Fill(((CatMuon*)myArrays[p]->At(o))->Pt());
	  histos[myArrayName[p]+"_eta"]->Fill(((CatMuon*)myArrays[p]->At(o))->Eta());
	  histos[myArrayName[p]+"_phi"]->Fill(((CatMuon*)myArrays[p]->At(o))->Phi());
	  histos[myArrayName[p]+"_chargedHadronIso"]->Fill(((CatMuon*)myArrays[p]->At(o))->chargedHadronIso());
	  histos[myArrayName[p]+"_photonIso"]->Fill(((CatMuon*)myArrays[p]->At(o))->photonIso());
	  histos[myArrayName[p]+"_neutralHadronIso"]->Fill(((CatMuon*)myArrays[p]->At(o))->neutralHadronIso());
	  histos[myArrayName[p]+"_nofValidPixelHits"]->Fill(((CatMuon*)myArrays[p]->At(o))->nofValidPixelHits());
	  histos[myArrayName[p]+"_nofMatchedStations"]->Fill(((CatMuon*)myArrays[p]->At(o))->nofMatchedStations());
	  
	}
	
	//
	// ELECTRON
	//
	
	else if (strstr(myArrayClass[p].c_str(),"cat::CatElectron")){
	  
	  if(o==0) histos[myArrayName[p]+"_size"]->Fill(myArrays[p]->GetEntries());
	  
	  // some local plots
	  hist="_pt";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),250,0,500);
	  hist="_eta";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),300,-6,6);	  
	  hist="_phi";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),250,-5,5);
	  hist="_chargedHadronIso";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),250,-5,25);
	  hist="_photonIso";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),250,-5,25);
	  hist="_neutralHadronIso";
	  if (histos.find(myArrayName[p]+hist) == histos.end()) histos[myArrayName[p]+hist]=new TH1F((myArrayName[p]+hist).c_str(),(myArrayName[p]).c_str(),250,-5,25);	  
	  
	  histos[myArrayName[p]+"_pt"]->Fill(((CatElectron*)myArrays[p]->At(o))->Pt());
	  histos[myArrayName[p]+"_eta"]->Fill(((CatElectron*)myArrays[p]->At(o))->Eta());
	  histos[myArrayName[p]+"_phi"]->Fill(((CatElectron*)myArrays[p]->At(o))->Phi());
	  histos[myArrayName[p]+"_chargedHadronIso"]->Fill(((CatElectron*)myArrays[p]->At(o))->chargedHadronIso());
	  histos[myArrayName[p]+"_photonIso"]->Fill(((CatElectron*)myArrays[p]->At(o))->photonIso());
	  histos[myArrayName[p]+"_neutralHadronIso"]->Fill(((CatElectron*)myArrays[p]->At(o))->neutralHadronIso());
	  
	}
	
      }
      
      //delete myArray;
    }
    
  } cout << endl; // end of event loop

  if(verbosity>1) cout <<"Writing the plots" << endl;

  mkdir("Output",0755);

  std::map<std::string, vector<TH1F*> > new_histos;

  for (std::map<std::string, TH1F* >::const_iterator it=histos.begin(); it != histos.end(); ++it) {

    new_histos[(string)it->second->GetTitle()].push_back(it->second);

  }

  vector<TCanvas*> canvas; // yes another vector but this is the last one, promise!
  vector<TFile*> out; // yes another vector but this is the last one, promise!

  for (std::map<std::string, vector<TH1F*> >::const_iterator it=new_histos.begin(); it != new_histos.end(); ++it) {

    //if (verbosity > 3) cout << myArrayName[p] << " " <<  it->second.size() << endl;

    if (it->second.size() == 0)
      continue;

    int dest = 1;

    int numCanvas = 0;
    
    for (unsigned int j=0; j<it->second.size(); j++) {

      if (j ==0 || dest % 5 == 0) {

	stringstream s; s << canvas.size()+1;

	string cName = "canvas_"+s.str();

	string cTitle = ("Output/"+(string)it->second[j]->GetTitle()).c_str(); // this to know in which dir we need to store it

	canvas.push_back(new TCanvas(cName.c_str(),cTitle.c_str(),800,600));

	mkdir(canvas[canvas.size()-1]->GetTitle(),0755);

	string rootFileName= (string)canvas[canvas.size()-1]->GetTitle()+"/"+(string)canvas[canvas.size()-1]->GetName()+".root";

	out.push_back(new TFile(rootFileName.c_str(),"RECREATE"));
  
	numCanvas = canvas.size()-1;

	if (it->second.size() == 1)
	  canvas[numCanvas]->Divide(1,1);
	else
	  canvas[numCanvas]->Divide(2,2);


	dest=1;

      }

      canvas[numCanvas]->cd(dest);

      gPad->SetGrid();

      if (it->second[j]->Integral() == 0)

	gPad->SetFillColor(kRed);
      
      it->second[j]->SetTitle(it->second[j]->GetName()); // put back title

      it->second[j]->Draw();

      out[out.size()-1]->cd();

      it->second[j]->Write();

      dest++;

    }

  }
  
  for (unsigned int c=0; c<canvas.size(); c++) {

    string saveTitle = (string)canvas[c]->GetTitle()+"/"+(string)canvas[c]->GetName()+".png";

    canvas[c]->SaveAs(saveTitle.c_str());

    out[c]->cd();

    canvas[c]->Write();

    out[c]->Close();

  }

  if(verbosity>1) cout<<"End of the Macro"<<endl;

  return 0;
}
