#include <iomanip>
#include "../interface/CatLepton.h"
#include "../interface/CatMuon.h"
#include "../interface/CatElectron.h"
#include "../interface/CatPhoton.h"
#include "../interface/CatJet.h"
#include "../interface/CatCaloJet.h"
#include "../interface/CatPFJet.h"
#include "../interface/CatMET.h"
#include "../interface/CatTrackMET.h"
#include "../interface/CatGenEvent.h"
#include "../interface/CatEvent.h"
#include "../interface/CatRun.h"
#include "../interface/CatParticle.h"
#include "../interface/CatMCParticle.h"
#include "../interface/CatVertex.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TBranch.h>
#include <TTree.h>
#include <string>
#include <algorithm>
#include <vector>


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


int main(int argc, char** argv){

  // RETRIEVING LIST OF FILENAMES TO CHECK

  if (argc != 3) {

    cout << "Usage: ./catContentDump --inputfiles file1;file2;fileN\n\n" << endl;

    exit(0);

  } else if (argc == 3 && !strstr(argv[1],"--inputfiles")) {

    cout << "Usage: ./catContentDump --inputfiles file1;file2;fileN\n\n" << endl;

    exit(0);

  }

  vector<string> fileNames;
  
  Tokenize(argv[2], fileNames, ";");


  // CHECKING THE FILECONTENT FOR FILE 0 AND COUNT EVENTS FOR ALL FILES

  unsigned int nEvents = 0; 

  for (int fileID=0; fileID < fileNames.size(); fileID++) {
  
    //cout << fileNames.at(fileID) << endl;

    TFile* f = TFile::Open(fileNames.at(fileID).c_str());
    
    TTree* runTree = (TTree*) f->Get("runTree");
    TTree* eventTree = (TTree*) f->Get("eventTree");
    
    TBranch* run_br = (TBranch *) runTree->GetBranch("runInfos");
    CatRun* runInfos = 0;
    run_br->SetAddress(&runInfos);
    
    TBranch* event_br = (TBranch *) eventTree->GetBranch("Event");
    CatEvent* event = 0;
    event_br->SetAddress(&event);
    
    nEvents += eventTree->GetEntriesFast();

    if (fileID == 0) {

      cout << "* Dumping the event content of the cat" << endl;
      
      for (int i=1; i<eventTree->GetListOfBranches()->GetEntriesFast(); i++) {
	
	TBranch * branch = (TBranch *)eventTree->GetListOfBranches()->At(i);
	
	TObject* obj = branch->GetListOfLeaves()->At(0);
	
	std::string ObjName = obj->GetName();
	
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
	else if (strstr(ObjName.c_str(), "TrackMET"))
	  className="cat::CatTrackMET";  
	else if (strstr(ObjName.c_str(),"MHT"))
	  className="cat::CatMHT";
	else if (strstr(ObjName.c_str(),"PrimaryVertex"))
	  className="cat::CatVertex";
	

	cout << "- " << className << setw(5) << " -> " << "\"" << ObjName.substr(0,position) << "\""  << endl;
   
	    }

     //runinfos

     runTree->GetEvent(0);

     if (runInfos->hltInputTag() != "")
       cout << "- " << "cat::CatRun" << setw(5) << " -> " << "\"" << runInfos->hltInputTag() << "\""  << endl;

    }

  }

  //cout << "\n* The cat file contains " << nEvents << " events\n" << endl;


}
  
