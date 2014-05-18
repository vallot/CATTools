#include <memory>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "./tinyxml/tinyxml.h"

#include "../interface/CatLepton.h"
#include "../interface/CatMuon.h"
#include "../interface/CatElectron.h"
#include "../interface/CatPhoton.h"
#include "../interface/CatJet.h"
#include "../interface/CatPFJet.h"
#include "../interface/CatCaloJet.h"
#include "../interface/CatJPTJet.h"
#include "../interface/CatGenJet.h"
#include "../interface/CatMET.h"
#include "../interface/CatCaloMET.h"
#include "../interface/CatPFMET.h"
#include "../interface/CatTrackMET.h"
#include "../interface/CatGenEvent.h"
#include "../interface/CatNPGenEvent.h"
#include "../interface/CatSpinCorrGen.h"
#include "../interface/CatEvent.h"
#include "../interface/CatRun.h"
#include "../interface/CatParticle.h"
#include "../interface/CatMCParticle.h"
#include "../interface/CatVertex.h"

#include "TFile.h"
#include "TBranch.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TObjString.h"

using namespace cat;

struct keepObjects
{
	TBranch* inBranch;
	TClonesArray* inArray;
	TClonesArray* outArray;
	string name;
	string type;
	float minPt;
	float maxEta;
	bool skipObjects;
	int minNObjects;
};

struct options
{
	bool skimOnHLT;
	string HLTPath1;
	string HLTPath2;
	bool HLTApplyAnd;
	string TriggerMenu;
	bool useJSON;
	string JSONFile;
};

vector<keepObjects> parseObjects(TiXmlDocument doc, TTree* outEventTree)
{
	vector<keepObjects> toKeep;
	TiXmlHandle hdl (&doc);
	TiXmlNode *node = 0;
	TiXmlElement *elem = 0;
	bool nodeFound = false;
	node = hdl.Node ();
	for( node = node->FirstChild(); node; node = node->NextSibling() )
	{
		if(node->Value () == string ("keepbranches"))
		{
			nodeFound = true;
			elem = node->FirstChildElement ();
			if(!elem)
			{
				cerr << "The node doesn't exist" << endl;
				delete node;
				delete elem;
			}
			while(elem)
			{
				keepObjects tempObj;
				tempObj.name = (TString) elem->Attribute("name");
				tempObj.type = (TString) elem->Attribute("type");
				
				tempObj.outArray = new TClonesArray((tempObj.type).c_str(), 1000);
				outEventTree->Branch ((tempObj.name).c_str(), "TClonesArray", &(tempObj.outArray) );
				
				elem->QueryFloatAttribute ("minPt", &(tempObj.minPt));
				elem->QueryFloatAttribute ("maxEta", &(tempObj.maxEta));
				
				int skipObj = 0, minNobj = 0;
				elem->QueryIntAttribute ("minNObjects", &(minNobj));
				tempObj.minNObjects = minNobj;
				elem->QueryIntAttribute ("skipObjects", &(skipObj));
				if(skipObj == 0)
					tempObj.skipObjects = false;
				else if(skipObj == 1)
					tempObj.skipObjects = true;
				else
					cerr << "Wrong skipObjects : " << skipObj << " for " << tempObj.name << endl;
				
				cout << "The skim will keep " << tempObj.name << ", of type " << tempObj.type << endl;
				cout << "With minPt = " << tempObj.minPt << " and maxEta = " << tempObj.maxEta << endl;
				if(tempObj.skipObjects)
					cout << "With skipObjects = true " << endl;
				else
					cout << "With skipObjects = false " << endl;

				toKeep.push_back(tempObj);
				
				elem = elem->NextSiblingElement ();	// iteration 
			}
		}
	}

	if(!nodeFound)
	{
		cerr << "The node doesn't exist" << endl;
		delete node;
		delete elem;
	}
	
	return toKeep;
}

vector<TString> parseFileName(TiXmlDocument doc, string name)
{
	vector<TString> inFileName;
	TiXmlHandle hdl (&doc);
	TiXmlNode *node = 0;
	TiXmlElement *elem = 0;
	bool nodeFound = false;
	node = hdl.Node ();
	for( node = node->FirstChild(); node; node = node->NextSibling() )
	{
		if (node->Value () == string (name))
		{
			nodeFound = true;
			elem = node->FirstChildElement ();
			if(!elem)
			{
				cerr << "The node doesn't exist" << endl;
				delete node;
				delete elem;
				exit (3);
			}
			while (elem)
			{
				inFileName.push_back( (TString) elem->Attribute("file") );
				elem = elem->NextSiblingElement ();	// iteration 
			}
		}
	}

	if(!nodeFound)
	{
		cerr << "The node doesn't exist" << endl;
		delete node;
		delete elem;
		exit (2);
	}
	return inFileName;
}

options parseOptions(TiXmlDocument doc)
{
	options tempOptions;
	TiXmlHandle hdl (&doc);
	TiXmlNode *node = 0;
	TiXmlElement *elem = 0;
	bool nodeFound = false;
	node = hdl.Node ();
	for( node = node->FirstChild(); node; node = node->NextSibling() )
	{
		if (node->Value () == string("options") )
		{
			nodeFound = true;
			elem = node->FirstChildElement ();
			if(!elem)
			{
				cerr << "The node doesn't exist" << endl;
				delete node;
				delete elem;
				exit (3);
			}
			while (elem)
			{
				int skimOnHLT = 0;
				elem->QueryIntAttribute("skimOnHLT", &skimOnHLT);
				if(skimOnHLT == 1)
				{
					tempOptions.skimOnHLT = true;
					tempOptions.HLTPath1 = elem->Attribute("HLTPath1");
					tempOptions.HLTPath2 = elem->Attribute("HLTPath2");
					tempOptions.TriggerMenu = elem->Attribute("TriggerMenu");
					int HLTApplyAnd = 0;
					elem->QueryIntAttribute("HLTApplyAnd", &HLTApplyAnd);
					if(HLTApplyAnd == 1) tempOptions.HLTApplyAnd = true;
					else tempOptions.HLTApplyAnd = false;
				}
				else tempOptions.skimOnHLT = false;
				
				int useJSON = 0;
				elem->QueryIntAttribute("useJSON", &useJSON);
				if(useJSON == 1)
				{
					tempOptions.useJSON = true;
					tempOptions.JSONFile = elem->Attribute("JSONFile");
				}
				else tempOptions.useJSON = false;
				
				elem = elem->NextSiblingElement ();	// iteration 
			}
		}
	}

	if(!nodeFound)
	{
		cerr << "The node doesn't exist" << endl;
		delete node;
		delete elem;
		exit (2);
	}
	
	return tempOptions;
}

int main()
{
	clock_t start = clock();

	// set verbosity equal to 0 (silent), 1 or 2 (debug)
	unsigned int verbosity = 0;

	// xml file
	char xmlfile[]="skim.xml";
	
	TiXmlDocument doc (xmlfile);
	
	if(!doc.LoadFile())
	{
		cerr << "Error while loading the xml file: " << xmlfile <<endl;
		cerr << " error #" << doc.ErrorId () << " : " << doc.ErrorDesc () << endl;
      return 1;
	}
	
	vector<TString> inFileName = parseFileName(doc, "inputdatasets");
	
	vector<TString> outFileName = parseFileName(doc, "outputfilename");
	
	options optionsToUse = parseOptions(doc);
	cout << "options:" << endl;
	cout << "HLTPath1 = " << optionsToUse.HLTPath1 << endl;
	cout << "HLTPath2 = " << optionsToUse.HLTPath2 << endl;
	cout << "HLTApplyAnd = " << optionsToUse.HLTApplyAnd << endl;
	cout << "JSONFile = " << optionsToUse.JSONFile << endl;

	cout << "output file: " << outFileName[0] << endl;

	TFile* outFile = new TFile(outFileName[0], "Recreate");
	outFile->cd();

	CatRun* outRunInfos = 0;
	TTree* outRunTree = new TTree("runTree", "Global Run Infos");
	outRunTree->Branch("runInfos", "cat::CatRun", &outRunInfos);
	
	CatEvent* outRootEvent = 0;
	TTree* outEventTree = new TTree("eventTree", "Event Infos");
	outEventTree->Branch("Event", "cat::CatEvent", &outRootEvent);
	
	cout << "Parsing objectsToKeep from xml file..." << endl;
	
	//	parse input objects from xml
	vector<keepObjects> objectsToKeep = parseObjects(doc, outEventTree);
	
	if( verbosity > 0 ) cout << "objectsToKeep.size() = " << objectsToKeep.size() << endl;
	
	// Prepare userInfo to be added to the outEventTree
	string info = "This cat was skimmed with the following properties:\n";

	for(unsigned int j=0; j !=objectsToKeep.size(); j++)
	{
		stringstream tmpMinPt, tmpMaxEta, tmpMinNObjects;
		tmpMinPt << objectsToKeep[j].minPt;
		tmpMaxEta << objectsToKeep[j].maxEta;
		tmpMinNObjects << objectsToKeep[j].minNObjects;

		info += "name = " + objectsToKeep[j].name + "  type = " + objectsToKeep[j].type;
		info += "  minPt = " + tmpMinPt.str() + "  maxEta = " + tmpMaxEta.str();
		if(objectsToKeep[j].skipObjects) info += "  skipObjects = 1";
		else info += "  skipObjects = 0";
		info += "  minNObjects = " + tmpMinNObjects.str() + "\n";
	}
	
	// Add UserInfo to outEventTree
	TObjString* userInfo = new TObjString();
	userInfo->SetString(info.c_str());
	
	cout << "userInfo that will be added to the outEventTree:\n" << userInfo->GetString() << endl;
	
	outEventTree->GetUserInfo()->Add( userInfo );
	
	vector< vector<int> > runLumiInfo; // To store the info from the JSON file
	
	if(optionsToUse.useJSON)
	{
		if(verbosity > 3) cout << "Reading in JSON file " << endl;

		string inputJSON;

		ifstream myfile ( (optionsToUse.JSONFile).c_str() );
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

	unsigned int nOutEvents = 0, nInEvents = 0, NHLTAccept = 0; 
	vector<unsigned int> hltAccept;
	
	for(unsigned int nFile = 0; nFile < inFileName.size(); nFile++)
	{
		cout << "input file[" << nFile << "] = " << inFileName[nFile] << endl;

		TFile* inFile = TFile::Open(inFileName[nFile]);
	
		TTree* inRunTree = (TTree*) inFile->Get("runTree");

		TBranch* inRunBranch = (TBranch *) inRunTree->GetBranch("runInfos");
		CatRun* inRunInfos = 0;
		inRunBranch->SetAddress(&inRunInfos);

		TTree* inEventTree = (TTree*) inFile->Get("eventTree");
	
		TBranch* inEventBranch = (TBranch *) inEventTree->GetBranch("Event");
		CatEvent* inEvent = 0;
		inEventBranch->SetAddress(&inEvent);
	
		if( verbosity > 0 ) cout << "inRunTree->GetEntries() = " << inRunTree->GetEntries()<<endl;
		inRunTree->GetEvent(0);
		
		if( verbosity > 1 ) cout << "outRunInfos->nHLTEvents() = " << outRunInfos->nHLTEvents() << endl;

		for(unsigned int j=0; j !=objectsToKeep.size(); j++)
		{
			objectsToKeep[j].inBranch = (TBranch *) inEventTree->GetBranch((objectsToKeep[j].name).c_str());
			objectsToKeep[j].inArray = new TClonesArray((objectsToKeep[j].type).c_str(), 0);
			objectsToKeep[j].inBranch->SetAddress( &(objectsToKeep[j].inArray) );
		}
	
		unsigned int nTempEvents = (int) inEventTree->GetEntries();
		nInEvents += nTempEvents;

		vector<cat::CatHLTInfo> tmpRunInfos = outRunInfos->copyHLTinfos();

		//loop over events
		for(unsigned int ievt=0; ievt<nTempEvents; ievt++)
		{
			if( verbosity > 1 ) cout << ">>> Trying to get event " << ievt << endl;

			inEventTree->GetEvent(ievt);
			
			// updating the outRunTree info
			if(ievt == 0 && nFile == 0)
			{
				outRunInfos->setHLTInputTag(inRunInfos->hltInputTag());
			}
			else
			{
				if(inRunInfos->hltInputTag() != outRunInfos->hltInputTag())
				{
					cout << "Different HLT inputTags!" << endl;
					cout << "inRunInfos->hltInputTag() = " << inRunInfos->hltInputTag() << "  outRunInfos->hltInputTag() = " << outRunInfos->hltInputTag() << endl;
					exit (4);
				}
			}
			
			// updating HLT info

			// The HLT info is stored per run in the CatRun. For file >= 0, we just copy the elements from the hltInfos vector and it's done...
			if (outRunInfos->getHLTinfo(inEvent->runId()).RunID() == 0) // if this run is not yet in the vector, add it.
			  tmpRunInfos.push_back(inRunInfos->getHLTinfo(inEvent->runId()));

			outRunInfos->setHLTinfos(tmpRunInfos); // put the new vector in CatRun

			if( verbosity > 1 ) cout << "outRunInfos->copyHLTinfos().size(): " << outRunInfos->copyHLTinfos().size() << endl;

			if( verbosity > 1 ) cout << ">>> Analyzing event " << ievt << endl;
			else if((int) ievt/10000 == (double) ievt/10000)  cout << ">>> Analyzing event " << ievt << endl;
		
			bool keepEvent = true;

			// apply JSON
			if(optionsToUse.useJSON)
			{
				bool goodEvent = false;
				for(unsigned int k=0; k<runLumiInfo.size(); k++)
				{
					if(inEvent->runId() == runLumiInfo[k][0] && inEvent->lumiBlockId() >= runLumiInfo[k][1] && inEvent->lumiBlockId() <= runLumiInfo[k][2])
						goodEvent = true;
				}
				if(goodEvent == false) keepEvent = false;
			}
			
			if(!keepEvent) continue;
			
			// apply selection based on HLT trigger
			if(optionsToUse.skimOnHLT)
			{
				int HLT1Bit = -1, HLT2Bit = -1;
				if(optionsToUse.TriggerMenu == inRunInfos->hltInputTag())
				{
					cat::CatHLTInfo hltInfo = inRunInfos->getHLTinfo(inEvent->runId());
					HLT1Bit = hltInfo.hltPath((optionsToUse.HLTPath1).c_str());
					if(optionsToUse.HLTPath2 != "")
					{
						HLT2Bit = hltInfo.hltPath((optionsToUse.HLTPath2).c_str());
						bool passed1 = inEvent->trigHLT(HLT1Bit);
						bool passed2 = inEvent->trigHLT(HLT2Bit);
						if(keepEvent && optionsToUse.HLTApplyAnd)
							keepEvent = passed1 && passed2;
						else if(keepEvent)
							keepEvent = passed1 || passed2;
					}
					else
						keepEvent = inEvent->trigHLT(HLT1Bit);
				}
				else
				{
					cerr << "Unknown HLT InputTag: " << optionsToUse.TriggerMenu << endl;
					exit(5);
				}
			}
			
			outFile->cd();
		
			outRootEvent = inEvent;
			
			for(unsigned int j=0; j !=objectsToKeep.size(); j++)
			{
				if( verbosity > 0 )
				{
					cout << "objectsToKeep[" << j << "] : " <<endl;
					cout << "name = " << objectsToKeep[j].name << " type = " << objectsToKeep[j].type << endl;
					cout << "minPt = " << objectsToKeep[j].minPt << " maxEta = " << objectsToKeep[j].maxEta << endl;
					if(objectsToKeep[j].skipObjects) cout << "skipObjects = true" << endl;
					else cout << "skipObjects = false" << endl;
				}
			}
		
			if(!keepEvent) continue;
		
			for(unsigned int j=0; j !=objectsToKeep.size(); j++)
			{
				if(objectsToKeep[j].type == "cat::CatVertex")
				{
					CatVertex* vertex;
					int verticesKept=0;

					for(int i=0; i<(objectsToKeep[j].inArray)->GetEntriesFast(); i++)
					{
						vertex = (CatVertex*) (objectsToKeep[j].inArray)->At(i);
						verticesKept++;					
						
						if( ! objectsToKeep[j].skipObjects )
							new( (*(objectsToKeep[j].outArray))[i] ) CatVertex(*vertex);
					}

					if(verticesKept < objectsToKeep[j].minNObjects)
					{
						keepEvent = false;
						if( verbosity > 1 ) cout << "Too small number of selected Primary Vertices: verticesKept = " << verticesKept << endl;
					}

					if( verbosity > 1 ) cout << "Processed " << objectsToKeep[j].name << endl;
					if( verbosity > 1 ) cout << "input = " << (objectsToKeep[j].inArray)->GetEntriesFast() << " output = " << (objectsToKeep[j].outArray)->GetEntriesFast() << endl;
				}

				else if(objectsToKeep[j].type == "cat::CatGenEvent")
				{
					CatGenEvent* genEvt;
					genEvt = (CatGenEvent*) (objectsToKeep[j].inArray)->At(0);
				
					if( ! objectsToKeep[j].skipObjects )
						new( (*(objectsToKeep[j].outArray))[0] ) CatGenEvent(*genEvt);

					if( verbosity > 1 ) cout << "Processed " << objectsToKeep[j].name << endl;
					if( verbosity > 1 ) cout << "input = " << (objectsToKeep[j].inArray)->GetEntriesFast() << " output = " << (objectsToKeep[j].outArray)->GetEntriesFast() << endl;
				}

				else if(objectsToKeep[j].type == "cat::CatNPGenEvent")
				{
					CatNPGenEvent* npGenEvt;
					npGenEvt = (CatNPGenEvent*) (objectsToKeep[j].inArray)->At(0);
				
					if( ! objectsToKeep[j].skipObjects )
						new( (*(objectsToKeep[j].outArray))[0] ) CatNPGenEvent(*npGenEvt);

					if( verbosity > 1 ) cout << "Processed " << objectsToKeep[j].name << endl;
					if( verbosity > 1 ) cout << "input = " << (objectsToKeep[j].inArray)->GetEntriesFast() << " output = " << (objectsToKeep[j].outArray)->GetEntriesFast() << endl;
				}

				else if(objectsToKeep[j].type == "cat::CatGenJet")
				{
					CatGenJet* genJet;
					int genJetsKept=0;
					for(int i=0; i<(objectsToKeep[j].inArray)->GetEntriesFast(); i++)
					{
						genJet = (CatGenJet*) (objectsToKeep[j].inArray)->At(i);
						bool keepGenJet = true;
					
						if(genJet->Pt() < objectsToKeep[j].minPt || fabs(genJet->Eta()) > objectsToKeep[j].maxEta)
						{
							keepGenJet = false;
							if( verbosity > 1 ) cout << "skip genJet with pT = " << genJet->Pt() << " and eta = " << genJet->Eta() << endl;
						}
						else genJetsKept++;
						
						if( ! objectsToKeep[j].skipObjects )
							new( (*(objectsToKeep[j].outArray))[i] ) CatGenJet(*genJet);
						else if(keepGenJet)
							new( (*(objectsToKeep[j].outArray))[genJetsKept-1] ) CatGenJet(*genJet);
					}
				
					if(genJetsKept < objectsToKeep[j].minNObjects)
					{
						keepEvent = false;
						if( verbosity > 1 ) cout << "Too small number of selected jets: genJetsKept = " << genJetsKept << endl;
					}
					
					if( verbosity > 1 ) cout << "Processed " << objectsToKeep[j].name << endl;
					if( verbosity > 1 ) cout << "input = " << (objectsToKeep[j].inArray)->GetEntriesFast() << " output = " << (objectsToKeep[j].outArray)->GetEntriesFast() << endl;
				}

				else if(objectsToKeep[j].type == "cat::CatPFJet")
				{
					CatPFJet* pfJet;
					int pfJetsKept=0;
					for(int i=0; i<(objectsToKeep[j].inArray)->GetEntriesFast(); i++)
					{
						pfJet = (CatPFJet*) (objectsToKeep[j].inArray)->At(i);
						bool keepPFJet = true;
						
						if(pfJet->Pt() < objectsToKeep[j].minPt || fabs(pfJet->Eta()) > objectsToKeep[j].maxEta)
						{
							keepPFJet = false;
							if( verbosity > 1 ) cout << "skip PFJet with pT = " << pfJet->Pt() << " and eta = " << pfJet->Eta() << endl;
						}
						else pfJetsKept++;
						
						if( ! objectsToKeep[j].skipObjects )
							new( (*(objectsToKeep[j].outArray))[i] ) CatPFJet(*pfJet);
						else if(keepPFJet)
							new( (*(objectsToKeep[j].outArray))[pfJetsKept-1] ) CatPFJet(*pfJet);			
					}
				
					if(pfJetsKept < objectsToKeep[j].minNObjects)
					{
						keepEvent = false;
						if( verbosity > 1 ) cout << "Too small number of selected jets: pfJetsKept = " << pfJetsKept << endl;
					}
				
					if( verbosity > 1 ) cout << "Processed " << objectsToKeep[j].name << endl;
					if( verbosity > 1 ) cout << "input = " << (objectsToKeep[j].inArray)->GetEntriesFast() << " output = " << (objectsToKeep[j].outArray)->GetEntriesFast() << endl;
				}

				else if(objectsToKeep[j].type == "cat::CatCaloJet")
				{
					CatCaloJet* caloJet;
					int caloJetsKept=0;
					for(int i=0; i<(objectsToKeep[j].inArray)->GetEntriesFast(); i++)
					{
						caloJet = (CatCaloJet*) (objectsToKeep[j].inArray)->At(i);
						bool keepCaloJet = true;
					
						if(caloJet->Pt() < objectsToKeep[j].minPt || fabs(caloJet->Eta()) > objectsToKeep[j].maxEta)
						{
							keepCaloJet = false;
							if( verbosity > 1 ) cout << "skip CaloJet with pT = " << caloJet->Pt() << " and eta = " << caloJet->Eta() << endl;
						}
						else caloJetsKept++;

						if( ! objectsToKeep[j].skipObjects )
							new( (*(objectsToKeep[j].outArray))[i] ) CatCaloJet(*caloJet);
						else if(keepCaloJet)
							new( (*(objectsToKeep[j].outArray))[caloJetsKept-1] ) CatCaloJet(*caloJet);			
					}
				
					if(caloJetsKept < objectsToKeep[j].minNObjects)
					{
						keepEvent = false;
						if( verbosity > 1 ) cout << "Too small number of selected caloJets: caloJetsKept = " << caloJetsKept << endl;
					}
				
					if( verbosity > 1 ) cout << "Processed " << objectsToKeep[j].name << endl;
					if( verbosity > 1 ) cout << "input = " << (objectsToKeep[j].inArray)->GetEntriesFast() << " output = " << (objectsToKeep[j].outArray)->GetEntriesFast() << endl;
				}

				else if(objectsToKeep[j].type == "cat::CatJPTJet")
				{
					CatJPTJet* JPTJet;
					int JPTJetsKept=0;
					for(int i=0; i<(objectsToKeep[j].inArray)->GetEntriesFast(); i++)
					{
						JPTJet = (CatJPTJet*) (objectsToKeep[j].inArray)->At(i);
						bool keepJPTJet = true;
					
						if(JPTJet->Pt() < objectsToKeep[j].minPt || fabs(JPTJet->Eta()) > objectsToKeep[j].maxEta)
						{
							keepJPTJet = false;
							if( verbosity > 1 ) cout << "skip JPTJet with pT = " << JPTJet->Pt() << " and eta = " << JPTJet->Eta() << endl;
						}
						else JPTJetsKept++;

						if( ! objectsToKeep[j].skipObjects )
							new( (*(objectsToKeep[j].outArray))[i] ) CatJPTJet(*JPTJet);			
						else if(keepJPTJet)
							new( (*(objectsToKeep[j].outArray))[JPTJetsKept-1] ) CatJPTJet(*JPTJet);			
					}
				
					if(JPTJetsKept < objectsToKeep[j].minNObjects)
					{
						keepEvent = false;
						if( verbosity > 1 ) cout << "Too small number of selected JPTJets: JPTJetsKept = " << JPTJetsKept << endl;
					}
				
					if( verbosity > 1 ) cout << "Processed " << objectsToKeep[j].name << endl;
					if( verbosity > 1 ) cout << "input = " << (objectsToKeep[j].inArray)->GetEntriesFast() << " output = " << (objectsToKeep[j].outArray)->GetEntriesFast() << endl;
				}
			
				else if(objectsToKeep[j].type == "cat::CatMCParticle")
				{
					CatMCParticle* mcparticle;
					int mcparticlesKept=0;
					for(int i=0; i<(objectsToKeep[j].inArray)->GetEntriesFast(); i++)
					{
						mcparticle = (CatMCParticle*) (objectsToKeep[j].inArray)->At(i);
						bool keepMCParticle = true;
	
						if(mcparticle->Pt() < objectsToKeep[j].minPt || fabs(mcparticle->Eta()) > objectsToKeep[j].maxEta)
						{
							keepMCParticle = false;
							if( verbosity > 1 ) cout << "skip MCparticle with pT = " << mcparticle->Pt() << " and eta = " << mcparticle->Eta() << endl;
						}
						else mcparticlesKept++;
					
						if( ! objectsToKeep[j].skipObjects )
							new( (*(objectsToKeep[j].outArray))[i] ) CatMCParticle(*mcparticle);
						else if(keepMCParticle)
							new( (*(objectsToKeep[j].outArray))[mcparticlesKept-1] ) CatMCParticle(*mcparticle);
					}

					if(mcparticlesKept < objectsToKeep[j].minNObjects)
					{
						keepEvent = false;
						if( verbosity > 1 ) cout << "Too small number of selected MCparticles: mcparticlesKept = " << mcparticlesKept << endl;
					}

					if( verbosity > 1 ) cout << "Processed " << objectsToKeep[j].name << endl;
					if( verbosity > 1 ) cout << "input = " << (objectsToKeep[j].inArray)->GetEntriesFast() << " output = " << (objectsToKeep[j].outArray)->GetEntriesFast() << endl;
				}
			
				else if(objectsToKeep[j].type == "cat::CatMET")
				{
					CatMET* met;
					int METKept=0;
					met = (CatMET*) (objectsToKeep[j].inArray)->At(0);
					bool keepMET = true;
				
					if(met->Pt() < objectsToKeep[j].minPt || fabs(met->Eta()) > objectsToKeep[j].maxEta)
					{
						keepMET = false;
						if( verbosity > 1 ) cout << "skip MET with pT = " << met->Pt() << " and eta = " << met->Eta() << endl;
					}
					else METKept++;
					
					if( ! objectsToKeep[j].skipObjects || keepMET )
						new( (*(objectsToKeep[j].outArray))[0] ) CatMET(*met);	

					if(METKept < objectsToKeep[j].minNObjects)
					{
						keepEvent = false;
						if( verbosity > 1 ) cout << "Too small number of selected MET: METKept = " << METKept << endl;
					}

					if( verbosity > 1 ) cout << "Processed " << objectsToKeep[j].name << endl;
					if( verbosity > 1 ) cout << "input = " << (objectsToKeep[j].inArray)->GetEntriesFast() << " output = " << (objectsToKeep[j].outArray)->GetEntriesFast() << endl;
				}

				else if(objectsToKeep[j].type == "cat::CatCaloMET")
				{
					CatCaloMET* met;
					int METKept=0;
					met = (CatCaloMET*) (objectsToKeep[j].inArray)->At(0);
					bool keepMET = true;
				
					if(met->Pt() < objectsToKeep[j].minPt || fabs(met->Eta()) > objectsToKeep[j].maxEta)
					{
						keepMET = false;
						if( verbosity > 1 ) cout << "skip CaloMET with pT = " << met->Pt() << " and eta = " << met->Eta() << endl;
					}
					else METKept++;
					
					if( ! objectsToKeep[j].skipObjects || keepMET )
						new( (*(objectsToKeep[j].outArray))[0] ) CatCaloMET(*met);

					if(METKept < objectsToKeep[j].minNObjects)
					{
						keepEvent = false;
						if( verbosity > 1 ) cout << "Too small number of selected CaloMET: METKept = " << METKept << endl;
					}

					if( verbosity > 1 ) cout << "Processed " << objectsToKeep[j].name << endl;
					if( verbosity > 1 ) cout << "input = " << (objectsToKeep[j].inArray)->GetEntriesFast() << " output = " << (objectsToKeep[j].outArray)->GetEntriesFast() << endl;
				}

				else if(objectsToKeep[j].type == "cat::CatPFMET")
				{
					CatPFMET* met;
					int METKept=0;
					met = (CatPFMET*) (objectsToKeep[j].inArray)->At(0);
					bool keepMET = true;
				
					if(met->Pt() < objectsToKeep[j].minPt || fabs(met->Eta()) > objectsToKeep[j].maxEta)
					{
						keepMET = false;
						if( verbosity > 1 ) cout << "skip PFMET with pT = " << met->Pt() << " and eta = " << met->Eta() << endl;
					}
					else METKept++;
				
					if( ! objectsToKeep[j].skipObjects || keepMET )
						new( (*(objectsToKeep[j].outArray))[0] ) CatMET(*met);

					if(METKept < objectsToKeep[j].minNObjects)
					{
						keepEvent = false;
						if( verbosity > 1 ) cout << "Too small number of selected PFMET: METKept = " << METKept << endl;
					}

					if( verbosity > 1 ) cout << "Processed " << objectsToKeep[j].name << endl;
					if( verbosity > 1 ) cout << "input = " << (objectsToKeep[j].inArray)->GetEntriesFast() << " output = " << (objectsToKeep[j].outArray)->GetEntriesFast() << endl;
				}
	                     else if(objectsToKeep[j].type == "cat::CatTrackMET")
                                {
                                        CatTrackMET* met;
                                        int METKept=0;
                                        met = (CatTrackMET*) (objectsToKeep[j].inArray)->At(0);
                                        bool keepMET = true;

                                        if(met->Pt() < objectsToKeep[j].minPt || fabs(met->Eta()) > objectsToKeep[j].maxEta)
                                        {
                                                keepMET = false;
                                                if( verbosity > 1 ) cout << "skip TrackMET with pT = " << met->Pt() << " and eta = " << met->Eta() << endl;
                                        }
                                        else METKept++;

                                        if( ! objectsToKeep[j].skipObjects || keepMET )
                                                new( (*(objectsToKeep[j].outArray))[0] ) CatMET(*met);

                                        if(METKept < objectsToKeep[j].minNObjects)
                                        {
                                                keepEvent = false;
                                                if( verbosity > 1 ) cout << "Too small number of selected TrackMET: METKept = " << METKept << endl;
                                        }

                                        if( verbosity > 1 ) cout << "Processed " << objectsToKeep[j].name << endl;
                                        if( verbosity > 1 ) cout << "input = " << (objectsToKeep[j].inArray)->GetEntriesFast() << " output = " << (objectsToKeep[j].outArray)->GetEntriesFast() << endl;
                                }
			
			
				else if(objectsToKeep[j].type == "cat::CatElectron")
				{
					CatElectron* electron;
					int electronsKept=0;
					for(int i=0; i<(objectsToKeep[j].inArray)->GetEntriesFast(); i++)
					{
						electron = (CatElectron*) (objectsToKeep[j].inArray)->At(i);
						bool keepElectron = true;
					
						if(electron->Pt() < objectsToKeep[j].minPt || fabs(electron->Eta()) > objectsToKeep[j].maxEta)
						{
							keepElectron = false;
							if( verbosity > 1 ) cout << "skip Electron with pT = " << electron->Pt() << " and eta = " << electron->Eta() << endl;
						}
						else electronsKept++;
						
						if( ! objectsToKeep[j].skipObjects )
							new( (*(objectsToKeep[j].outArray))[i] ) CatElectron(*electron);
						else if(keepElectron)
							new( (*(objectsToKeep[j].outArray))[electronsKept-1] ) CatElectron(*electron);
					}

					if(electronsKept < objectsToKeep[j].minNObjects)
					{
						keepEvent = false;
						if( verbosity > 1 ) cout << "Too small number of selected Electrons: electronsKept = " << electronsKept << endl;
					}

					if( verbosity > 1 ) cout << "Processed " << objectsToKeep[j].name << endl;
					if( verbosity > 1 ) cout << "input = " << (objectsToKeep[j].inArray)->GetEntriesFast() << " output = " << (objectsToKeep[j].outArray)->GetEntriesFast() << endl;
				}
			
				else if(objectsToKeep[j].type == "cat::CatMuon")
				{
					CatMuon* muon;
					int muonsKept=0;
					for(int i=0; i<(objectsToKeep[j].inArray)->GetEntriesFast(); i++)
					{
						muon = (CatMuon*) (objectsToKeep[j].inArray)->At(i);
						bool keepMuon = true;
				
						if(muon->Pt() < objectsToKeep[j].minPt || fabs(muon->Eta()) > objectsToKeep[j].maxEta)
						{
							keepMuon = false;
							if( verbosity > 1 ) cout << "skip Muon with pT = " << muon->Pt() << " and eta = " << muon->Eta() << endl;
						}
						else muonsKept++;
						
						if( ! objectsToKeep[j].skipObjects )
							new( (*(objectsToKeep[j].outArray))[i] ) CatMuon(*muon);
						else if(keepMuon)
							new( (*(objectsToKeep[j].outArray))[muonsKept-1] ) CatMuon(*muon);
					}

					if(muonsKept < objectsToKeep[j].minNObjects)
					{
						keepEvent = false;
						if( verbosity > 1 ) cout << "Too small number of selected Muons: muonsKept = " << muonsKept << endl;
					}

					if( verbosity > 1 ) cout << "Processed " << objectsToKeep[j].name << endl;
					if( verbosity > 1 ) cout << "input = " << (objectsToKeep[j].inArray)->GetEntriesFast() << " output = " << (objectsToKeep[j].outArray)->GetEntriesFast() << endl;
				}
				
				else
					cerr << "Unknown type: " << objectsToKeep[j].type << endl;
			}
		
			if(keepEvent)
			{
				nOutEvents++; // nr of output events
				
//				Calculate HLT stuff
				bool passHLT = false;
				for(unsigned int k=0; k!=hltAccept.size() ;k++)
				{
					if(outRootEvent->trigHLT(k))
					{
						hltAccept[k]++;
						passHLT = true;
					}
				}
				if(passHLT) NHLTAccept++;

				if( verbosity > 1 ) cout << "Filling the outEventTree" << endl;
			
				outEventTree->Fill();
				
				if( verbosity > 1 ) cout << "outEventTree is filled" << endl;
				
			}
		
			for(vector<keepObjects>::const_iterator iter=objectsToKeep.begin(); iter!=objectsToKeep.end(); iter++)
			{
//				if( verbosity > 1 ) cout << "(iter->outArray)->GetEntriesFast() = "  << (iter->outArray)->GetEntriesFast() << endl;
//				if( verbosity > 1 ) cout << "iter->outArray: " << iter->outArray << endl;

				( *(iter->outArray) ).Delete();

				if( verbosity > 1 ) cout << "Deleted the output TClonesArray" << endl;

			}

			if( verbosity > 1 ) cout << "Analyzing event is " << ievt << " finished!" << endl;
		
		} // loop over events
	
		if( verbosity > 1 ) cout << "Analyzing input file " << inFileName[nFile] << " finished!" << endl;

		if(inRunInfos) inRunInfos->Delete();
		inFile->Close();
		if(inFile) inFile->Delete();

	} // loop over input files
       

	cout << "Filling outRunTree" << endl;
	
	outRunTree->Fill();

	cout << "Writing to output file" << endl;

	outFile->Write();
	
	cout << "Closing output file" << endl;
	
	outFile->Close();

	delete outFile;
	delete outRunInfos;
	delete outRootEvent;

//	WARNING! Don't remove or modify the next line!!! The Automatic cat Producer depends on it!
	cout << "--> Skimmed " << nOutEvents << " out of a total of " << nInEvents << " events" << endl;

	if (((double)clock() - start) / CLOCKS_PER_SEC < 60)
		cout << "--> The skimming took " << ((double)clock() - start) / CLOCKS_PER_SEC << " seconds." << endl;
	
	else
		cout << "--> The skimming took " << (((double)clock() - start) / CLOCKS_PER_SEC)/60 << " minutes." << endl;
	
	cout << "Code running was succesfull!" << endl;
     
	return(0);
}
