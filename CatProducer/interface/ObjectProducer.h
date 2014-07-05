#ifndef ObjectProducer_h
#define ObjectProducer_h

// system include files
#include <memory>
#include <string>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "../interface/HLTAnalyzer.h"
#include "../interface/MCAnalyzer.h"
#include "../interface/MCAssociator.h"
#include "../interface/VertexAnalyzer.h"
#include "../interface/JetAnalyzer.h"
#include "../interface/GenJetAnalyzer.h"
#include "../interface/PFJetAnalyzer.h"
#include "../interface/MuonAnalyzer.h"
#include "../interface/ElectronAnalyzer.h"
#include "../interface/PhotonAnalyzer.h"
#include "../interface/METAnalyzer.h"
#include "../interface/PFMETAnalyzer.h"
#include "../interface/GenEventAnalyzer.h"
#include "../interface/NPGenEventAnalyzer.h"
#include "../interface/SpinCorrGenAnalyzer.h"

#include "CATTools/DataFormats/interface/Run.h"
#include "CATTools/DataFormats/interface/Event.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "CATTools/DataFormats/interface/MCParticle.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/GenJet.h"
#include "CATTools/DataFormats/interface/PFJet.h"
#include "CATTools/DataFormats/interface/Lepton.h"
#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/MET.h"
#include "CATTools/DataFormats/interface/PFMET.h"
#include "CATTools/DataFormats/interface/GenEvent.h"
#include "CATTools/DataFormats/interface/NPGenEvent.h"
#include "CATTools/DataFormats/interface/SpinCorrGen.h"
#include "CATTools/DataFormats/interface/Vertex.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TClonesArray.h"

//using namespace cat;

class ObjectProducer : public edm::EDAnalyzer {
public:
	explicit ObjectProducer(const edm::ParameterSet&);
	~ObjectProducer();
	
	
private:
	virtual void beginJob() ;
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endLuminosityBlock(const edm::LuminosityBlock&, const EventSetup&);
	virtual void endJob() ;

	edm::ParameterSet myConfig_;
	edm::ParameterSet producersNames_;
	
	int verbosity;
	std::string rootFileName_ ;
	TFile* rootFile_ ;
	TTree* eventTree_;
	TTree* runTree_;
        TH1F* tmp_;
	bool doHLT;
	bool doMC;
	bool doPDFInfo;
	bool doSignalMuMuGamma;
	bool doSignalTopTop;
	bool doPrimaryVertex;
	bool runGeneralTracks;
	bool doGenJet;
	bool doPFJet;
	bool doMuon;
	bool doElectron;
	bool doPhoton;
	bool doPFMET;
	bool doGenEvent;
	bool doNPGenEvent;
	bool doSpinCorrGen;
	bool drawMCTree;
	std::vector<std::string> vGenJetProducer;
	std::vector<std::string> vPFJetProducer;
	std::vector<std::string> vMuonProducer;
	std::vector<std::string> vElectronProducer;
	std::vector<std::string> vPhotonProducer;
        std::vector<std::string> vPFmetProducer;
        std::vector<std::string> vTrackmetProducer; 
	int nTotEvt_;
	HLTAnalyzer* hltAnalyzer_;
	cat::Run* runInfos_;
	cat::Event* rootEvent;
	TClonesArray* mcParticles;
	TClonesArray* tracks;
	std::vector<TClonesArray*> vcaloJets;
	std::vector<TClonesArray*> vgenJets;
	std::vector<TClonesArray*> vpfJets;
	std::vector<TClonesArray*> vjptJets;
	std::vector<TClonesArray*> vmuons;
	std::vector<TClonesArray*> velectrons;
	std::vector<TClonesArray*> vphotons;
	TClonesArray* CALOmet;
	std::vector<TClonesArray*> vPFmets;
	std::vector<TClonesArray*> vTrackmets;
	TClonesArray* TCmet;
	TClonesArray* genEvent;
	TClonesArray* NPgenEvent;
	TClonesArray* spinCorrGen;
	TClonesArray* primaryVertex;

        bool useEventCounter_;
        std::vector<std::string> filters_;

        bool isRealData_;
};

#endif
