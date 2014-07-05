#ifndef MCAnalyzer_h
#define MCAnalyzer_h

// system include files
#include <memory>
#include <string>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruthFinder.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruth.h"
#include "RecoEgamma/EgammaMCTools/interface/ElectronMCTruth.h"
//#include "DataFormats/EgammaCandidates/interface/ConvertedPhoton.h"

#include "../interface/ParticleTreeDrawer.h"
#include "CATTools/DataFormats/interface/Event.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "CATTools/DataFormats/interface/MCParticle.h"
#include "CATTools/DataFormats/interface/Jet.h"

#include "TClonesArray.h"

using namespace cat;

class MCAnalyzer{
	
public:
	MCAnalyzer();
	MCAnalyzer(const edm::ParameterSet& config, const edm::ParameterSet& producersNames);
	~MCAnalyzer();
	void SetVerbosity(int verbosity) {verbosity_ = verbosity; };
	void DrawMCTree(const edm::Event& iEvent, const edm::EventSetup& iSetup, const edm::ParameterSet& config, const edm::ParameterSet& producersNames);
	void PDFInfo(const edm::Event& iEvent, Event* rootEvent);
  void ProcessMCParticle(const edm::Event& iEvent, TClonesArray* rootMCParticles);	

private:

	int verbosity_;
	
	bool doElectronMC_;
	double electronMC_etaMax_;
	double electronMC_ptMin_;
	bool doMuonMC_;
	double muonMC_etaMax_;
	double muonMC_ptMin_;
	bool doJetMC_;
	double jetMC_etaMax_;
	double jetMC_ptMin_;
	bool doUnstablePartsMC_;
	bool doMETMC_;

	std::string signalGenerator_;
        edm::InputTag genParticlesProducer_;
};

#endif
