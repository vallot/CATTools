#ifndef NPGenEventAnalyzer_h
#define NPGenEventAnalyzer_h

// system include files
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "CATTools/DataFormats/interface/NPGenEvent.h"

#include "TClonesArray.h"
#include "TLorentzVector.h"

using namespace cat;

class NPGenEventAnalyzer{
	
public:
	NPGenEventAnalyzer(const edm::ParameterSet& producersNames);
	NPGenEventAnalyzer(const edm::ParameterSet& producersNames, int verbosity);
	NPGenEventAnalyzer(const edm::ParameterSet& producersNames, const edm::ParameterSet& myConfig, int verbosity);
	~NPGenEventAnalyzer();
	void SetVerbosity(int verbosity) {verbosity_ = verbosity; };
	void Process(const edm::Event& iEvent, TClonesArray* rootGenEvent);
        TLorentzVector P4toTLV (reco::Particle::LorentzVector a);
        MCParticle ConvertMCPart(reco::GenParticleCollection::const_iterator t);
        MCParticle ConvertMCPart(reco::GenParticle::const_iterator t);
  
private:
	int verbosity_;
	edm::InputTag genParticlesProducer_;
 ///
};

#endif
