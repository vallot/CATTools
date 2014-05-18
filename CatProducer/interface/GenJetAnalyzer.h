#ifndef GenJetAnalyzer_h
#define GenJetAnalyzer_h

// system include files
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"

#include "../interface/CatGenJet.h"

#include "TClonesArray.h"


class GenJetAnalyzer {

public:
	GenJetAnalyzer(const edm::ParameterSet& producersNames);
	GenJetAnalyzer(const edm::ParameterSet& producersNames, int verbosity);
	GenJetAnalyzer(const edm::ParameterSet& producersNames, const edm::ParameterSet& myConfig, int verbosity);
	GenJetAnalyzer(const edm::ParameterSet& producersNames, int iter, const edm::ParameterSet& myConfig, int verbosity);
	~GenJetAnalyzer();
	void SetVerbosity(int verbosity) { verbosity_ = verbosity; };
	void Process(const edm::Event& iEvent, TClonesArray* rootGenJets);

private:
	int verbosity_;
	edm::InputTag genJetProducer_;
	edm::InputTag mcProducer_;
	std::vector<std::string> vGenJetProducer;

        std::vector<const reco::Candidate *> getAncestors(const reco::Candidate &c);
        bool hasBottom(const reco::Candidate &c);
        bool hasCharm(const reco::Candidate &c);
        bool decayFromBHadron(const reco::Candidate &c);
        bool decayFromCHadron(const reco::Candidate &c);
        const reco::Candidate* lastBHadron(const reco::Candidate &c);
        const reco::Candidate* lastCHadron(const reco::Candidate &c);

};

#endif
