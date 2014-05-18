#ifndef CaloJetAnalyzer_h
#define CaloJetAnalyzer_h

// system include files
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"

#include "../interface/CatCaloJet.h"
#include "../interface/JetAnalyzer.h"

#include "TClonesArray.h"


class CaloJetAnalyzer : JetAnalyzer {

public:
	CaloJetAnalyzer(const edm::ParameterSet& producersNames);
	CaloJetAnalyzer(const edm::ParameterSet& producersNames, int verbosity);
	CaloJetAnalyzer(const edm::ParameterSet& producersNames, const edm::ParameterSet& myConfig, int verbosity);
	CaloJetAnalyzer(const edm::ParameterSet& producersNames, int iter, const edm::ParameterSet& myConfig, int verbosity);
	~CaloJetAnalyzer();
	void SetVerbosity(int verbosity) {verbosity_ = verbosity; };
	void Process(const edm::Event& iEvent, TClonesArray* rootCaloJets, const edm::EventSetup& iSetup);

private:
	int verbosity_;
	edm::InputTag caloJetProducer_;
	edm::InputTag mcProducer_;
	bool doCaloJetId_;
	std::vector<std::string> vCaloJetProducer;
	JetAnalyzer* myJetAnalyzer; // FIXME: Handle the deletion of the JetAnalyzer

};

#endif
