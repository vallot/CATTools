#ifndef JPTJetAnalyzer_h
#define JPTJetAnalyzer_h

// system include files
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/JetReco/interface/JPTJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"

#include "../interface/CatJPTJet.h"
#include "../interface/JetAnalyzer.h"

#include "TClonesArray.h"


class JPTJetAnalyzer : JetAnalyzer {

public:
	JPTJetAnalyzer(const edm::ParameterSet& producersNames);
	JPTJetAnalyzer(const edm::ParameterSet& producersNames, int verbosity);
	JPTJetAnalyzer(const edm::ParameterSet& producersNames, const edm::ParameterSet& myConfig, int verbosity);
	JPTJetAnalyzer(const edm::ParameterSet& producersNames, int iter, const edm::ParameterSet& myConfig, int verbosity);
	~JPTJetAnalyzer();
	void SetVerbosity(int verbosity) {verbosity_ = verbosity; };
	void Process(const edm::Event& iEvent, TClonesArray* rootJPTJets, const edm::EventSetup& iSetup);

private:
	int verbosity_;
	edm::InputTag JPTJetProducer_;
	edm::InputTag mcProducer_;
	bool doJPTJetId_;
	std::vector<std::string> vJPTJetProducer;
	JetAnalyzer* myJetAnalyzer; // FIXME: Handle the deletion of the JetAnalyzer

};

#endif
