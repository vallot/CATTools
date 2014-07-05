#ifndef GenEventAnalyzer_h
#define GenEventAnalyzer_h

// system include files
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"

#include "CATTools/DataFormats/interface/CatGenEvent.h"

#include "TClonesArray.h"
#include "TLorentzVector.h"

class GenEventAnalyzer{
	
public:
	GenEventAnalyzer(const edm::ParameterSet& producersNames);
	GenEventAnalyzer(const edm::ParameterSet& producersNames, int verbosity);
	GenEventAnalyzer(const edm::ParameterSet& producersNames, const edm::ParameterSet& myConfig, int verbosity);
	~GenEventAnalyzer();
	void SetVerbosity(int verbosity) {verbosity_ = verbosity; };
	void Process(const edm::Event& iEvent, TClonesArray* rootGenEvent);

private:
	int verbosity_;
	edm::InputTag genEventProducer_;
 ///
};

#endif
