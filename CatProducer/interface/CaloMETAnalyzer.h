#ifndef CaloMETAnalyzer_h
#define CaloMETAnalyzer_h

// system include files
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

//#include "DataFormats/METReco/interface/MET.h"
//#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "../interface/CatEvent.h"
#include "../interface/CatCaloMET.h"

#include "../interface/METAnalyzer.h"

#include "TClonesArray.h"


class CaloMETAnalyzer{
	
public:
	CaloMETAnalyzer(const edm::ParameterSet& producersNames);
	CaloMETAnalyzer(const edm::ParameterSet& producersNames, const edm::ParameterSet& myConfig, int verbosity);
	~CaloMETAnalyzer();
	void SetVerbosity(int verbosity) {verbosity_ = verbosity; };
	void Process(const edm::Event& iEvent, TClonesArray* rootMET);

private:
	int verbosity_;
	edm::InputTag metProducer_;
	bool useMC_;

	METAnalyzer* myMETAnalyzer;
};

#endif

