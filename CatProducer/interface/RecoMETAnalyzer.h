#ifndef RecoMETAnalyzer_h
#define RecoMETAnalyzer_h

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
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "../interface/CatEvent.h"
#include "../interface/CatMET.h"

#include "TClonesArray.h"


class RecoMETAnalyzer{
	
public:
	RecoMETAnalyzer(const edm::ParameterSet& producersNames);
	RecoMETAnalyzer(const edm::ParameterSet& producersNames, const edm::ParameterSet& myConfig, int verbosity);
	~RecoMETAnalyzer();
	void SetVerbosity(int verbosity) {verbosity_ = verbosity; };
	cat::CatMET Process(const reco::Candidate* met);

private:
	int verbosity_;
	edm::InputTag metProducer_;
	bool useMC_;
};

#endif
