#ifndef TrackMETAnalyzer_h
#define TrackMETAnalyzer_h

// system include files
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

//#include "DataFormats/METReco/interface/MET.h"
//#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "../interface/CatEvent.h"
#include "../interface/CatTrackMET.h"

#include "../interface/RecoMETAnalyzer.h"

#include "TClonesArray.h"


class TrackMETAnalyzer{
	
public:
	TrackMETAnalyzer(const edm::ParameterSet& producersNames);
	TrackMETAnalyzer(const edm::ParameterSet& producersNames, const edm::ParameterSet& myConfig, int verbosity);
	TrackMETAnalyzer(const edm::ParameterSet& producersNames, int iter, const edm::ParameterSet& myConfig, int verbosity);
	~TrackMETAnalyzer();
	void SetVerbosity(int verbosity) {verbosity_ = verbosity; };
	void Process(const edm::Event& iEvent, TClonesArray* rootMET);

private:
	int verbosity_;
	edm::InputTag metProducer_;
	bool useMC_;
  std::vector<std::string> vTrackMETProducer;
	RecoMETAnalyzer* myMETAnalyzer;
};

#endif

