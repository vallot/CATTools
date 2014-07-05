#ifndef JetAnalyzer_h
#define JetAnalyzer_h

// system include files
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"

#include "RecoBTag/Records/interface/BTagPerformanceRecord.h"
#include "CondFormats/PhysicsToolsObjects/interface/BinningPointByMap.h"
#include "RecoBTag/PerformanceDB/interface/BtagPerformance.h"

#include "CATTools/DataFormats/interface/Jet.h"

using namespace cat;

class JetAnalyzer {
	
public:
	JetAnalyzer();
	JetAnalyzer(int verbosity);
	JetAnalyzer(const edm::ParameterSet& myConfig, int verbosity);
	~JetAnalyzer();
	void SetVerbosity(int verbosity) {verbosity_ = verbosity; };
	Jet Process(const reco::Jet* jet, const edm::EventSetup& iSetup);

private:
	int verbosity_;
	bool useMC_;
	bool isData_;

};

#endif
