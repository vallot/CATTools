#ifndef SpinCorrGenAnalyzer_h
#define SpinCorrGenAnalyzer_h

// system include files
#include <iostream>
#include <Math/VectorUtil.h>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/Vector3D.h"

#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"

#include "../interface/CatSpinCorrGen.h"

#include "TClonesArray.h"
#include "TLorentzVector.h"

class SpinCorrGenAnalyzer{
	
public:
	SpinCorrGenAnalyzer(const edm::ParameterSet& producersNames);
	SpinCorrGenAnalyzer(const edm::ParameterSet& producersNames, int verbosity);
	SpinCorrGenAnalyzer(const edm::ParameterSet& producersNames, const edm::ParameterSet& myConfig, int verbosity);
	~SpinCorrGenAnalyzer();
	void SetVerbosity(int verbosity) {verbosity_ = verbosity; };
	void Process(const edm::Event& iEvent, TClonesArray* rootSpinCorrGen);

private:
	int verbosity_;
	edm::InputTag genEventProducer_;
 ///
};

#endif
