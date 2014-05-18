#ifndef ElectronAnalyzer_h
#define ElectronAnalyzer_h

// system include files
#include <iostream>
#include <TLorentzVector.h>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
//#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
//#include "TrackingTools/IPTools/interface/IPTools.h"
//#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "../interface/CatElectron.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
//#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"

#include "TClonesArray.h"

using namespace edm;
using namespace reco;
class ElectronAnalyzer{
	
public:
	ElectronAnalyzer(const edm::ParameterSet& producersNames);
	ElectronAnalyzer(const edm::ParameterSet& producersNames, const edm::ParameterSet& myConfig, int verbosity);
	ElectronAnalyzer(const edm::ParameterSet& producersNames, int iter, const edm::ParameterSet& myConfig, int verbosity);
	~ElectronAnalyzer();
	void SetVerbosity(int verbosity) {verbosity_ = verbosity; };
	void Process(const edm::Event& iEvent, TClonesArray* rootElectrons, const edm::EventSetup& iSetup);
  typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;

private:
	int verbosity_;
	edm::InputTag electronProducer_;
	std::vector<std::string> vElectronProducer;	
	edm::InputTag primaryVertexProducer_;
	edm::InputTag TrackLabel_;
	bool useMC_;
	bool runSuperCluster_;
	bool newId_;
	bool doPrimaryVertex_;
	bool isData_;

};

#endif
