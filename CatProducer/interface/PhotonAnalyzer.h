#ifndef PhotonAnalyzer_h
#define PhotonAnalyzer_h

// system include files
#include <iostream>
#include <Math/VectorUtil.h>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "CATTools/DataFormats/interface/Photon.h"

#include "TClonesArray.h"


class PhotonAnalyzer
{
	
public:
	PhotonAnalyzer(const edm::ParameterSet& producersNames);
	PhotonAnalyzer(const edm::ParameterSet& producersNames, const edm::ParameterSet& myConfig, int verbosity);
	PhotonAnalyzer(const edm::ParameterSet& producersNames, int iter, const edm::ParameterSet& myConfig, int verbosity);
	~PhotonAnalyzer();
	void SetVerbosity(int verbosity) { verbosity_ = verbosity; };
	void Process(const edm::Event& iEvent, TClonesArray* rootMuons, const edm::EventSetup& iSetup);

        typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > > IsoDepositMaps;
        typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;

private:
	int verbosity_;
	edm::InputTag photonProducer_;
	std::vector<std::string> vPhotonProducer;
	bool useMC_;
 
        //photon isolation
        std::vector<edm::InputTag> inputTagIsoDepPhotons_;
        std::vector<edm::InputTag> inputTagIsoValPhotons_; 
};

#endif
