/**
  \class    cat::CATMuonProducer CATMuonProducer.h "CATTools/CatProducer/interface/CATMuonProducer.h"
  \brief    CAT Muon 
*/


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "CATTools/DataFormats/interface/Muon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "FWCore/Utilities/interface/isFinite.h"

namespace cat {

  class CATMuonProducer : public edm::EDProducer {
    public:
      explicit CATMuonProducer(const edm::ParameterSet & iConfig);
      virtual ~CATMuonProducer() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    private:
      edm::EDGetTokenT<edm::View<pat::Muon> > src_;
      edm::EDGetTokenT<edm::View<reco::Vertex> > vertexLabel_;

  };

} // namespace

cat::CATMuonProducer::CATMuonProducer(const edm::ParameterSet & iConfig) :
    src_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("src"))),
    vertexLabel_(consumes<edm::View<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertexLabel")))
{
    produces<std::vector<cat::Muon> >();
}

void 
cat::CATMuonProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;
    using namespace std;

    Handle<View<pat::Muon> > src;
    iEvent.getByToken(src_, src);

    Handle<View<reco::Vertex> > recVtxs;
    iEvent.getByToken(vertexLabel_,recVtxs);

    reco::Vertex pv = recVtxs->at(0);
    
    auto_ptr<vector<cat::Muon> >  out(new vector<cat::Muon>());

    for (View<pat::Muon>::const_iterator it = src->begin(), ed = src->end(); it != ed; ++it) {
      unsigned int idx = it - src->begin();
      const pat::Muon & aPatMuon = src->at(idx);
      cat::Muon aMuon(aPatMuon);
  
      aMuon.setChargedHadronIso( aPatMuon.chargedHadronIso() );
      aMuon.setNeutralHadronIso( aPatMuon.neutralHadronIso() );
      aMuon.setPhotonIso( aPatMuon.photonIso() );
      aMuon.setPUChargedHadronIso( aPatMuon.puChargedHadronIso() );

      aMuon.setIsTightMuon( aPatMuon.isTightMuon(pv) );
      aMuon.setIsLooseMuon( aPatMuon.isLooseMuon() );
      aMuon.setIsSoftMuon( aPatMuon.isSoftMuon(pv) );
    
      out->push_back(aMuon);

    }

    iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATMuonProducer);
