#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "CATTools/DataFormats/interface/Photon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "FWCore/Utilities/interface/isFinite.h"

using namespace edm;
using namespace std;

namespace cat {

  class CATPhotonProducer : public edm::stream::EDProducer<> {
  public:
    explicit CATPhotonProducer(const edm::ParameterSet & iConfig);
    virtual ~CATPhotonProducer() { }

    void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

  private:
    edm::EDGetTokenT<pat::PhotonCollection> src_;
  };

} // namespace

cat::CATPhotonProducer::CATPhotonProducer(const edm::ParameterSet & iConfig) :
  src_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("src")))
{
  produces<std::vector<cat::Photon> >();
}

void 
cat::CATPhotonProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup)
{
  bool runOnMC_ = !iEvent.isRealData();

  Handle<pat::PhotonCollection> src;
  iEvent.getByToken(src_, src);

  auto_ptr<vector<cat::Photon> >  out(new vector<cat::Photon>());

  for (const pat::Photon & aPatPhoton : *src){
    cat::Photon aPhoton(aPatPhoton);
    if (runOnMC_){
      aPhoton.setGenParticleRef(aPatPhoton.genParticleRef());
    }
    out->push_back(aPhoton);
  }

  iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATPhotonProducer);
