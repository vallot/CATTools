#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "CATTools/DataFormats/interface/Tau.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

using namespace edm;
using namespace std;

namespace cat {

  class CATTauProducer : public edm::stream::EDProducer<> {
  public:
    explicit CATTauProducer(const edm::ParameterSet & iConfig);

    void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

  private:
    edm::EDGetTokenT<pat::TauCollection> src_;
  };

} // namespace

cat::CATTauProducer::CATTauProducer(const edm::ParameterSet & iConfig) :
  src_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("src")))
{
  produces<cat::TauCollection>();
}

void
cat::CATTauProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup)
{
  bool runOnMC_ = !iEvent.isRealData();

  Handle<pat::TauCollection> src;
  iEvent.getByToken(src_, src);

  unique_ptr<cat::TauCollection>  out(new cat::TauCollection());

  for (const pat::Tau & aPatTau : *src){
    cat::Tau aTau(aPatTau);
    if (runOnMC_){
      aTau.setGenParticleRef(aPatTau.genParticleRef());
    }
    out->push_back(aTau);
  }

  iEvent.put(std::move(out));
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATTauProducer);
