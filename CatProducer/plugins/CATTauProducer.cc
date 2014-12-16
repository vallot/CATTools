#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "CATTools/DataFormats/interface/Tau.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

using namespace edm;
using namespace std;

namespace cat {

  class CATTauProducer : public edm::EDProducer {
  public:
    explicit CATTauProducer(const edm::ParameterSet & iConfig);
    virtual ~CATTauProducer() { }

    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

  private:
    edm::EDGetTokenT<pat::TauCollection> src_;
  };

} // namespace

cat::CATTauProducer::CATTauProducer(const edm::ParameterSet & iConfig) :
  src_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("src")))
{
  produces<std::vector<cat::Tau> >();
}

void 
cat::CATTauProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup)
{
  Handle<pat::TauCollection> src;
  iEvent.getByToken(src_, src);

  auto_ptr<vector<cat::Tau> >  out(new vector<cat::Tau>());

  for (const pat::Tau & aPatTau : *src){
    cat::Tau aTau(aPatTau);
    out->push_back(aTau);
  }

  iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATTauProducer);
