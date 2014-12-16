#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "CATTools/DataFormats/interface/Photon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

using namespace edm;
using namespace std;

namespace cat {

  class CATPhotonProducer : public edm::EDProducer {
  public:
    explicit CATPhotonProducer(const edm::ParameterSet & iConfig);
    virtual ~CATPhotonProducer() { }

    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

  private:
    edm::InputTag src_;
  };

} // namespace

cat::CATPhotonProducer::CATPhotonProducer(const edm::ParameterSet & iConfig) :
  src_(iConfig.getParameter<edm::InputTag>( "src" ))
{
  produces<std::vector<cat::Photon> >();
}

void 
cat::CATPhotonProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup)
{
  Handle<pat::PhotonCollection> src;
  iEvent.getByLabel(src_, src);

  auto_ptr<vector<cat::Photon> >  out(new vector<cat::Photon>());

  for (const pat::Photon & aPatPhoton : *src){
    cat::Photon aPhoton(aPatPhoton);
    out->push_back(aPhoton);
  }

  iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATPhotonProducer);
