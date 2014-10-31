/**
   \class    cat::CATPhotonProducer CATPhotonProducer.h "CATTools/CatProducer/interface/CATPhotonProducer.h"
   \brief    CAT Photon 
*/


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
//#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
//#include "FWCore/Utilities/interface/isFinite.h"

namespace cat {

  class CATPhotonProducer : public edm::EDProducer {
  public:
    explicit CATPhotonProducer(const edm::ParameterSet & iConfig);
    virtual ~CATPhotonProducer() { }

    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

  private:
    //   edm::EDGetTokenT<edm::View<pat::Photon> > src_;
    edm::InputTag src;

  };

} // namespace

cat::CATPhotonProducer::CATPhotonProducer(const edm::ParameterSet & iConfig) :
  //  src_(consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("src")))
  src_(iConfig.getParameter<edm::InputTag>( "src" ))
{
  produces<std::vector<cat::Photon> >();
}

void 
cat::CATPhotonProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
  using namespace edm;
  using namespace std;

  Handle<View<pat::Photon> > src;
  iEvent.getByLabel(src_, src);

  auto_ptr<vector<cat::Photon> >  out(new vector<cat::Photon>());

  for (View<pat::Photon>::const_iterator it = src->begin(), ed = src->end(); it != ed; ++it) {
    unsigned int idx = it - src->begin();
    const pat::Photon & aPatPhoton = src->at(idx);
    cat::Photon aPhoton(aPatPhoton);
    out->push_back(aPhoton);

  }

  iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATPhotonProducer);
