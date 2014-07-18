/**
  \class    cat::CATJetProducer CATJetProducer.h "CATTools/CatProducer/interface/CATJetProducer.h"
  \brief    CAT Jet 
*/


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "FWCore/Utilities/interface/isFinite.h"

namespace cat {

  class CATJetProducer : public edm::EDProducer {
    public:
      explicit CATJetProducer(const edm::ParameterSet & iConfig);
      virtual ~CATJetProducer() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    private:
      edm::EDGetTokenT<edm::View<pat::Jet> > src_;

  };

} // namespace

cat::CATJetProducer::CATJetProducer(const edm::ParameterSet & iConfig) :
    src_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("src")))
{
    produces<std::vector<cat::Jet> >();
}

void 
cat::CATJetProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;
    using namespace std;

    Handle<View<pat::Jet> > src;
    iEvent.getByToken(src_, src);

    auto_ptr<vector<cat::Jet> >  out(new vector<cat::Jet>());

    for (View<pat::Jet>::const_iterator it = src->begin(), ed = src->end(); it != ed; ++it) {
      unsigned int idx = it - src->begin();
      const pat::Jet & aPatJet = src->at(idx);
      cat::Jet aJet(aPatJet);
      out->push_back(aJet);

    }

    iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATJetProducer);
