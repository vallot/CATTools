/**
  \class    cat::CATMETProducer CATMETProducer.h "CATTools/CatProducer/interface/CATMETProducer.h"
  \brief    CAT MET 
*/


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "CATTools/DataFormats/interface/MET.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
//#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
//#include "FWCore/Utilities/interface/isFinite.h"

namespace cat {

  class CATMETProducer : public edm::EDProducer {
    public:
      explicit CATMETProducer(const edm::ParameterSet & iConfig);
      virtual ~CATMETProducer() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    private:
      edm::EDGetTokenT<edm::View<pat::MET> > src_;

  };

} // namespace

cat::CATMETProducer::CATMETProducer(const edm::ParameterSet & iConfig) :
    src_(consumes<edm::View<pat::MET> >(iConfig.getParameter<edm::InputTag>("src")))
{
    produces<std::vector<cat::MET> >();
}

void 
cat::CATMETProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;
    using namespace std;

    Handle<View<pat::MET> > src;
    iEvent.getByToken(src_, src);

    auto_ptr<vector<cat::MET> >  out(new vector<cat::MET>());

    for (View<pat::MET>::const_iterator it = src->begin(), ed = src->end(); it != ed; ++it) {
      unsigned int idx = it - src->begin();
      const pat::MET & aPatMET = src->at(idx);
      cat::MET aMET(aPatMET);
      out->push_back(aMET);

    }

    iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATMETProducer);
