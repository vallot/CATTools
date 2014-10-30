/**
  \class    cat::CATElectronProducer CATElectronProducer.h "CATTools/CatProducer/interface/CATElectronProducer.h"
  \brief    CAT Electron 
*/


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
//#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
//#include "FWCore/Utilities/interface/isFinite.h"

namespace cat {

  class CATElectronProducer : public edm::EDProducer {
    public:
      explicit CATElectronProducer(const edm::ParameterSet & iConfig);
      virtual ~CATElectronProducer() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    private:
      edm::EDGetTokenT<edm::View<pat::Electron> > src_;

  };

} // namespace

cat::CATElectronProducer::CATElectronProducer(const edm::ParameterSet & iConfig) :
    src_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("src")))
{
    produces<std::vector<cat::Electron> >();
}

void 
cat::CATElectronProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;
    using namespace std;

    Handle<View<pat::Electron> > src;
    iEvent.getByToken(src_, src);

    auto_ptr<vector<cat::Electron> >  out(new vector<cat::Electron>());

    for (View<pat::Electron>::const_iterator it = src->begin(), ed = src->end(); it != ed; ++it) {
      unsigned int idx = it - src->begin();
      const pat::Electron & aPatElectron = src->at(idx);
      cat::Electron aElectron(aPatElectron);
      out->push_back(aElectron);

    }

    iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATElectronProducer);
