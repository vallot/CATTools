#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "CATTools/DataFormats/interface/MCParticle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
//#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
//#include "FWCore/Utilities/interface/isFinite.h"

using namespace edm;
using namespace std;

namespace cat {

  class CATMCParticleProducer : public edm::EDProducer {
  public:
    explicit CATMCParticleProducer(const edm::ParameterSet & iConfig);
    virtual ~CATMCParticleProducer() { }

    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

  private:
    edm::InputTag src_;

    const double pt_;
    const double eta_;

  };

} // namespace

cat::CATMCParticleProducer::CATMCParticleProducer(const edm::ParameterSet & iConfig) :
  src_(iConfig.getParameter<edm::InputTag>( "src" )),
  pt_(iConfig.getParameter<double>("pt")),
  eta_(iConfig.getParameter<double>("eta"))
{
  produces<std::vector<cat::MCParticle> >();
}

void 
cat::CATMCParticleProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) 
{
  Handle<View<reco::GenParticle> > genParticles;
  iEvent.getByLabel(src_,genParticles);
    
  auto_ptr<vector<cat::MCParticle> >  out(new vector<cat::MCParticle>());

  for (View<reco::GenParticle>::const_iterator it = genParticles->begin(), ed = genParticles->end(); it != ed; ++it) {
    unsigned int idx = it - genParticles->begin();
    const reco::GenParticle & aGenParticle = genParticles->at(idx);

    if ( aGenParticle.pt() < pt_ || fabs(aGenParticle.eta()) > eta_ ) continue;  

    cat::MCParticle aMCParticle(aGenParticle);

    out->push_back(aMCParticle);

  }

  iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATMCParticleProducer);
