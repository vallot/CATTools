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
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

using namespace edm;
using namespace std;

namespace cat {

  class CATMETProducer : public edm::EDProducer {
  public:
    explicit CATMETProducer(const edm::ParameterSet & iConfig);
    virtual ~CATMETProducer() { }

    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

  private:
    edm::InputTag src_;

  };

} // namespace

cat::CATMETProducer::CATMETProducer(const edm::ParameterSet & iConfig) :
  src_(iConfig.getParameter<edm::InputTag>( "src" ))
{
  produces<std::vector<cat::MET> >();
}

void 
cat::CATMETProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup)
{
  Handle<View<pat::MET> > src;
  iEvent.getByLabel(src_, src);

  auto_ptr<vector<cat::MET> >  out(new vector<cat::MET>());

  for (const pat::MET & aPatMET : *src) {
    cat::MET aMET(aPatMET);
    out->push_back(aMET);
  }

  iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATMETProducer);
