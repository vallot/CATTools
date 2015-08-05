#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "CATTools/DataFormats/interface/VertexInfo.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

using namespace edm;
using namespace std;
using namespace reco;

namespace cat {

  class CATVertexInfoProducer : public edm::EDProducer {
  public:
    explicit CATVertexInfoProducer(const edm::ParameterSet & iConfig);
    virtual ~CATVertexInfoProducer() { }

    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

  private:
    cat::VertexInfo *out_;

    edm::EDGetTokenT<reco::VertexCollection> vertexLabel_;

  };

} // namespace

cat::CATVertexInfoProducer::CATVertexInfoProducer(const edm::ParameterSet & iConfig) :
  vertexLabel_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexLabel")))
{
  produces<cat::VertexInfo >();
}

void 
cat::CATVertexInfoProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) 
{
  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(vertexLabel_,recVtxs);

  out_ = new cat::VertexInfo();

  out_->setnumberOfVertex(recVtxs->size());
  if (!recVtxs->empty()){
    reco::Vertex pv = recVtxs->at(0);
    out_->setposition(pv.position());
    out_->setvalidity(pv.isValid());
    out_->setisFake(pv.isFake());
    out_->setndof(pv.ndof());
    out_->setchi2(pv.chi2());
  }
  auto_ptr<cat::VertexInfo> out(out_);
  iEvent.put(out); 
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATVertexInfoProducer);
