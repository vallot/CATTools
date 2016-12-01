#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

namespace cat {

  class CATVertexProducer : public edm::stream::EDProducer<> {
  public:
    explicit CATVertexProducer(const edm::ParameterSet & iConfig);
    virtual ~CATVertexProducer() { }

    void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

  private:
    const edm::EDGetTokenT<reco::VertexCollection> vertexLabel_;
    const unsigned int minNDOF_;
    const double maxAbsZ_;
    const double maxd0_;

  };

} // namespace

cat::CATVertexProducer::CATVertexProducer(const edm::ParameterSet & iConfig) :
  vertexLabel_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexLabel"))),
  minNDOF_(iConfig.getParameter<unsigned int>("minimumNDOF")),
  maxAbsZ_(iConfig.getParameter<double>("maxAbsZ")),
  maxd0_(iConfig.getParameter<double>("maxd0"))
{
  produces<reco::VertexCollection >();
  produces<int >("nGoodPV");
  produces<int >("nPV");
}

void cat::CATVertexProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup)
{
  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(vertexLabel_,recVtxs);

  const int nPV = recVtxs->size();
  int nGoodPV = 0;
  std::auto_ptr<reco::VertexCollection> out(new reco::VertexCollection());

  // only save the first good vertex!
  for ( auto &vtx : *recVtxs ){
    if ( vtx.ndof() > minNDOF_ &&
        ( (maxAbsZ_ <=0 ) || std::abs(vtx.z()) <= maxAbsZ_ ) &&
        ( (maxd0_ <=0 ) || std::abs(vtx.position().rho()) <= maxd0_ ) &&
        !(vtx.isFake() ) ){
      if (nGoodPV == 0) {
        out->push_back(vtx);
        out->back().removeTracks();
      }
      ++nGoodPV;
    }
  }

  iEvent.put(std::auto_ptr<int>(new int (nGoodPV)), "nGoodPV");
  iEvent.put(std::auto_ptr<int>(new int (nPV)), "nPV");
  iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATVertexProducer);
