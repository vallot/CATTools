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

  // only save the first vertex!
  // If the first vertex fails quality cut, put it in the 2nd place.
  int igoodpv0 = -1;
  for ( int i=0; i<nPV; ++i ) {
    const auto& vtx = recVtxs->at(i);
    if ( vtx.ndof() <= minNDOF_ ) continue;
    if ( maxAbsZ_ > 0 and std::abs(vtx.z()) > maxAbsZ_ ) continue;
    if ( maxd0_ > 0 and std::abs(vtx.position().rho()) > maxd0_ ) continue;
    if ( vtx.isFake() ) continue;

    if ( igoodpv0 < 0 ) igoodpv0 = i;
    ++nGoodPV;
  }

  std::unique_ptr<reco::VertexCollection> out(new reco::VertexCollection());
  if ( igoodpv0 >= 0 ) {
    out->push_back(recVtxs->at(igoodpv0));
    out->back().removeTracks();
  }
  if ( nPV > 0 and igoodpv0 != 0 ) {
    out->push_back(recVtxs->front());
    out->back().removeTracks();
  }

  iEvent.put(std::move(std::unique_ptr<int>(new int (nGoodPV))), "nGoodPV");
  iEvent.put(std::move(std::unique_ptr<int>(new int (nPV))), "nPV");
  iEvent.put(std::move(out));
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATVertexProducer);
