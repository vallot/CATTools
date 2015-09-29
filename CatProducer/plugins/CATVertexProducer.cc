#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

using namespace edm;
using namespace std;
using namespace reco;

namespace cat {

  class CATVertexProducer : public edm::stream::EDProducer<> {
  public:
    explicit CATVertexProducer(const edm::ParameterSet & iConfig);
    virtual ~CATVertexProducer() { }

    void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

  private:
    reco::VertexCollection *out_;

    edm::EDGetTokenT<reco::VertexCollection> vertexLabel_;
    unsigned int minNDOF;
    double maxAbsZ;
    double maxd0;

  };

} // namespace

cat::CATVertexProducer::CATVertexProducer(const edm::ParameterSet & iConfig) :
  vertexLabel_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexLabel")))
{
  minNDOF = iConfig.getParameter<unsigned int>("minimumNDOF");
  maxAbsZ = iConfig.getParameter<double>("maxAbsZ");
  maxd0 = iConfig.getParameter<double>("maxd0");
  produces<reco::VertexCollection >();
  produces<int >("nGoodPV");
  produces<int >("nPV");
}

void
cat::CATVertexProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup)
{
  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(vertexLabel_,recVtxs);

  int nPV = recVtxs->size();
  int nGoodPV = 0;
  out_ = new reco::VertexCollection();

  // only save the first good vertex!
  for (auto &vtx : *recVtxs){
    if ( vtx.ndof() > minNDOF &&
	 ( (maxAbsZ <=0 ) || std::abs(vtx.z()) <= maxAbsZ ) &&
	 ( (maxd0 <=0 ) || std::abs(vtx.position().rho()) <= maxd0 ) &&
	 !(vtx.isFake() ) ){
      if (nGoodPV == 0){
	out_->push_back(vtx);
	out_->back().removeTracks();
      }
      nGoodPV++;
    }
  }

  iEvent.put(std::auto_ptr<int>(new int (nGoodPV)), "nGoodPV");
  iEvent.put(std::auto_ptr<int>(new int (nPV)), "nPV");
  auto_ptr<reco::VertexCollection> out(out_);
  iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATVertexProducer);
