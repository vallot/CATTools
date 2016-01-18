#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "CATTools/CommonTools/interface/TTbarModeDefs.h"

using namespace std;
using namespace cat;

class TopPtWeightProducer: public edm::stream::EDProducer<>
{
public:
  TopPtWeightProducer(const edm::ParameterSet& pset);
  void produce(edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  typedef std::vector<int> vint;

  edm::EDGetTokenT<reco::GenParticleCollection> srcToken_;
  edm::EDGetTokenT<vint> modesToken_;
};

TopPtWeightProducer::TopPtWeightProducer(const edm::ParameterSet& pset)
{
  const auto srcLabel = pset.getParameter<edm::InputTag>("src");
  srcToken_ = consumes<reco::GenParticleCollection>(srcLabel);
  modesToken_ = consumes<vint>(edm::InputTag(srcLabel.label(), "modes"));

  produces<float>();
}

void TopPtWeightProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<reco::GenParticleCollection> srcHandle;
  event.getByToken(srcToken_, srcHandle);

  edm::Handle<vint> modes;
  event.getByToken(modesToken_, modes);

  std::vector<const reco::GenParticle*> tQuarks;
  for ( auto& p : *srcHandle ) {
    if ( abs(p.pdgId()) == 6 ) tQuarks.push_back(&p);
  }

  double ptWeight = 1;
  if ( tQuarks.size() == 2 and modes->size() == 2 ) {

    const double pt1 = tQuarks.at(0)->pt();
    const double pt2 = tQuarks.at(1)->pt();

    const int mode1 = modes->at(0), mode2 = modes->at(1);
    const bool isLep1 = (mode1 == CH_MUON or mode2 == CH_ELECTRON);
    const bool isLep2 = (mode2 == CH_MUON or mode2 == CH_ELECTRON);

    double a = 0.156, b = -0.00137;
    if ( isLep1 and isLep2 ) { a = 0.148; b = -0.00129; }
    else if ( isLep1 or isLep2 ) { a = 0.159; b = -0.00141; }

    ptWeight = sqrt(exp(a+b*pt1)*exp(a+b*pt2));
  }

  event.put(std::auto_ptr<float>(new float(ptWeight)), "ptWeight");
}
