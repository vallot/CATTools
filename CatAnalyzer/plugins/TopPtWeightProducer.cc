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

  std::vector<double> paramsAny_, paramsLL_, paramsLJ_;
};

TopPtWeightProducer::TopPtWeightProducer(const edm::ParameterSet& pset)
{
  const auto srcLabel = pset.getParameter<edm::InputTag>("src");
  srcToken_ = consumes<reco::GenParticleCollection>(srcLabel);
  modesToken_ = consumes<vint>(edm::InputTag(srcLabel.label(), "modes"));
  paramsAny_ = pset.getParameter<std::vector<double>>("paramsAny");
  paramsLL_ = pset.getParameter<std::vector<double>>("paramsLL");
  paramsLJ_ = pset.getParameter<std::vector<double>>("paramsLJ");

  produces<float>();
}

void TopPtWeightProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<reco::GenParticleCollection> srcHandle;
  edm::Handle<vint> modes;
  if ( event.isRealData() or
      !event.getByToken(srcToken_, srcHandle) or !event.getByToken(modesToken_, modes) ) {
    event.put(std::auto_ptr<float>(new float(1)));
    return;
  }

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

    const auto& params = (isLep1 and isLep2) ? paramsLL_ : (isLep1 or isLep2) ? paramsLJ_ : paramsAny_;
    const double a = params[0], b = params[1];
    ptWeight = sqrt(exp(a+b*pt1)*exp(a+b*pt2));
  }

  event.put(std::auto_ptr<float>(new float(ptWeight)));
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TopPtWeightProducer);
