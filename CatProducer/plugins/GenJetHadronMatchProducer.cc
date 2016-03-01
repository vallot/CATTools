#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CATTools/DataFormats/interface/GenJet.h"
#include "CATTools/DataFormats/interface/MCParticle.h"

#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"

using namespace edm;
using namespace std;

namespace cat {

class GenJetHadronMatchProducer : public edm::stream::EDProducer<>
{
public:
  GenJetHadronMatchProducer(const edm::ParameterSet& pset);
  virtual ~GenJetHadronMatchProducer() { }

  void produce(edm::Event & event, const edm::EventSetup&) override;

private:
  void collectAncestorPartons(const reco::Candidate* cand, std::set<int>& partonIds) const;

  //edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<reco::GenJetCollection> genJetToken_;
  edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> genJetFlavourToken_;

};

} // namespace

using namespace cat;

typedef std::vector<int> vint;
typedef std::vector<vint> vvint;

GenJetHadronMatchProducer::GenJetHadronMatchProducer(const edm::ParameterSet& pset)
{
  //genParticlesToken_ = consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("genParticle"));
  genJetToken_ = consumes<reco::GenJetCollection>(pset.getParameter<edm::InputTag>("genJet"));
  genJetFlavourToken_ = consumes<reco::JetFlavourInfoMatchingCollection>(pset.getParameter<edm::InputTag>("genJetFlavour"));

  produces<vint>("flavour");
  produces<vvint>("ancestorIds");
}

void GenJetHadronMatchProducer::produce(edm::Event& event, const edm::EventSetup&)
{
  edm::Handle<reco::GenJetCollection> genJetHandle;
  event.getByToken(genJetToken_, genJetHandle);

  edm::Handle<reco::JetFlavourInfoMatchingCollection> genJetFlavourHandle;
  event.getByToken(genJetFlavourToken_, genJetFlavourHandle);

  // Collect HF, LF jets
  std::map<int, int> jetToFlavour;
  std::map<int, std::vector<int> > jetToPartonIds;
  for ( auto& jetFlavPair : *genJetFlavourHandle ) {
    const auto& jetRefBase = jetFlavPair.first;
    const auto& flavInfo = jetFlavPair.second;
    const int jetIndex = jetRefBase.key();

    const auto& bHadrons = flavInfo.getbHadrons();
    const auto& cHadrons = flavInfo.getcHadrons();
    //const auto& partons = flavInfo.getPartons();

    if      ( !bHadrons.empty() ) jetToFlavour[jetIndex] = 5;
    else if ( !cHadrons.empty() ) jetToFlavour[jetIndex] = 4;
    else jetToFlavour[jetIndex] = 0;

    std::set<int> partonIds;
    for ( auto& x : jetRefBase->getJetConstituents() ) {
      if ( x.isNull() ) continue;
      collectAncestorPartons(x.get(), partonIds);
    }

    jetToPartonIds.insert(std::make_pair(jetIndex, std::vector<int>()));
    auto& vPartonIds = jetToPartonIds[jetIndex];
    vPartonIds.reserve(partonIds.size());
    std::copy(partonIds.begin(), partonIds.end(), std::back_inserter(vPartonIds));
  }

  std::auto_ptr<vint> out_flavours(new vint);
  std::auto_ptr<vvint> out_ancestorIds(new vvint);
  for ( int i=0, n=genJetHandle->size(); i<n; ++i ) {
    auto fitr = jetToFlavour.find(i);
    if ( fitr == jetToFlavour.end() ) out_flavours->push_back(0);
    else out_flavours->push_back(fitr->second);

    auto itr = jetToPartonIds.find(i);
    if ( itr == jetToPartonIds.end() ) out_ancestorIds->push_back(vint());
    else out_ancestorIds->push_back(itr->second);
  }
  event.put(out_flavours, "flavour");
  event.put(out_ancestorIds, "ancestorIds");

}

void GenJetHadronMatchProducer::collectAncestorPartons(const reco::Candidate* cand, std::set<int>& partonIds) const {
  if ( !cand ) return;

  for ( int i=0, n=cand->numberOfMothers(); i<n; ++i ) {
    const reco::Candidate* p = cand->mother(i);
    if ( !p or p->status() == 4 ) continue;
    const int aid = abs(p->pdgId());
    if ( aid == 6 or aid == 24 or aid == 23 or aid == 25 ) partonIds.insert(p->pdgId());

    collectAncestorPartons(p, partonIds);
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(GenJetHadronMatchProducer);
