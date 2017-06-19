#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"

using namespace edm;
using namespace std;

class GenJetHadronMatchProducer : public edm::stream::EDProducer<>
{
public:
  GenJetHadronMatchProducer(const edm::ParameterSet& pset);

  void produce(edm::Event & event, const edm::EventSetup&) override;

private:
  void collectAncestorPartons(const reco::Candidate* cand, std::set<int>& partonIds) const;

  //edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<reco::GenJetCollection> genJetToken_;
  edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> genJetFlavourToken_;

};

typedef std::vector<int> vint;
typedef std::vector<vint> vvint;

GenJetHadronMatchProducer::GenJetHadronMatchProducer(const edm::ParameterSet& pset)
{
  //genParticlesToken_ = consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("genParticle"));
  genJetToken_ = consumes<reco::GenJetCollection>(pset.getParameter<edm::InputTag>("genJet"));
  genJetFlavourToken_ = consumes<reco::JetFlavourInfoMatchingCollection>(pset.getParameter<edm::InputTag>("genJetFlavour"));

  produces<vint>("nBHadron");
  produces<vint>("nCHadron");
  produces<vvint>("ancestors");
  produces<vvint>("hadronAncestors");
}

void GenJetHadronMatchProducer::produce(edm::Event& event, const edm::EventSetup&)
{
  edm::Handle<reco::GenJetCollection> genJetHandle;
  event.getByToken(genJetToken_, genJetHandle);

  edm::Handle<reco::JetFlavourInfoMatchingCollection> genJetFlavourHandle;
  event.getByToken(genJetFlavourToken_, genJetFlavourHandle);

  // Collect HF, LF jets
  std::map<int, int> jetToNBHadron, jetToNCHadron;
  std::map<int, vint> jetToPartonIds;
  std::map<int, vint> jetToHPartonIds;
  for ( auto& jetFlavPair : *genJetFlavourHandle ) {
    const auto& jetRefBase = jetFlavPair.first;
    const auto& flavInfo = jetFlavPair.second;
    const int jetIndex = jetRefBase.key();

    const auto& bHadrons = flavInfo.getbHadrons();
    const auto& cHadrons = flavInfo.getcHadrons();
    //const auto& partons = flavInfo.getPartons();

    jetToNBHadron[jetIndex] = bHadrons.size();
    jetToNCHadron[jetIndex] = cHadrons.size();

    std::set<int> hPartonIds;
    for ( auto& x : bHadrons ) collectAncestorPartons(x.get(), hPartonIds);
    for ( auto& x : cHadrons ) collectAncestorPartons(x.get(), hPartonIds);

    jetToHPartonIds.insert(std::make_pair(jetIndex, vint()));
    auto& vHPartonIds = jetToHPartonIds[jetIndex];
    vHPartonIds.reserve(hPartonIds.size());
    std::copy(hPartonIds.begin(), hPartonIds.end(), std::back_inserter(vHPartonIds));

    std::set<int> partonIds;
    for ( auto& x : jetRefBase->getJetConstituents() ) {
      if ( x.isNull() ) continue;
      collectAncestorPartons(x.get(), partonIds);
    }

    jetToPartonIds.insert(std::make_pair(jetIndex, vint()));
    auto& vPartonIds = jetToPartonIds[jetIndex];
    vPartonIds.reserve(partonIds.size());
    std::copy(partonIds.begin(), partonIds.end(), std::back_inserter(vPartonIds));
  }

  std::unique_ptr<vint> out_nBHadron(new vint), out_nCHadron(new vint);
  std::unique_ptr<vvint> out_ancestors(new vvint);
  std::unique_ptr<vvint> out_hAncestors(new vvint);
  for ( int i=0, n=genJetHandle->size(); i<n; ++i ) {
    auto nbitr = jetToNBHadron.find(i);
    if ( nbitr == jetToNBHadron.end() ) out_nBHadron->push_back(0);
    else out_nBHadron->push_back(nbitr->second);

    auto ncitr = jetToNCHadron.find(i);
    if ( ncitr == jetToNCHadron.end() ) out_nCHadron->push_back(0);
    else out_nCHadron->push_back(ncitr->second);

    auto pitr = jetToPartonIds.find(i);
    if ( pitr == jetToPartonIds.end() ) out_ancestors->push_back(vint());
    else out_ancestors->push_back(pitr->second);

    auto hitr = jetToHPartonIds.find(i);
    if ( hitr == jetToHPartonIds.end() ) out_hAncestors->push_back(vint());
    else out_hAncestors->push_back(hitr->second);
  }
  event.put(std::move(out_nBHadron), "nBHadron");
  event.put(std::move(out_nCHadron), "nCHadron");
  event.put(std::move(out_ancestors), "ancestors");
  event.put(std::move(out_hAncestors), "hadronAncestors");

}

void GenJetHadronMatchProducer::collectAncestorPartons(const reco::Candidate* cand, std::set<int>& partonIds) const
{
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
DEFINE_FWK_MODULE(GenJetHadronMatchProducer);
