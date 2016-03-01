#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CATTools/DataFormats/interface/GenJet.h"
#include "CATTools/DataFormats/interface/MCParticle.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"

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
  int getFlavour(const int pdgId);

private:
  edm::EDGetTokenT<reco::GenJetCollection> genJetToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;

  std::shared_ptr<fastjet::JetDefinition> fjDef_;

};

} // namespace

using namespace cat;

GenJetHadronMatchProducer::GenJetHadronMatchProducer(const edm::ParameterSet& pset)
{
  genJetToken_ = consumes<reco::GenJetCollection>(pset.getParameter<edm::InputTag>("genJet"));
  genParticlesToken_ = consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("genParticle"));

  fjDef_ = std::shared_ptr<fastjet::JetDefinition>(new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4));

  produces<edm::ValueMap<int> >();
}

void GenJetHadronMatchProducer::produce(edm::Event& event, const edm::EventSetup&)
{
  // Collect heavy-flavour hadrons
  edm::Handle<reco::GenParticleCollection> genParticlesHandle;
  event.getByToken(genParticlesToken_, genParticlesHandle);
  std::set<size_t> bHadrons, cHadrons;
  for ( size_t i=0, n=genParticlesHandle->size(); i<n; ++i ) {
    const auto& x = genParticlesHandle->at(i);
    if ( x.status() == 1 or x.status() == 4 ) continue; // Skip final states and incident beams

    const int aid = abs(x.pdgId());
    if ( aid < 100 ) continue; // Skip partons

    const int nDau = x.numberOfDaughters();
    if ( nDau == 0 ) continue; // Skip final states, just for confirmation
    if ( !x.isLastCopy() ) continue; // Consider the last copy only

    const int flav = getFlavour(aid);
    if ( flav < 4 ) continue;

    bool isLast = true;
    for ( int j=0; j<nDau; ++j ) {
      const int dauFlav = getFlavour(x.daughter(j)->pdgId());
      if ( flav == dauFlav ) { isLast = false; break; }
    }
    if ( !isLast ) continue; // Keep last hadrons only

    if      ( flav == 5 ) bHadrons.insert(i);
    else if ( flav == 4 ) cHadrons.insert(i);
  }
  const size_t nBHadrons = bHadrons.size();
  const size_t nCHadrons = cHadrons.size();

  edm::Handle<reco::GenJetCollection> genJetHandle;
  event.getByToken(genJetToken_, genJetHandle);

  std::vector<fastjet::PseudoJet> fjInputs;
  fjInputs.reserve(nBHadrons+nCHadrons+genJetHandle->size()*100); // nJetx100 from very crude estimation for the size of all jet constituents
  for ( const reco::GenJet & aGenJet : *genJetHandle ) {
    for ( auto& p : aGenJet.getJetConstituents() ) {
      if ( p.isNull() ) continue;
      fjInputs.push_back(fastjet::PseudoJet(p->px(), p->py(), p->pz(), p->energy()));
      fjInputs.back().set_user_index(fjInputs.size()); // User index for usual particles. This user index is forced to set >=1, to avoid overlap with B hadrons at index=0
    }
  }
  // Do the reclustering including the hadrons
  for ( auto& ip : bHadrons ) {
    const auto& p = genParticlesHandle->at(ip);
    const double p0 = p.p4().P();
    const double fP = 1E-7/p0; // Rescale to negligible factor, hadron will have 1e-7 GeV in momentum
    const double e = std::hypot(p.mass(), p0*fP); // re-calculate energy to reserve particle mass
    fjInputs.push_back(fastjet::PseudoJet(p.px()*fP, p.py()*fP, p.pz()*fP, e));
    fjInputs.back().set_user_index(-ip);
  }
  for ( auto& ip : cHadrons ) {
    const auto& p = genParticlesHandle->at(ip);
    const double p0 = p.p4().P();
    const double fP = 1E-7/p0; // Rescale to negligible factor, hadron will have 1e-7 GeV in momentum
    const double e = std::hypot(p.mass(), p0*fP); // re-calculate energy to reserve particle mass
    fjInputs.push_back(fastjet::PseudoJet(p.px()*fP, p.py()*fP, p.pz()*fP, e));
    fjInputs.back().set_user_index(-ip);
  }

  // Run the FastJet
  fastjet::ClusterSequence fjClusterSeq(fjInputs, *fjDef_);
  std::vector<fastjet::PseudoJet> fjJets = fastjet::sorted_by_pt(fjClusterSeq.inclusive_jets(5.0));

  // Collect heavy flavour jets
  std::map<int, int> jetToFlav;
  //std::vector<std::vector<int> > jetsToHadrons;
  for ( size_t i=0, n=fjJets.size(); i<n; ++i ) {
    const auto& jet = fjJets.at(i);
    const auto& fjCons = jet.constituents();
    //std::vector<int> matchedBHadrons, matchedCHadrons;

    int flav = 0;
    for ( auto& con : fjCons ) {
      const int index = con.user_index();
      if ( index >= 1 ) continue; // We are not intestested in the constituents from original jet constituents

      //if      ( bHadrons.find(index) != bHadrons.end() ) matchedBHadrons.push_back(-index);
      //else if ( cHadrons.find(index) != cHadrons.end() ) matchedCHadrons.push_back(-index);
      if      ( bHadrons.find(-index) != bHadrons.end() ) flav = 5;
      else if ( cHadrons.find(-index) != cHadrons.end() ) flav = 4;
    }

    jetToFlav[i] = flav;
    //std::sort(matchedBHadrons.begin(), matchedBHadrons.end(), [&](int a, int b){return genParticlesHandle->at(a).pt() > genParticlesHandle->at(b).pt();});
    //std::sort(matchedCHadrons.begin(), matchedCHadrons.end(), [&](int a, int b){return genParticlesHandle->at(a).pt() > genParticlesHandle->at(b).pt();});
    //jetsToHadrons.push_back(matchedBHadrons);
    //jetsToHadrons.back().insert(jetsToHadrons.back().end(), matchedCHadrons.begin(), matchedCHadrons.end());
  }

  // Start main loop to prepare GenJets
  std::auto_ptr<edm::ValueMap<int> > out(new edm::ValueMap<int>);
  std::vector<int> flavours;
  flavours.reserve(genJetHandle->size());
  for ( size_t i=0, n=genJetHandle->size(); i<n; ++i ) {
    //const reco::GenJet& genJet = genJetHandle->at(i);
    // Assume that the genJets indices are unchanged by the reclustering
    int flav = 0;
    if ( jetToFlav.find(i) != jetToFlav.end() ) flav = jetToFlav[i];

    flavours.push_back(flav);
  }
  edm::ValueMap<int>::Filler outFiller(*out);
  outFiller.insert(genJetHandle, flavours.begin(), flavours.end());
  outFiller.fill();
  event.put(out);

}

int GenJetHadronMatchProducer::getFlavour(const int pdgId) {
  const int aid = abs(pdgId);
  const int code1 = (aid/ 100)%10;
  const int code2 = (aid/1000)%10;

  return std::max(code1, code2);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(GenJetHadronMatchProducer);
