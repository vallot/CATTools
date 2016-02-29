#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CATTools/DataFormats/interface/GenJet.h"
#include "CATTools/DataFormats/interface/MCParticle.h"

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"

using namespace edm;
using namespace std;

namespace cat {

  class CATGenJetProducer : public edm::stream::EDProducer<> {
  public:
    explicit CATGenJetProducer(const edm::ParameterSet & iConfig);
    virtual ~CATGenJetProducer() { }

    void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

  private:
    edm::EDGetTokenT<reco::GenJetCollection> src_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
    const double pt_;
    const double eta_;
    std::shared_ptr<fastjet::JetDefinition> fjDef_;

    std::vector<const reco::Candidate *> getAncestors(const reco::Candidate &c);
    bool hasBottom(const reco::Candidate &c);
    bool hasCharm(const reco::Candidate &c);
    bool decayFromBHadron(const reco::Candidate &c);
    bool decayFromCHadron(const reco::Candidate &c);
    const reco::Candidate* lastBHadron(const reco::Candidate &c);
    const reco::Candidate* lastCHadron(const reco::Candidate &c);
    int getFlavour(const int pdgId);

  };

} // namespace

cat::CATGenJetProducer::CATGenJetProducer(const edm::ParameterSet & iConfig) :
  src_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("src"))),
  genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  pt_(iConfig.getParameter<double>("pt")),
  eta_(iConfig.getParameter<double>("eta")),
  fjDef_(std::shared_ptr<fastjet::JetDefinition>(new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4)))
{
  produces<std::vector<cat::GenJet> >();
}

void cat::CATGenJetProducer::produce(edm::Event & iEvent, const edm::EventSetup&)
{
  auto_ptr<vector<cat::GenJet> >  out(new vector<cat::GenJet>());

  // Collect heavy-flavour hadrons
  edm::Handle<reco::GenParticleCollection> genParticlesHandle;
  iEvent.getByToken(genParticlesToken_, genParticlesHandle);
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
    for ( int i=0; i<nDau; ++i ) {
      const int dauFlav = getFlavour(x.daughter(i)->pdgId());
      if ( flav == dauFlav ) { isLast = false; break; }
    }
    if ( !isLast ) continue; // Keep last hadrons only

    if      ( flav == 5 ) bHadrons.insert(i);
    else if ( flav == 4 ) cHadrons.insert(i);
  }
  const size_t nBHadrons = bHadrons.size();
  const size_t nCHadrons = cHadrons.size();

  Handle<reco::GenJetCollection> src;
  iEvent.getByToken(src_, src);

  std::vector<fastjet::PseudoJet> fjInputs;
  fjInputs.reserve(nBHadrons+nCHadrons+src->size()*100); // nJetx100 from very crude estimation for the size of all jet constituents
  for ( const reco::GenJet & aGenJet : *src ) {
    for ( auto& p : aGenJet.getJetConstituents() ) {
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
  std::vector<fastjet::PseudoJet> fjJets = fastjet::sorted_by_pt(fjClusterSeq.inclusive_jets(pt_));

  // Collect heavy flavour jets
  std::vector<std::vector<int> > jetsToHadrons;
  for ( auto& jet : fjJets ) {
    if ( std::abs(jet.eta()) > eta_ ) continue;
    const auto& fjCons = jet.constituents();
    std::vector<int> matchedBHadrons, matchedCHadrons;
    for ( auto& con : fjCons ) {
      const int index = con.user_index();
      if ( index >= 1 ) continue; // We are not intestested in the constituents from original jet constituents

      if ( bHadrons.find(index) != bHadrons.end() ) matchedBHadrons.push_back(-index);
      else if ( cHadrons.find(index) != cHadrons.end() ) matchedCHadrons.push_back(-index);
    }
    std::sort(matchedBHadrons.begin(), matchedBHadrons.end(), [&](int a, int b){return genParticlesHandle->at(a).pt() > genParticlesHandle->at(b).pt();});
    std::sort(matchedCHadrons.begin(), matchedCHadrons.end(), [&](int a, int b){return genParticlesHandle->at(a).pt() > genParticlesHandle->at(b).pt();});
    jetsToHadrons.push_back(matchedBHadrons);
    jetsToHadrons.back().insert(jetsToHadrons.back().end(), matchedCHadrons.begin(), matchedCHadrons.end());
  }

  // Start main loop to prepare cat::GenJets
  for ( size_t i=0, n=src->size(); i<n; ++i ) {
    const reco::GenJet& aGenJet = src->at(i);
    if ( aGenJet.pt() < pt_ || std::abs(aGenJet.eta()) > eta_ ) continue;

    cat::GenJet aCatGenJet(aGenJet);

    // Hadron matching based on particle decay history
    cat::MCParticle matched;
    reco::Jet::Constituents jc = aGenJet.getJetConstituents();
    //if B-Hadron matched, always assign B-Hadron
    for ( reco::Jet::Constituents::const_iterator itr = jc.begin(); itr != jc.end(); ++itr ){
      if (itr->isAvailable()){
        const reco::Candidate* mcpart = dynamic_cast<const reco::Candidate*>(itr->get());
        const reco::Candidate* lastB = lastBHadron(*mcpart);
        if (lastB){
          matched = cat::MCParticle(*lastB);
          break;
        }
      }
    }
    if (std::abs(matched.pdgId()) != 5){
      //if only no B-Hadron matched, assign C-Hadron
      for ( reco::Jet::Constituents::const_iterator itr = jc.begin(); itr != jc.end(); ++itr ){
        if (itr->isAvailable()){
          const reco::Candidate* mcpart = dynamic_cast<const reco::Candidate*>(itr->get());
          const reco::Candidate* lastC = lastCHadron(*mcpart);
          if (lastC){
            matched = cat::MCParticle(*lastC);
            break;
          }
        }
      }
    }

    // Another hadron matching based on ghost tagging
    // Assume that the genJets indices are unchanged by the reclustering
    cat::MCParticle matchedGhost;
    const auto& matchedHadrons = jetsToHadrons.at(i);
    if ( !matchedHadrons.empty() ) matchedGhost = cat::MCParticle(genParticlesHandle->at(matchedHadrons.at(0)));

    aCatGenJet.setHadron(matched);
    aCatGenJet.setPdgId(matched.pdgId());
    // int partonFlavour = aGenJet.partonFlavour();
    // int partonPdgId = aGenJet.genParton();
    // temp - find better way to match flavour and id
    aCatGenJet.setPartonFlavour(matched.pdgId());
    aCatGenJet.setPartonPdgId(matched.pdgId());

    aCatGenJet.setGhost(matchedGhost);

    out->push_back(aCatGenJet);

  }

  iEvent.put(out);
}

std::vector<const reco::Candidate *> cat::CATGenJetProducer::getAncestors(const reco::Candidate &c)
{
  vector<const reco::Candidate *> moms;
  if( c.numberOfMothers() == 1 ) {
    const reco::Candidate * dau = &c;
    const reco::Candidate * mom = c.mother();
    while ( dau->numberOfMothers() == 1) {
      moms.push_back( dau );
      dau = mom ;
      mom = dau->mother();
    }
  }
  return moms;
}

bool cat::CATGenJetProducer::hasBottom(const reco::Candidate &c)
{
  int code1;
  int code2;
  bool tmpHasBottom = false;
  code1 = (int)( ( abs(c.pdgId() ) / 100)%10 );
  code2 = (int)( ( abs(c.pdgId() ) /1000)%10 );
  if ( code1 == 5 || code2 == 5) tmpHasBottom = true;
  return tmpHasBottom;
}

bool cat::CATGenJetProducer::hasCharm(const reco::Candidate &c)
{
  int code1;
  int code2;
  bool tmpHasCharm = false;
  code1 = (int)( ( abs(c.pdgId() ) / 100)%10 );
  code2 = (int)( ( abs(c.pdgId() ) /1000)%10 );
  if ( code1 == 4 || code2 == 4) tmpHasCharm = true;
  return tmpHasCharm;
}

bool cat::CATGenJetProducer::decayFromBHadron(const reco::Candidate & c)
{
  bool isFromB = false;
  vector<const reco::Candidate *> allParents = getAncestors( c );
  for( vector<const reco::Candidate *>::const_iterator aParent = allParents.begin();
       aParent != allParents.end();
       aParent ++ )
    {
      if( hasBottom(**aParent) ) isFromB = true;
      /*
	cout << " particle Parent is " << (*aParent)->status()
	<< " type " << (*aParent)->pdgId()
	<< " pt= " << (*aParent)->pt()
	<< " isB = " << isFromB
	<< endl;
      */
    }
  return isFromB;
}

bool cat::CATGenJetProducer::decayFromCHadron(const reco::Candidate & c)
{
  bool isFromC = false;
  vector<const reco::Candidate *> allParents = getAncestors( c );
  for( vector<const reco::Candidate *>::const_iterator aParent = allParents.begin();
       aParent != allParents.end();
       aParent ++ )
    {
      if( hasCharm(**aParent) ) isFromC = true;
      /*
	cout << " particle Parent is " << (*aParent)->status()
	<< " type " << (*aParent)->pdgId()
	<< " pt=" << (*aParent)->pt()
	<< " isC = " << isFromC
	<< endl;
      */
    }
  return isFromC;
}


const reco::Candidate* cat::CATGenJetProducer::lastBHadron(const reco::Candidate & c)
{
  const reco::Candidate * out = 0;
  vector<const reco::Candidate *> allParents = getAncestors( c );
  for( vector<const reco::Candidate *>::const_iterator aParent = allParents.begin();
       aParent != allParents.end();
       aParent ++ )
    {
      if( hasBottom(**aParent) ) out = *aParent;
    }
  return out;
}

const reco::Candidate* cat::CATGenJetProducer::lastCHadron(const reco::Candidate & c)
{
  const reco::Candidate * out = 0;
  vector<const reco::Candidate *> allParents = getAncestors( c );
  for( vector<const reco::Candidate *>::const_iterator aParent = allParents.begin();
       aParent != allParents.end();
       aParent ++ )
    {
      if( hasCharm(**aParent) ) out = *aParent;
    }

  return out;
}

int cat::CATGenJetProducer::getFlavour(const int pdgId) {
  const int aid = abs(pdgId);
  const int code1 = (aid/ 100)%10;
  const int code2 = (aid/1000)%10;

  return std::max(code1, code2);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATGenJetProducer);
