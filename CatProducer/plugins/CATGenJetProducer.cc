#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CATTools/DataFormats/interface/GenJet.h"
#include "CATTools/DataFormats/interface/MCParticle.h"

using namespace edm;
using namespace std;

namespace cat {

  class CATGenJetProducer : public edm::stream::EDProducer<> {
  public:
    explicit CATGenJetProducer(const edm::ParameterSet & iConfig);

    void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

  private:
    edm::EDGetTokenT<reco::GenJetCollection> src_;
    const double pt_;
    const double eta_;

    std::vector<const reco::Candidate *> getAncestors(const reco::Candidate &c);
    bool hasBottom(const reco::Candidate &c);
    bool hasCharm(const reco::Candidate &c);
    bool decayFromBHadron(const reco::Candidate &c);
    bool decayFromCHadron(const reco::Candidate &c);
    const reco::Candidate* lastBHadron(const reco::Candidate &c);
    const reco::Candidate* lastCHadron(const reco::Candidate &c);

  };

} // namespace

cat::CATGenJetProducer::CATGenJetProducer(const edm::ParameterSet & iConfig) :
  src_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("src"))),
  pt_(iConfig.getParameter<double>("pt")),
  eta_(iConfig.getParameter<double>("eta"))
{
  produces<cat::GenJetCollection>();
}

void
cat::CATGenJetProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup)
{
  Handle<reco::GenJetCollection> src;
  iEvent.getByToken(src_, src);

  unique_ptr<cat::GenJetCollection>  out(new cat::GenJetCollection());

  for (const reco::GenJet & aGenJet : *src) {
    if ( aGenJet.pt() < pt_ || std::abs(aGenJet.eta()) > eta_ ) continue;

    cat::GenJet aCatGenJet(aGenJet);

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

    aCatGenJet.setHadron(matched);
    aCatGenJet.setPdgId(matched.pdgId());
    // int partonFlavour = aGenJet.partonFlavour();
    // int partonPdgId = aGenJet.genParton();
    // temp - find better way to match flavour and id
    aCatGenJet.setPartonFlavour(matched.pdgId());
    aCatGenJet.setPartonPdgId(matched.pdgId());

    out->push_back(aCatGenJet);

  }

  iEvent.put(std::move(out));
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

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATGenJetProducer);
