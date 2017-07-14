#include "CATTools/CommonTools/interface/GenParticleHelper.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

using namespace std;

std::vector<reco::GenParticleRef> cat::selectGenParticles(const edm::Handle<reco::GenParticleCollection>& srcColl,
                                                          const unsigned int pdgId, const std::vector<unsigned int>& motherIds)
{
  std::vector<reco::GenParticleRef> outColl;
  for ( size_t i=0, n=srcColl->size(); i<n; ++i ) {
    const auto& p = srcColl->at(i);
    if ( std::abs(p.pdgId()) != pdgId ) continue;

    for ( auto mother = p.mother(); mother!=0; mother = mother->mother() ) {
      if ( std::find(motherIds.begin(), motherIds.end(), mother->pdgId()) != motherIds.end() ) {
        reco::GenParticleRef pRef(srcColl, i);
        outColl.push_back(pRef);
        break;
      }
    }
  }
  return outColl;
}

bool cat::isMatchedByDRDPt(const reco::Candidate::LorentzVector& a, const reco::Candidate::LorentzVector& b, bool exact=false)
{
  const double dR = deltaR(a, b);
  const double dPt = a.pt() > 0 ? std::abs(a.pt()-b.pt())/a.pt() : 999;

  // If we are comparing two objects for which the candidates should
  // be exactly the same, cut hard. Otherwise take cuts from user.
  if ( exact ) return (dR < 1e-3 and dPt < 1e-3);
  return (dR < 0.025 and dPt < 0.025);
}

bool cat::isLast(const reco::GenParticle& p)
{
  const unsigned int nDau = p.numberOfDaughters();
  if ( nDau == 0 ) return true;

  const int id = p.pdgId();

  for ( size_t i=0; i<nDau; ++i ) {
    const reco::Candidate* dau = p.daughter(i);
    if ( !dau ) continue;
    if( dau->pdgId() == id) return false;
  }

  return true;
}
