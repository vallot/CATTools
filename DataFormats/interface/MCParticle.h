#ifndef CATTools_MCParticle_H
#define CATTools_MCParticle_H

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// Define typedefs for convenience
namespace cat {
  class MCParticle;
  typedef std::vector<MCParticle>              MCParticleCollection;
  typedef edm::Ref<MCParticleCollection>       MCParticleRef;
  typedef edm::RefVector<MCParticleCollection> MCParticleRefVector;
}

namespace cat {

  class MCParticle : public reco::LeafCandidate{
  public:
    MCParticle();
    //    MCParticle(const reco::GenParticle & aMCParticle);
    MCParticle(const reco::Candidate & aMCParticle);
    virtual ~MCParticle();

  private:

  };
}

#endif
