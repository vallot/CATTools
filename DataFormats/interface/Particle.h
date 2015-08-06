#ifndef CATTools_Particle_H
#define CATTools_Particle_H 

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "TLorentzVector.h"
// Define typedefs for convenience
namespace cat {
  class Particle;
  typedef std::vector<Particle>              ParticleCollection;
  typedef edm::Ref<ParticleCollection>       ParticleRef;
  typedef edm::RefVector<ParticleCollection> ParticleRefVector;
}

namespace cat {

  class Particle : public reco::LeafCandidate{
  public:
    Particle();
    Particle(const reco::LeafCandidate & aParticle); 
    virtual ~Particle();
    TLorentzVector tlv() const {return TLorentzVector(this->px(), this->py(),this->pz(),this->energy());} 
  private:

  };
}

#endif
