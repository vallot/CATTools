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
  typedef math::XYZPoint Point;
}

namespace cat {

  class Particle : public reco::LeafCandidate{
  public:
    Particle();
    Particle(const reco::LeafCandidate & aParticle);
    virtual ~Particle();

    const reco::GenParticle * genParticle() const {return genParticleRef_.get();}
    void setGenParticleRef(reco::GenParticleRef gj){ genParticleRef_ = gj;}
    bool hasGenParticle() const { return genParticleRef_.isNonnull(); }

    TLorentzVector tlv() const {return TLorentzVector(this->px(), this->py(),this->pz(),this->energy());}

  private:

    reco::GenParticleRef  genParticleRef_;

  };
}

#endif
