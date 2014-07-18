#include "CATTools/DataFormats/interface/Particle.h"

using namespace cat;

/// default constructor
Particle::Particle(){
}

Particle::Particle(const reco::LeafCandidate & aParticle) : reco::LeafCandidate(aParticle) {
}

/// destructor
Particle::~Particle() {
}

