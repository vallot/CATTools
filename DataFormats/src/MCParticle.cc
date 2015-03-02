#include "CATTools/DataFormats/interface/MCParticle.h"

using namespace cat;

/// default constructor
MCParticle::MCParticle(){
}

// MCParticle::MCParticle(const reco::GenParticle & aMCParticle) : reco::LeafCandidate(aMCParticle) {
// }

MCParticle::MCParticle(const reco::Candidate & aMCParticle) : reco::LeafCandidate(aMCParticle) {
}

/// destructor
MCParticle::~MCParticle() {
}

