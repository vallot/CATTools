#include "CATTools/DataFormats/interface/Jet.h"

using namespace cat;

/// default constructor
Jet::Jet() {
}

Jet::Jet(const reco::LeafCandidate & aJet) : Particle( aJet ) {
}

/// destructor
Jet::~Jet() {
}
