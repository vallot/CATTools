#include "CATTools/DataFormats/interface/Jet.h"
#include <unordered_map>
#include <algorithm>

using namespace cat;

/// default constructor
Jet::Jet() {
}

Jet::Jet(const reco::LeafCandidate & aJet) : Particle( aJet ) {
}

/// destructor
Jet::~Jet() {
}
