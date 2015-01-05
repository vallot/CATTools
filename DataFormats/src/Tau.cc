#include "CATTools/DataFormats/interface/Tau.h"

using namespace cat;

/// default constructor
Tau::Tau() {
}

Tau::Tau(const reco::LeafCandidate & aTau) : Particle( aTau ) {
}

/// destructor
Tau::~Tau() {
}
