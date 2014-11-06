#include "CATTools/DataFormats/interface/SecVertex.h"

using namespace cat;

/// default constructor
SecVertex::SecVertex() {
}

SecVertex::SecVertex(const reco::LeafCandidate & aSecVertex) : Particle( aSecVertex ) {
}

/// destructor
SecVertex::~SecVertex() {
}
