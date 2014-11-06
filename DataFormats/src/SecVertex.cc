#include "CATTools/DataFormats/interface/SecVertex.h"

using namespace cat;

/// default constructor
SecVertex::SecVertex() {
}

SecVertex::SecVertex(reco::VertexCompositeCandidate & aSecVertex) : reco::VertexCompositeCandidate( aSecVertex ) {
}

/// destructor
SecVertex::~SecVertex() {
}
