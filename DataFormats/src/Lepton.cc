#include "CATTools/DataFormats/interface/Lepton.h"

using namespace cat;

/// default constructor
Lepton::Lepton() {
}

Lepton::Lepton(const reco::LeafCandidate & aLepton) : Particle( aLepton ) {
}

/// destructor
Lepton::~Lepton() {
}
