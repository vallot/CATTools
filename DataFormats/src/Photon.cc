#include "CATTools/DataFormats/interface/Photon.h"

using namespace cat;

/// default constructor
Photon::Photon() {
}

Photon::Photon(const reco::LeafCandidate & aPhoton) : Particle( aPhoton ) {
}

/// destructor
Photon::~Photon() {
}
