#include "CATTools/DataFormats/interface/Electron.h"

using namespace cat;

/// default constructor
Electron::Electron() {
}

Electron::Electron(const reco::LeafCandidate & aElectron) : Particle( aElectron ) {
}

/// destructor
Electron::~Electron() {
}
