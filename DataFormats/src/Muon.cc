#include "CATTools/DataFormats/interface/Muon.h"

using namespace cat;

/// default constructor
Muon::Muon() {
}

Muon::Muon(const reco::LeafCandidate & aMuon) : Particle( aMuon ) {
}

/// destructor
Muon::~Muon() {
}
