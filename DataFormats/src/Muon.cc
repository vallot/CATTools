#include "CATTools/DataFormats/interface/Muon.h"

using namespace cat;

/// default constructor
Muon::Muon() {
}

Muon::Muon(const reco::LeafCandidate & aMuon) : Lepton( aMuon ) {
}

/// destructor
Muon::~Muon() {
}
