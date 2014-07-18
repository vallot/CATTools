#include "CATTools/DataFormats/interface/MET.h"

using namespace cat;

/// default constructor
MET::MET() {
}

MET::MET(const reco::LeafCandidate & aMET) : Particle( aMET ) {
}

/// destructor
MET::~MET() {
}
