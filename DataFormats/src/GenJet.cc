#include "CATTools/DataFormats/interface/GenJet.h"
#include <unordered_map>
#include <algorithm>

using namespace cat;

/// default constructor
GenJet::GenJet() {
}

GenJet::GenJet(const reco::GenJet & aGenJet) : Particle( aGenJet ) {
}

/// destructor
GenJet::~GenJet() {
}
