#include "CATTools/DataFormats/interface/Lepton.h"

using namespace cat;

/// default constructor
Lepton::Lepton() {
}

Lepton::Lepton(const reco::LeafCandidate & aLepton) :
  Particle( aLepton ),
  isPF_(false),
  isTight_(false),
  isMedium_(false),
  isLoose_(false),
  dz_(0),
  dxy_(0),
  chargedHadronIso03_(0),
  puChargedHadronIso03_(0),
  neutralHadronIso03_(0),
  photonIso03_(0),
  chargedHadronIso04_(0),
  puChargedHadronIso04_(0),
  neutralHadronIso04_(0),
  photonIso04_(0),
  mcMatched_(false)
{}

/// destructor
Lepton::~Lepton() {
}
