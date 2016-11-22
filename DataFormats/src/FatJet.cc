#include "CATTools/DataFormats/interface/FatJet.h"
#include <unordered_map>
#include <algorithm>

using namespace cat;

/// default constructor
FatJet::FatJet():
    Jet(),
    tau1_(-1), tau2_(-1), tau3_(-1),
    prunedMass_(0), softdropMass_(0),
    puppi_pt_(0), 
    puppi_eta_(0), 
    puppi_phi_(0), 
    puppi_mass_(0),
    puppi_tau1_(-1), puppi_tau2_(-1), puppi_tau3_(-1)
{
}

FatJet::FatJet(const reco::LeafCandidate & aJet) :
  Jet( aJet ),
    tau1_(-1), tau2_(-1), tau3_(-1),
    prunedMass_(0), softdropMass_(0),
    puppi_pt_(0), 
    puppi_eta_(0), 
    puppi_phi_(0), 
    puppi_mass_(0),
    puppi_tau1_(-1), puppi_tau2_(-1), puppi_tau3_(-1)
{}

/// destructor
FatJet::~FatJet() {
}
