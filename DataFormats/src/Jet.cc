#include "CATTools/DataFormats/interface/Jet.h"
#include <unordered_map>
#include <algorithm>

using namespace cat;

/// default constructor
Jet::Jet():
  fJER_(1), fJERUp_(1), fJERDown_(1) {
}

Jet::Jet(const reco::LeafCandidate & aJet) :
  Particle( aJet ),
  tightJetID_(false),
  tightLepVetoJetID_(false),
  pileupJetId_(0),
  chargedEmEnergyFraction_(0),
  vtxMass_(0),
  vtxNtracks_(0),
  vtx3DVal_(0),
  vtx3DSig_(0),
  partonFlavour_(0),
  hadronFlavour_(0),
  partonPdgId_(0),
  shiftedEnDown_(1),
  shiftedEnUp_(1),					     
  fJER_(1), fJERUp_(1), fJERDown_(1),
  qgLikelihood_(-2)
{}

/// destructor
Jet::~Jet() {
}

/// get b discriminant from label name
float Jet::bDiscriminator(const std::string & aLabel) const {
  float discriminator = -1000.;
  const std::string & theLabel = ((aLabel == "" || aLabel == "default")) ? "trackCountingHighEffBJetTags" : aLabel;
  for(unsigned int i=0; i!=pairDiscriVector_.size(); i++){
    if(pairDiscriVector_[i].first == theLabel){
      discriminator = pairDiscriVector_[i].second;
    }
  }
  return discriminator;
}

/// print all bjet Discriminators
int Jet::printBDiscriminator() const {
  for(unsigned int i=0; i!=pairDiscriVector_.size(); i++){
    std::cout << pairDiscriVector_[i].first << " = " << pairDiscriVector_[i].second << std::endl;
  }
  return 1;
}

float Jet::smearedRes(int direction, int era) const {
  // The era-based JER is going to be removed
  if ( era == 0 ) return direction == 0 ? fJER_ : direction > 0 ? fJERUp_ : fJERDown_;

  const auto aGenJet = this->genJet();
  if ( !aGenJet ) return 1; // No JER

  const double absEta = std::abs(this->eta());
  if (absEta >=5.0 ) return 1; // No JER
    
  //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
  std::vector<double> etaBins = {5.0};
  std::vector<double> cJERs = {1.0};
  std::vector<double> cJERsUp = {1.0};
  std::vector<double> cJERsDn = {1.0};
  if ( era == 2018 ) { //This need to be updated
    etaBins = {0.5, 0.8, 1.1, 1.3, 1.7, 1.9, 2.1, 2.3, 2.5, 2.8, 3.0, 3.2, 5.0};
    cJERs   = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
    cJERsUp = {1.109+0.008, 1.138+0.013, 1.114+0.013, 1.123+0.024, 1.084+0.011, 1.082+0.035, 1.140+0.047, 1.067+0.053, 1.177+0.041, 1.364+0.039, 1.857+0.071, 1.328+0.022, 1.16+0.029};
    cJERsDn = {1.109-0.008, 1.138-0.013, 1.114-0.013, 1.123-0.024, 1.084-0.011, 1.082-0.035, 1.140-0.047, 1.067-0.053, 1.177-0.041, 1.364-0.039, 1.857-0.071, 1.328-0.022, 1.16-0.029};
  }
  // call lower_bound to find bin location.
  const size_t bin = std::lower_bound(etaBins.begin(), etaBins.end(), absEta) - etaBins.begin();
  const double jetPt = this->pt();
  const double genJetPt = aGenJet->pt();
  const double dPt = jetPt-genJetPt;

  double cJER = 0;
  if      ( direction == 0 ) cJER = cJERs[bin];
  else if ( direction >  0 ) cJER = cJERsUp[bin];
  else  cJER = cJERsDn[bin];
  
  const double fJER = std::max(0., (genJetPt+dPt*cJER)/jetPt);
  return fJER;
}
