#include "CATTools/DataFormats/interface/Jet.h"
#include <unordered_map>
#include <algorithm>

using namespace cat;

/// default constructor
Jet::Jet() {
}

Jet::Jet(const reco::LeafCandidate & aJet) : Particle( aJet ) {
}

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
void Jet::bDiscriminatorPrint() const {
  for(unsigned int i=0; i!=pairDiscriVector_.size(); i++){
    std::cout << pairDiscriVector_[i].first << " = " << pairDiscriVector_[i].second << std::endl;
  }
}

float Jet::smearedRes(int direction) const {
  const auto aGenJet = this->genJet();
  if ( !aGenJet ) return 1; // No JER

  const double absEta = std::abs(this->eta());

  // 2012 values from https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
  // need to update for run2
  const int nbins = 7;
  const double etaBins[nbins] = {0.5, 1.1, 1.7, 2.3, 2.8, 3.2, 5.0};
  const double cJERs[nbins]   = {1.079, 1.099, 1.121, 1.208, 1.254, 1.395, 1.056};
  const double cJERsUp[nbins] = {1.105, 1.127, 1.150, 1.254, 1.316, 1.458, 1.247};
  const double cJERsDn[nbins] = {1.053, 1.071, 1.092, 1.162, 1.192, 1.332, 0.865};
  // call lower_bound to find bin location.
  // lower_bound returns 0 if absEta < etaBins[0]
  // lower_bound returns nbins if absEta > etaBins[nbins-1]. Extrapolate JER factor for higher eta
  const size_t bin = std::min(int(std::lower_bound(etaBins, etaBins+nbins, absEta)-etaBins), nbins);

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
