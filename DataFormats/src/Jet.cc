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

float Jet::jecFactor(const std::string & aLabel) const {
  float jec = 0.;
  for(unsigned int i=0; i!=jecFactor_.size(); i++){
    if(jecFactor_[i].first == aLabel){
      jec = jecFactor_[i].second;
    }
  }
  return jec;
}

Jet Jet::correctedJet(const std::string &theLabel) const {
  float jec = jecFactor(theLabel);
  Jet correctedJet(*this);
  //rescale p4 of the jet
  correctedJet.setP4(jec*p4());
  return correctedJet;
}
