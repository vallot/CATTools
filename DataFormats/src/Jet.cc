#include "CATTools/DataFormats/interface/Jet.h"
#include <unordered_map>
#include <algorithm>

using namespace cat;

namespace {
  std::unordered_map<std::string,cat::Jet::BTagWP> _btag_wp_map;
  const std::unordered_map<std::string,cat::Jet::BTagWP> & fill_btag_wp_map() {
    if (_btag_wp_map.empty()) {
      _btag_wp_map["TCHPT"] = cat::Jet::TCHPT;
      _btag_wp_map["JPL"] = cat::Jet::JPL;
      _btag_wp_map["JPM"] = cat::Jet::JPM;
      _btag_wp_map["JPT"] = cat::Jet::JPT;
      _btag_wp_map["CSVL"] = cat::Jet::CSVL;
      _btag_wp_map["CSVM"] = cat::Jet::CSVM;
      _btag_wp_map["CSVT"] = cat::Jet::CSVT;
      _btag_wp_map["CSVV1L"] = cat::Jet::CSVV1L;
      _btag_wp_map["CSVV1M"] = cat::Jet::CSVV1M;
      _btag_wp_map["CSVV1T"] = cat::Jet::CSVV1T;
      _btag_wp_map["CSVSLV1L"] = cat::Jet::CSVSLV1L;
      _btag_wp_map["CSVSLV1M"] = cat::Jet::CSVSLV1M;
      _btag_wp_map["CSVSLV1T"] = cat::Jet::CSVSLV1T;
      _btag_wp_map["CSVIVFV2L"] = cat::Jet::CSVIVFV2L;
      _btag_wp_map["CSVIVFV2M"] = cat::Jet::CSVIVFV2M;
      _btag_wp_map["CSVIVFV2T"] = cat::Jet::CSVIVFV2T;
    }
    return _btag_wp_map;
  }
}

/// default constructor
Jet::Jet() {
}

Jet::Jet(const reco::LeafCandidate & aJet) : Particle( aJet ) {
}

/// destructor
Jet::~Jet() {
}


double cat::Jet::btag(const char* s) const{
  if(!btag_.size()) return -9.9;
  const std::string name(s);
  const TagNameArray::const_iterator index = std::find( btagNames_.begin(),btagNames_.end(), name);
  if( index == btagNames_.end() ) return -9.9;
  return btag_.at( index - btagNames_.begin() );
}

bool cat::Jet::btagWP(BTagWP wp) const {
  switch(wp) {
    case TCHPT: return bDiscriminator("trackCountingHighPurBJetTags") > 3.41;
    case JPL: return bDiscriminator("jetProbabilityBJetTags") > 0.275;
    case JPM: return bDiscriminator("jetProbabilityBJetTags") > 0.545;
    case JPT: return bDiscriminator("jetProbabilityBJetTags") > 0.790;
    case CSVL: return bDiscriminator("combinedSecondaryVertexBJetTags") > 0.244;
    case CSVM: return bDiscriminator("combinedSecondaryVertexBJetTags") > 0.679;
    case CSVT: return bDiscriminator("combinedSecondaryVertexBJetTags") > 0.898;
    case CSVV1L: return bDiscriminator("combinedSecondaryVertexV1BJetTags") > 0.405;
    case CSVV1M: return bDiscriminator("combinedSecondaryVertexV1BJetTags") > 0.783;
    case CSVV1T: return bDiscriminator("combinedSecondaryVertexV1BJetTags") > 0.920;
    case CSVSLV1L: return bDiscriminator("combinedSecondaryVertexSoftPFLeptonV1BJetTags") > 0.527;
    case CSVSLV1M: return bDiscriminator("combinedSecondaryVertexSoftPFLeptonV1BJetTags") > 0.756;
    case CSVSLV1T: return bDiscriminator("combinedSecondaryVertexSoftPFLeptonV1BJetTags") > 0.859;
    case CSVIVFV2L: return bDiscriminator("combinedSecondaryVertexIVFV2BJetTags") > 0.423;
    case CSVIVFV2M: return bDiscriminator("combinedSecondaryVertexIVFV2BJetTags") > 0.814;
    case CSVIVFV2T: return bDiscriminator("combinedSecondaryVertexIVFV2BJetTags") > 0.941;
  }
  throw cms::Exception("InvalidArgument", "Unrecognized BTag WP choice");
}

bool cat::Jet::btagWP(const std::string &wp) const {
  const std::unordered_map<std::string,cat::Jet::BTagWP> & wp_map = fill_btag_wp_map();
  auto match = wp_map.find(wp);
  if (match == wp_map.end()) throw cms::Exception("InvalidArgument", wp+" btag WP not known");
  return btagWP(match->second);
}

bool cat::Jet::btagWP(const char *wp) const
{
  return btagWP(std::string(wp));
}
