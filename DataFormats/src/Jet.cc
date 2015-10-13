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

float Jet::scaleFactorCSVv2(Jet::BTAGCSV_CUT cutType, int syst, JETFLAV flav) const {
  if (std::abs(this->eta()) > 2.4 ) return -1; // reject jets out of eta range
  //based on https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation74X50ns
  const double pt = this->pt();
  if ( pt < 30 || pt >= 670 ) return -1;

  if ( flav == JETFLAV_LIGHT ) {
    switch ( cutType ) {
      case BTAGCSV_LOOSE : return 1+syst*0.15; break;
      case BTAGCSV_MEDIUM: return 1+syst*0.20; break;
      case BTAGCSV_TIGHT : return 1+syst*0.40; break;
      default: return -1;
    }
  }

  if ( cutType == BTAGCSV_LOOSE ) {
    const double sf0 = 0.901434*(1.+(0.0852659*pt))/(1.+(0.0769021*pt));
    if ( syst == 0 ) return sf0;

    if ( flav == JETFLAV_C ) {
      if      ( pt <  50 ) return sf0+syst*0.10247411578893661;
      else if ( pt <  70 ) return sf0+syst*0.09483686089515686;
      else if ( pt < 100 ) return sf0+syst*0.075944989919662476;
      else if ( pt < 140 ) return sf0+syst*0.065169334411621094;
      else if ( pt < 200 ) return sf0+syst*0.12829481065273285;
      else if ( pt < 300 ) return sf0+syst*0.21497605741024017;
      else                return sf0+syst*0.22071903944015503;
    }
    else if ( flav == JETFLAV_B ) {
      if      ( pt <  50 ) return sf0+syst*0.051237057894468307;
      else if ( pt <  70 ) return sf0+syst*0.04741843044757843;
      else if ( pt < 100 ) return sf0+syst*0.037972494959831238;
      else if ( pt < 140 ) return sf0+syst*0.032584667205810547;
      else if ( pt < 200 ) return sf0+syst*0.064147405326366425;
      else if ( pt < 300 ) return sf0+syst*0.10748802870512009;
      else                 return sf0+syst*0.11035951972007751;
    }
  }
  else if ( cutType == BTAGCSV_MEDIUM ) {
    const double sf0 = 0.968546;
    if (syst == 0 ) return sf0;

    if ( flav == JETFLAV_C ) {
      if      ( pt <  50 ) return sf0+syst*0.093810759484767914;
      else if ( pt <  70 ) return sf0+syst*0.10516738146543503;
      else if ( pt < 100 ) return sf0+syst*0.095603741705417633;
      else if ( pt < 140 ) return sf0+syst*0.11572366207838058;
      else if ( pt < 200 ) return sf0+syst*0.19687847793102264;
      else if ( pt < 300 ) return sf0+syst*0.16733628511428833;
      else                 return sf0+syst*0.25004029273986816;
    }
    else if ( flav == JETFLAV_B ) {
      if      ( pt <  50 ) return sf0+syst*0.046905379742383957;
      else if ( pt <  70 ) return sf0+syst*0.052583690732717514;
      else if ( pt < 100 ) return sf0+syst*0.047801870852708817;
      else if ( pt < 140 ) return sf0+syst*0.057861831039190292;
      else if ( pt < 200 ) return sf0+syst*0.098439238965511322;
      else if ( pt < 300 ) return sf0+syst*0.083668142557144165;
      else                 return sf0+syst*0.12502014636993408;
    }
  }
  else if ( cutType == BTAGCSV_TIGHT ) {
    const double sf0 = 0.924703;
    if ( syst == 0 ) return sf0;

    if ( flav == JETFLAV_C ) {
      if      ( pt <  50 ) return sf0+syst*0.077449031174182892;
      else if ( pt <  70 ) return sf0+syst*0.10176990181207657;
      else if ( pt < 100 ) return sf0+syst*0.079285271465778351;
      else if ( pt < 140 ) return sf0+syst*0.19293521344661713;
      else if ( pt < 200 ) return sf0+syst*0.1810869574546814;
      else if ( pt < 300 ) return sf0+syst*0.5206335186958313;
      else                 return sf0+syst*0.42109844088554382;
    }
    else if ( flav == JETFLAV_B ) {
      if      ( pt <  50 ) return sf0+syst*0.038724515587091446;
      else if ( pt <  70 ) return sf0+syst*0.050884950906038284;
      else if ( pt < 100 ) return sf0+syst*0.039642635732889175;
      else if ( pt < 140 ) return sf0+syst*0.096467606723308563;
      else if ( pt < 200 ) return sf0+syst*0.090543478727340698;
      else if ( pt < 300 ) return sf0+syst*0.26031675934791565;
      else                 return sf0+syst*0.21054922044277191;
    }
  }
  return -1;
}
