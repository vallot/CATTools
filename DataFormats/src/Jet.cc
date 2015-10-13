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

float Jet::scaleFactorCSVv2(int op, int sys, int flav) const {
  //based on https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation74X50ns
  if (op == 0 && sys == 0 && flav == 1 && this->pt() > 30 && this->pt() < 670 ) return 0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt())));
  if (op == 0 && sys == 0 && flav == 0 && this->pt() > 30 && this->pt() < 670 ) return 0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt())));
  if (op == 0 && sys == -1 && flav == 1 && this->pt() > 30 && this->pt() < 50 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))-0.10247411578893661;
  if (op == 0 && sys == -1 && flav == 1 && this->pt() > 50 && this->pt() < 70 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))-0.09483686089515686;
  if (op == 0 && sys == -1 && flav == 1 && this->pt() > 70 && this->pt() < 100 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))-0.075944989919662476;
  if (op == 0 && sys == -1 && flav == 1 && this->pt() > 100 && this->pt() < 140 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))-0.065169334411621094;
  if (op == 0 && sys == -1 && flav == 1 && this->pt() > 140 && this->pt() < 200 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))-0.12829481065273285;
  if (op == 0 && sys == -1 && flav == 1 && this->pt() > 200 && this->pt() < 300 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))-0.21497605741024017;
  if (op == 0 && sys == -1 && flav == 1 && this->pt() > 300 && this->pt() < 670 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))-0.22071903944015503;
  if (op == 0 && sys == -1 && flav == 0 && this->pt() > 30 && this->pt() < 50 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))-0.051237057894468307;
  if (op == 0 && sys == -1 && flav == 0 && this->pt() > 50 && this->pt() < 70 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))-0.04741843044757843;
  if (op == 0 && sys == -1 && flav == 0 && this->pt() > 70 && this->pt() < 100 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))-0.037972494959831238;
  if (op == 0 && sys == -1 && flav == 0 && this->pt() > 100 && this->pt() < 140 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))-0.032584667205810547;
  if (op == 0 && sys == -1 && flav == 0 && this->pt() > 140 && this->pt() < 200 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))-0.064147405326366425;
  if (op == 0 && sys == -1 && flav == 0 && this->pt() > 200 && this->pt() < 300 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))-0.10748802870512009;
  if (op == 0 && sys == -1 && flav == 0 && this->pt() > 300 && this->pt() < 670 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))-0.11035951972007751;
  if (op == 0 && sys == 1 && flav == 1 && this->pt() > 30 && this->pt() < 50 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))+0.10247411578893661;
  if (op == 0 && sys == 1 && flav == 1 && this->pt() > 50 && this->pt() < 70 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))+0.09483686089515686;
  if (op == 0 && sys == 1 && flav == 1 && this->pt() > 70 && this->pt() < 100 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))+0.075944989919662476;
  if (op == 0 && sys == 1 && flav == 1 && this->pt() > 100 && this->pt() < 140 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))+0.065169334411621094;
  if (op == 0 && sys == 1 && flav == 1 && this->pt() > 140 && this->pt() < 200 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))+0.12829481065273285;
  if (op == 0 && sys == 1 && flav == 1 && this->pt() > 200 && this->pt() < 300 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))+0.21497605741024017;
  if (op == 0 && sys == 1 && flav == 1 && this->pt() > 300 && this->pt() < 670 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))+0.22071903944015503;
  if (op == 0 && sys == 1 && flav == 0 && this->pt() > 30 && this->pt() < 50 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))+0.051237057894468307;
  if (op == 0 && sys == 1 && flav == 0 && this->pt() > 50 && this->pt() < 70 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))+0.04741843044757843;
  if (op == 0 && sys == 1 && flav == 0 && this->pt() > 70 && this->pt() < 100 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))+0.037972494959831238;
  if (op == 0 && sys == 1 && flav == 0 && this->pt() > 100 && this->pt() < 140 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))+0.032584667205810547;
  if (op == 0 && sys == 1 && flav == 0 && this->pt() > 140 && this->pt() < 200 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))+0.064147405326366425;
  if (op == 0 && sys == 1 && flav == 0 && this->pt() > 200 && this->pt() < 300 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))+0.10748802870512009;
  if (op == 0 && sys == 1 && flav == 0 && this->pt() > 300 && this->pt() < 670 ) return (0.901434*((1.+(0.0852659*this->pt()))/(1.+(0.0769021*this->pt()))))+0.11035951972007751;
  if (op == 1 && sys == 0 && flav == 1 && this->pt() > 30 && this->pt() < 670 ) return 0.968546;
  if (op == 1 && sys == 0 && flav == 0 && this->pt() > 30 && this->pt() < 670 ) return 0.968546;
  if (op == 1 && sys == -1 && flav == 1 && this->pt() > 30 && this->pt() < 50 ) return 0.968546-0.093810759484767914;
  if (op == 1 && sys == -1 && flav == 1 && this->pt() > 50 && this->pt() < 70 ) return 0.968546-0.10516738146543503;
  if (op == 1 && sys == -1 && flav == 1 && this->pt() > 70 && this->pt() < 100 ) return 0.968546-0.095603741705417633;
  if (op == 1 && sys == -1 && flav == 1 && this->pt() > 100 && this->pt() < 140 ) return 0.968546-0.11572366207838058;
  if (op == 1 && sys == -1 && flav == 1 && this->pt() > 140 && this->pt() < 200 ) return 0.968546-0.19687847793102264;
  if (op == 1 && sys == -1 && flav == 1 && this->pt() > 200 && this->pt() < 300 ) return 0.968546-0.16733628511428833;
  if (op == 1 && sys == -1 && flav == 1 && this->pt() > 300 && this->pt() < 670 ) return 0.968546-0.25004029273986816;
  if (op == 1 && sys == -1 && flav == 0 && this->pt() > 30 && this->pt() < 50 ) return 0.968546-0.046905379742383957;
  if (op == 1 && sys == -1 && flav == 0 && this->pt() > 50 && this->pt() < 70 ) return 0.968546-0.052583690732717514;
  if (op == 1 && sys == -1 && flav == 0 && this->pt() > 70 && this->pt() < 100 ) return 0.968546-0.047801870852708817;
  if (op == 1 && sys == -1 && flav == 0 && this->pt() > 100 && this->pt() < 140 ) return 0.968546-0.057861831039190292;
  if (op == 1 && sys == -1 && flav == 0 && this->pt() > 140 && this->pt() < 200 ) return 0.968546-0.098439238965511322;
  if (op == 1 && sys == -1 && flav == 0 && this->pt() > 200 && this->pt() < 300 ) return 0.968546-0.083668142557144165;
  if (op == 1 && sys == -1 && flav == 0 && this->pt() > 300 && this->pt() < 670 ) return 0.968546-0.12502014636993408;
  if (op == 1 && sys == 1 && flav == 1 && this->pt() > 30 && this->pt() < 50 ) return 0.968546+0.093810759484767914;
  if (op == 1 && sys == 1 && flav == 1 && this->pt() > 50 && this->pt() < 70 ) return 0.968546+0.10516738146543503;
  if (op == 1 && sys == 1 && flav == 1 && this->pt() > 70 && this->pt() < 100 ) return 0.968546+0.095603741705417633;
  if (op == 1 && sys == 1 && flav == 1 && this->pt() > 100 && this->pt() < 140 ) return 0.968546+0.11572366207838058;
  if (op == 1 && sys == 1 && flav == 1 && this->pt() > 140 && this->pt() < 200 ) return 0.968546+0.19687847793102264;
  if (op == 1 && sys == 1 && flav == 1 && this->pt() > 200 && this->pt() < 300 ) return 0.968546+0.16733628511428833;
  if (op == 1 && sys == 1 && flav == 1 && this->pt() > 300 && this->pt() < 670 ) return 0.968546+0.25004029273986816;
  if (op == 1 && sys == 1 && flav == 0 && this->pt() > 30 && this->pt() < 50 ) return 0.968546+0.046905379742383957;
  if (op == 1 && sys == 1 && flav == 0 && this->pt() > 50 && this->pt() < 70 ) return 0.968546+0.052583690732717514;
  if (op == 1 && sys == 1 && flav == 0 && this->pt() > 70 && this->pt() < 100 ) return 0.968546+0.047801870852708817;
  if (op == 1 && sys == 1 && flav == 0 && this->pt() > 100 && this->pt() < 140 ) return 0.968546+0.057861831039190292;
  if (op == 1 && sys == 1 && flav == 0 && this->pt() > 140 && this->pt() < 200 ) return 0.968546+0.098439238965511322;
  if (op == 1 && sys == 1 && flav == 0 && this->pt() > 200 && this->pt() < 300 ) return 0.968546+0.083668142557144165;
  if (op == 1 && sys == 1 && flav == 0 && this->pt() > 300 && this->pt() < 670 ) return 0.968546+0.12502014636993408;
  if (op == 2 && sys == 0 && flav == 1 && this->pt() > 30 && this->pt() < 670 ) return 0.924703;
  if (op == 2 && sys == 0 && flav == 0 && this->pt() > 30 && this->pt() < 670 ) return 0.924703;
  if (op == 2 && sys == -1 && flav == 1 && this->pt() > 30 && this->pt() < 50 ) return 0.924703-0.077449031174182892;
  if (op == 2 && sys == -1 && flav == 1 && this->pt() > 50 && this->pt() < 70 ) return 0.924703-0.10176990181207657;
  if (op == 2 && sys == -1 && flav == 1 && this->pt() > 70 && this->pt() < 100 ) return 0.924703-0.079285271465778351;
  if (op == 2 && sys == -1 && flav == 1 && this->pt() > 100 && this->pt() < 140 ) return 0.924703-0.19293521344661713;
  if (op == 2 && sys == -1 && flav == 1 && this->pt() > 140 && this->pt() < 200 ) return 0.924703-0.1810869574546814;
  if (op == 2 && sys == -1 && flav == 1 && this->pt() > 200 && this->pt() < 300 ) return 0.924703-0.5206335186958313;
  if (op == 2 && sys == -1 && flav == 1 && this->pt() > 300 && this->pt() < 670 ) return 0.924703-0.42109844088554382;
  if (op == 2 && sys == -1 && flav == 0 && this->pt() > 30 && this->pt() < 50 ) return 0.924703-0.038724515587091446;
  if (op == 2 && sys == -1 && flav == 0 && this->pt() > 50 && this->pt() < 70 ) return 0.924703-0.050884950906038284;
  if (op == 2 && sys == -1 && flav == 0 && this->pt() > 70 && this->pt() < 100 ) return 0.924703-0.039642635732889175;
  if (op == 2 && sys == -1 && flav == 0 && this->pt() > 100 && this->pt() < 140 ) return 0.924703-0.096467606723308563;
  if (op == 2 && sys == -1 && flav == 0 && this->pt() > 140 && this->pt() < 200 ) return 0.924703-0.090543478727340698;
  if (op == 2 && sys == -1 && flav == 0 && this->pt() > 200 && this->pt() < 300 ) return 0.924703-0.26031675934791565;
  if (op == 2 && sys == -1 && flav == 0 && this->pt() > 300 && this->pt() < 670 ) return 0.924703-0.21054922044277191;
  if (op == 2 && sys == 1 && flav == 1 && this->pt() > 30 && this->pt() < 50 ) return 0.924703+0.077449031174182892;
  if (op == 2 && sys == 1 && flav == 1 && this->pt() > 50 && this->pt() < 70 ) return 0.924703+0.10176990181207657;
  if (op == 2 && sys == 1 && flav == 1 && this->pt() > 70 && this->pt() < 100 ) return 0.924703+0.079285271465778351;
  if (op == 2 && sys == 1 && flav == 1 && this->pt() > 100 && this->pt() < 140 ) return 0.924703+0.19293521344661713;
  if (op == 2 && sys == 1 && flav == 1 && this->pt() > 140 && this->pt() < 200 ) return 0.924703+0.1810869574546814;
  if (op == 2 && sys == 1 && flav == 1 && this->pt() > 200 && this->pt() < 300 ) return 0.924703+0.5206335186958313;
  if (op == 2 && sys == 1 && flav == 1 && this->pt() > 300 && this->pt() < 670 ) return 0.924703+0.42109844088554382;
  if (op == 2 && sys == 1 && flav == 0 && this->pt() > 30 && this->pt() < 50 ) return 0.924703+0.038724515587091446;
  if (op == 2 && sys == 1 && flav == 0 && this->pt() > 50 && this->pt() < 70 ) return 0.924703+0.050884950906038284;
  if (op == 2 && sys == 1 && flav == 0 && this->pt() > 70 && this->pt() < 100 ) return 0.924703+0.039642635732889175;
  if (op == 2 && sys == 1 && flav == 0 && this->pt() > 100 && this->pt() < 140 ) return 0.924703+0.096467606723308563;
  if (op == 2 && sys == 1 && flav == 0 && this->pt() > 140 && this->pt() < 200 ) return 0.924703+0.090543478727340698;
  if (op == 2 && sys == 1 && flav == 0 && this->pt() > 200 && this->pt() < 300 ) return 0.924703+0.26031675934791565;
  if (op == 2 && sys == 1 && flav == 0 && this->pt() > 300 && this->pt() < 670 ) return 0.924703+0.21054922044277191;
  if (op == 0 && sys == 0 && flav == 2 && this->pt() > 30 && this->pt() < 670 ) return 1.;
  if (op == 0 && sys == 1 && flav == 2 && this->pt() > 30 && this->pt() < 670 ) return 1.+0.15;
  if (op == 0 && sys == -1 && flav == 2 && this->pt() > 30 && this->pt() < 670 ) return 1.-0.15;
  if (op == 1 && sys == 0 && flav == 2 && this->pt() > 30 && this->pt() < 670 ) return 1.;
  if (op == 1 && sys == 1 && flav == 2 && this->pt() > 30 && this->pt() < 670 ) return 1.+0.20;
  if (op == 1 && sys == -1 && flav == 2 && this->pt() > 30 && this->pt() < 670 ) return 1.-0.20;
  if (op == 2 && sys == 0 && flav == 2 && this->pt() > 30 && this->pt() < 670 ) return 1.;
  if (op == 2 && sys == 1 && flav == 2 && this->pt() > 30 && this->pt() < 670 ) return 1.+0.40;
  if (op == 2 && sys == -1 && flav == 2 && this->pt() > 30 && this->pt() < 670 ) return 1.-0.40;
  return -1;
}
