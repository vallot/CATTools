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

float Jet::smearedRes(int direction, int era) const {
  // The era-based JER is going to be removed
  if ( era == 0 ) return direction == 0 ? fJER_ : direction > 0 ? fJERUp_ : fJERDown_;

  const auto aGenJet = this->genJet();
  if ( !aGenJet ) return 1; // No JER

  const double absEta = std::abs(this->eta());
  if (absEta >=5.0 ) return 1; // No JER
    
  //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
  std::vector<double> etaBins = {0.5, 0.8, 1.1, 1.3, 1.7, 1.9, 2.1, 2.3, 2.5, 2.8, 3.0, 3.2, 5.0};
  std::vector<double> cJERs   = {1.095, 1.120, 1.097, 1.103, 1.118, 1.100, 1.162, 1.160, 1.161, 1.209, 1.564, 1.384, 1.216};
  std::vector<double> cJERsUp = {1.095+0.018, 1.120+0.028, 1.097+0.017, 1.103+0.033, 1.118+0.014,
                                 1.100+0.033, 1.162+0.044, 1.160+0.048, 1.161+0.060, 1.209+0.059,
                                 1.564+0.321, 1.384+0.033, 1.216+0.050};
  std::vector<double> cJERsDn = {1.095-0.018, 1.120-0.028, 1.097-0.017, 1.103-0.033, 1.118-0.014,
                                 1.100-0.033, 1.162-0.044, 1.160-0.048, 1.161-0.060, 1.209-0.059,
                                 1.564-0.321, 1.384-0.033, 1.216-0.050};
  if (era == 2012){
    // 2012 values from https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
    etaBins = {0.5, 1.1, 1.7, 2.3, 2.8, 3.2, 5.0};
    cJERs   = {1.079, 1.099, 1.121, 1.208, 1.254, 1.395, 1.056};
    cJERsUp = {1.105, 1.127, 1.150, 1.254, 1.316, 1.458, 1.247};
    cJERsDn = {1.053, 1.071, 1.092, 1.162, 1.192, 1.332, 0.865};
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

float Jet::scaleFactorCSVv2(Jet::BTAGCSV_CUT cutType, int syst) const {
  if (std::abs(this->eta()) > 2.4 ) return -1; // reject jets out of eta range
  //based on https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation74X50ns
  const double pt = this->pt();
  if ( pt < 20 || pt > 1000 ) return -1;

  int flav = 2;
  if (std::abs(this->partonFlavour()) == 5) flav = 0;
  if (std::abs(this->partonFlavour()) == 4) flav = 1;
  
  if (cutType == 0 && syst == 0 && flav == 1 && pt >= 30 && pt < 670 ) return 0.908299+(2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144)))))));
  if (cutType == 0 && syst == 0 && flav == 0 && pt >= 30 && pt < 670 ) return 0.908299+(2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144)))))));
  if (cutType == 0 && syst == -1 && flav == 1 && pt >= 30 && pt < 50 ) return 0.908299+((2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144)))))))-0.044655226171016693);
  if (cutType == 0 && syst == -1 && flav == 1 && pt >= 50 && pt < 70 ) return 0.908299+((2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144)))))))-0.030660966411232948);
  if (cutType == 0 && syst == -1 && flav == 1 && pt >= 70 && pt < 100 ) return 0.908299+((2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144)))))))-0.048987984657287598);
  if (cutType == 0 && syst == -1 && flav == 1 && pt >= 100 && pt < 140 ) return 0.908299+((2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144)))))))-0.041866477578878403);
  if (cutType == 0 && syst == -1 && flav == 1 && pt >= 140 && pt < 200 ) return 0.908299+((2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144)))))))-0.058439217507839203);
  if (cutType == 0 && syst == -1 && flav == 1 && pt >= 200 && pt < 300 ) return 0.908299+((2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144)))))))-0.079142965376377106);
  if (cutType == 0 && syst == -1 && flav == 1 && pt >= 300 && pt < 670 ) return 0.908299+((2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144)))))))-0.094659518897533417);
  if (cutType == 0 && syst == -1 && flav == 0 && pt >= 30 && pt < 50 ) return 0.908299+((2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144)))))))-0.022327613085508347);
  if (cutType == 0 && syst == -1 && flav == 0 && pt >= 50 && pt < 70 ) return 0.908299+((2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144)))))))-0.015330483205616474);
  if (cutType == 0 && syst == -1 && flav == 0 && pt >= 70 && pt < 100 ) return 0.908299+((2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144)))))))-0.024493992328643799);
  if (cutType == 0 && syst == -1 && flav == 0 && pt >= 100 && pt < 140 ) return 0.908299+((2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144)))))))-0.020933238789439201);
  if (cutType == 0 && syst == -1 && flav == 0 && pt >= 140 && pt < 200 ) return 0.908299+((2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144)))))))-0.029219608753919601);
  if (cutType == 0 && syst == -1 && flav == 0 && pt >= 200 && pt < 300 ) return 0.908299+((2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144)))))))-0.039571482688188553);
  if (cutType == 0 && syst == -1 && flav == 0 && pt >= 300 && pt < 670 ) return 0.908299+((2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144)))))))-0.047329759448766708);
  if (cutType == 0 && syst == +1 && flav == 1 && pt >= 30 && pt < 50 ) return (0.908299+(2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144))))))))+0.044655226171016693;
  if (cutType == 0 && syst == +1 && flav == 1 && pt >= 50 && pt < 70 ) return (0.908299+(2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144))))))))+0.030660966411232948;
  if (cutType == 0 && syst == +1 && flav == 1 && pt >= 70 && pt < 100 ) return (0.908299+(2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144))))))))+0.048987984657287598;
  if (cutType == 0 && syst == +1 && flav == 1 && pt >= 100 && pt < 140 ) return (0.908299+(2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144))))))))+0.041866477578878403;
  if (cutType == 0 && syst == +1 && flav == 1 && pt >= 140 && pt < 200 ) return (0.908299+(2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144))))))))+0.058439217507839203;
  if (cutType == 0 && syst == +1 && flav == 1 && pt >= 200 && pt < 300 ) return (0.908299+(2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144))))))))+0.079142965376377106;
  if (cutType == 0 && syst == +1 && flav == 1 && pt >= 300 && pt < 670 ) return (0.908299+(2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144))))))))+0.094659518897533417;
  if (cutType == 0 && syst == +1 && flav == 0 && pt >= 30 && pt < 50 ) return (0.908299+(2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144))))))))+0.022327613085508347;
  if (cutType == 0 && syst == +1 && flav == 0 && pt >= 50 && pt < 70 ) return (0.908299+(2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144))))))))+0.015330483205616474;
  if (cutType == 0 && syst == +1 && flav == 0 && pt >= 70 && pt < 100 ) return (0.908299+(2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144))))))))+0.024493992328643799;
  if (cutType == 0 && syst == +1 && flav == 0 && pt >= 100 && pt < 140 ) return (0.908299+(2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144))))))))+0.020933238789439201;
  if (cutType == 0 && syst == +1 && flav == 0 && pt >= 140 && pt < 200 ) return (0.908299+(2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144))))))))+0.029219608753919601;
  if (cutType == 0 && syst == +1 && flav == 0 && pt >= 200 && pt < 300 ) return (0.908299+(2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144))))))))+0.039571482688188553;
  if (cutType == 0 && syst == +1 && flav == 0 && pt >= 300 && pt < 670 ) return (0.908299+(2.70877e-06*(log(pt+370.144)*(log(pt+370.144)*(3-(-(104.614*log(pt+370.144))))))))+0.047329759448766708;
  if (cutType == 1 && syst == 0 && flav == 1 && pt >= 30 && pt < 670 ) return -(0.0443172)+(0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85)))))));
  if (cutType == 1 && syst == 0 && flav == 0 && pt >= 30 && pt < 670 ) return -(0.0443172)+(0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85)))))));
  if (cutType == 1 && syst == -1 && flav == 1 && pt >= 30 && pt < 50 ) return -(0.0443172)+((0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85)))))))-0.063294470310211182);
  if (cutType == 1 && syst == -1 && flav == 1 && pt >= 50 && pt < 70 ) return -(0.0443172)+((0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85)))))))-0.043231822550296783);
  if (cutType == 1 && syst == -1 && flav == 1 && pt >= 70 && pt < 100 ) return -(0.0443172)+((0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85)))))))-0.065539278090000153);
  if (cutType == 1 && syst == -1 && flav == 1 && pt >= 100 && pt < 140 ) return -(0.0443172)+((0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85)))))))-0.04837958887219429);
  if (cutType == 1 && syst == -1 && flav == 1 && pt >= 140 && pt < 200 ) return -(0.0443172)+((0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85)))))))-0.087311208248138428);
  if (cutType == 1 && syst == -1 && flav == 1 && pt >= 200 && pt < 300 ) return -(0.0443172)+((0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85)))))))-0.12093273550271988);
  if (cutType == 1 && syst == -1 && flav == 1 && pt >= 300 && pt < 670 ) return -(0.0443172)+((0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85)))))))-0.129528530836105347);
  if (cutType == 1 && syst == -1 && flav == 0 && pt >= 30 && pt < 50 ) return -(0.0443172)+((0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85)))))))-0.031647235155105591);
  if (cutType == 1 && syst == -1 && flav == 0 && pt >= 50 && pt < 70 ) return -(0.0443172)+((0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85)))))))-0.021615911275148392);
  if (cutType == 1 && syst == -1 && flav == 0 && pt >= 70 && pt < 100 ) return -(0.0443172)+((0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85)))))))-0.032769639045000076);
  if (cutType == 1 && syst == -1 && flav == 0 && pt >= 100 && pt < 140 ) return -(0.0443172)+((0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85)))))))-0.024189794436097145);
  if (cutType == 1 && syst == -1 && flav == 0 && pt >= 140 && pt < 200 ) return -(0.0443172)+((0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85)))))))-0.043655604124069214);
  if (cutType == 1 && syst == -1 && flav == 0 && pt >= 200 && pt < 300 ) return -(0.0443172)+((0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85)))))))-0.06046636775135994);
  if (cutType == 1 && syst == -1 && flav == 0 && pt >= 300 && pt < 670 ) return -(0.0443172)+((0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85)))))))-0.064764265418052673);
  if (cutType == 1 && syst == +1 && flav == 1 && pt >= 30 && pt < 50 ) return (-(0.0443172)+(0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85))))))))+0.063294470310211182;
  if (cutType == 1 && syst == +1 && flav == 1 && pt >= 50 && pt < 70 ) return (-(0.0443172)+(0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85))))))))+0.043231822550296783;
  if (cutType == 1 && syst == +1 && flav == 1 && pt >= 70 && pt < 100 ) return (-(0.0443172)+(0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85))))))))+0.065539278090000153;
  if (cutType == 1 && syst == +1 && flav == 1 && pt >= 100 && pt < 140 ) return (-(0.0443172)+(0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85))))))))+0.04837958887219429;
  if (cutType == 1 && syst == +1 && flav == 1 && pt >= 140 && pt < 200 ) return (-(0.0443172)+(0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85))))))))+0.087311208248138428;
  if (cutType == 1 && syst == +1 && flav == 1 && pt >= 200 && pt < 300 ) return (-(0.0443172)+(0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85))))))))+0.12093273550271988;
  if (cutType == 1 && syst == +1 && flav == 1 && pt >= 300 && pt < 670 ) return (-(0.0443172)+(0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85))))))))+0.129528530836105347;
  if (cutType == 1 && syst == +1 && flav == 0 && pt >= 30 && pt < 50 ) return (-(0.0443172)+(0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85))))))))+0.031647235155105591;
  if (cutType == 1 && syst == +1 && flav == 0 && pt >= 50 && pt < 70 ) return (-(0.0443172)+(0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85))))))))+0.021615911275148392;
  if (cutType == 1 && syst == +1 && flav == 0 && pt >= 70 && pt < 100 ) return (-(0.0443172)+(0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85))))))))+0.032769639045000076;
  if (cutType == 1 && syst == +1 && flav == 0 && pt >= 100 && pt < 140 ) return (-(0.0443172)+(0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85))))))))+0.024189794436097145;
  if (cutType == 1 && syst == +1 && flav == 0 && pt >= 140 && pt < 200 ) return (-(0.0443172)+(0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85))))))))+0.043655604124069214;
  if (cutType == 1 && syst == +1 && flav == 0 && pt >= 200 && pt < 300 ) return (-(0.0443172)+(0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85))))))))+0.06046636775135994;
  if (cutType == 1 && syst == +1 && flav == 0 && pt >= 300 && pt < 670 ) return (-(0.0443172)+(0.00496634*(log(pt+1267.85)*(log(pt+1267.85)*(3-(-(0.110428*log(pt+1267.85))))))))+0.064764265418052673;
  if (cutType == 2 && syst == 0 && flav == 1 && pt >= 30 && pt < 670 ) return -(5.1345)+(0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1)))))));
  if (cutType == 2 && syst == 0 && flav == 0 && pt >= 30 && pt < 670 ) return -(5.1345)+(0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1)))))));
  if (cutType == 2 && syst == -1 && flav == 1 && pt >= 30 && pt < 50 ) return -(5.1345)+((0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1)))))))-0.063704296946525574);
  if (cutType == 2 && syst == -1 && flav == 1 && pt >= 50 && pt < 70 ) return -(5.1345)+((0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1)))))))-0.047892197966575623);
  if (cutType == 2 && syst == -1 && flav == 1 && pt >= 70 && pt < 100 ) return -(5.1345)+((0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1)))))))-0.077270857989788055);
  if (cutType == 2 && syst == -1 && flav == 1 && pt >= 100 && pt < 140 ) return -(5.1345)+((0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1)))))))-0.062878459692001343);
  if (cutType == 2 && syst == -1 && flav == 1 && pt >= 140 && pt < 200 ) return -(5.1345)+((0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1)))))))-0.098963312804698944);
  if (cutType == 2 && syst == -1 && flav == 1 && pt >= 200 && pt < 300 ) return -(5.1345)+((0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1)))))))-0.14805065095424652);
  if (cutType == 2 && syst == -1 && flav == 1 && pt >= 300 && pt < 670 ) return -(5.1345)+((0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1)))))))-0.149765393555164337);
  if (cutType == 2 && syst == -1 && flav == 0 && pt >= 30 && pt < 50 ) return -(5.1345)+((0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1)))))))-0.031852148473262787);
  if (cutType == 2 && syst == -1 && flav == 0 && pt >= 50 && pt < 70 ) return -(5.1345)+((0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1)))))))-0.023946098983287811);
  if (cutType == 2 && syst == -1 && flav == 0 && pt >= 70 && pt < 100 ) return -(5.1345)+((0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1)))))))-0.038635428994894028);
  if (cutType == 2 && syst == -1 && flav == 0 && pt >= 100 && pt < 140 ) return -(5.1345)+((0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1)))))))-0.031439229846000671);
  if (cutType == 2 && syst == -1 && flav == 0 && pt >= 140 && pt < 200 ) return -(5.1345)+((0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1)))))))-0.049481656402349472);
  if (cutType == 2 && syst == -1 && flav == 0 && pt >= 200 && pt < 300 ) return -(5.1345)+((0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1)))))))-0.07402532547712326);
  if (cutType == 2 && syst == -1 && flav == 0 && pt >= 300 && pt < 670 ) return -(5.1345)+((0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1)))))))-0.074882696777582169);
  if (cutType == 2 && syst == +1 && flav == 1 && pt >= 30 && pt < 50 ) return (-(5.1345)+(0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1))))))))+0.063704296946525574;
  if (cutType == 2 && syst == +1 && flav == 1 && pt >= 50 && pt < 70 ) return (-(5.1345)+(0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1))))))))+0.047892197966575623;
  if (cutType == 2 && syst == +1 && flav == 1 && pt >= 70 && pt < 100 ) return (-(5.1345)+(0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1))))))))+0.077270857989788055;
  if (cutType == 2 && syst == +1 && flav == 1 && pt >= 100 && pt < 140 ) return (-(5.1345)+(0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1))))))))+0.062878459692001343;
  if (cutType == 2 && syst == +1 && flav == 1 && pt >= 140 && pt < 200 ) return (-(5.1345)+(0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1))))))))+0.098963312804698944;
  if (cutType == 2 && syst == +1 && flav == 1 && pt >= 200 && pt < 300 ) return (-(5.1345)+(0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1))))))))+0.14805065095424652;
  if (cutType == 2 && syst == +1 && flav == 1 && pt >= 300 && pt < 670 ) return (-(5.1345)+(0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1))))))))+0.149765393555164337;
  if (cutType == 2 && syst == +1 && flav == 0 && pt >= 30 && pt < 50 ) return (-(5.1345)+(0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1))))))))+0.031852148473262787;
  if (cutType == 2 && syst == +1 && flav == 0 && pt >= 50 && pt < 70 ) return (-(5.1345)+(0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1))))))))+0.023946098983287811;
  if (cutType == 2 && syst == +1 && flav == 0 && pt >= 70 && pt < 100 ) return (-(5.1345)+(0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1))))))))+0.038635428994894028;
  if (cutType == 2 && syst == +1 && flav == 0 && pt >= 100 && pt < 140 ) return (-(5.1345)+(0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1))))))))+0.031439229846000671;
  if (cutType == 2 && syst == +1 && flav == 0 && pt >= 140 && pt < 200 ) return (-(5.1345)+(0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1))))))))+0.049481656402349472;
  if (cutType == 2 && syst == +1 && flav == 0 && pt >= 200 && pt < 300 ) return (-(5.1345)+(0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1))))))))+0.07402532547712326;
  if (cutType == 2 && syst == +1 && flav == 0 && pt >= 300 && pt < 670 ) return (-(5.1345)+(0.000820101*(log(pt+11518.1)*(log(pt+11518.1)*(3-(-(8.66128*log(pt+11518.1))))))))+0.074882696777582169;
  if (cutType == 0 && syst == 0 && flav == 2 && pt >= 20 && pt < 1000 ) return ((1.07278+(0.000535714*pt))+(-1.14886e-06*(pt*pt)))+(7.0636e-10*(pt*(pt*pt)));
  if (cutType == 0 && syst == -1 && flav == 2 && pt >= 20 && pt < 1000 ) return ((1.01637+(0.000265653*pt))+(-4.22531e-07*(pt*pt)))+(2.23396e-10*(pt*(pt*pt)));
  if (cutType == 0 && syst == +1 && flav == 2 && pt >= 20 && pt < 1000 ) return ((1.12921+(0.000804962*pt))+(-1.87332e-06*(pt*pt)))+(1.18864e-09*(pt*(pt*pt)));
  if (cutType == 1 && syst == 0 && flav == 2 && pt >= 20 && pt < 1000 ) return 1.14022;
  if (cutType == 1 && syst == -1 && flav == 2 && pt >= 20 && pt < 1000 ) return 0.94022;
  if (cutType == 1 && syst == +1 && flav == 2 && pt >= 20 && pt < 1000 ) return 1.34022;
  if (cutType == 2 && syst == 0 && flav == 2 && pt >= 20 && pt < 1000 ) return 0.907317;
  if (cutType == 2 && syst == -1 && flav == 2 && pt >= 20 && pt < 1000 ) return 0.557317;
  if (cutType == 2 && syst == +1 && flav == 2 && pt >= 20 && pt < 1000 ) return 1.257317;
  return -1;
}
