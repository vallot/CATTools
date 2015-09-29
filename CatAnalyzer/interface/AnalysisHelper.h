#ifndef CATTools_CatAnalysis_AnalysisHelper_H
#define CATTools_CatAnalysis_AnalysisHelper_H

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TString.h"
#include "TLorentzVector.h"

class AnalysisHelper {

 public:
  static TLorentzVector leafToTLorentzVector(const reco::LeafCandidate & leaf){return TLorentzVector(leaf.px(), leaf.py(),leaf.pz(),leaf.energy());}
  static bool ptSorting(const reco::LeafCandidate & s1, const reco::LeafCandidate & s2) { return ( s1.pt() > s2.pt() ); }

 private:
  
  edm::TriggerNames triggerNames_;
  edm::Handle<edm::TriggerResults> triggerResults_;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  bool triggerInfoSet_;
  
 public:
  
  AnalysisHelper(){triggerInfoSet_ = false;}
  AnalysisHelper(const edm::TriggerNames& triggerNames, edm::Handle<edm::TriggerResults> triggerResults, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects){triggerNames_ = triggerNames; triggerResults_ = triggerResults; triggerObjects_ = triggerObjects; triggerInfoSet_ = true;}
  ~AnalysisHelper(){}

  bool triggerNotSet();
  bool triggerFired(TString trigname);
  bool triggerMatched(TString trigname, cat::Particle & recoObj, float dR = 0.1);

};

namespace cat {

math::XYZTLorentzVector getLVFromPtPhi(const double pt, const double phi);

template<typename T1, typename T2>
math::XYZTLorentzVector computeMETVariation(const math::XYZTLorentzVector& metP4, T1 corrBegin, T1 corrEnd, T2 rawBegin, T2 rawEnd)
{
  // Calculate effect of energy correction
  // correctedPx - rawPx, correctedPy - rawPy of visible particle are calculated
  double dpx = 0, dpy = 0;
  for ( auto x = corrBegin; x != corrEnd; ++x ) { dpx += x->px(); dpy += x->py(); }
  for ( auto x = rawBegin; x != rawEnd; ++x ) { dpx -= x->px(); dpy -= x->py(); }

  // MET is shifted by -(dpx, dpy)
  const double px = metP4.px() - dpx;
  const double py = metP4.py() - dpy;
  const double energy = std::hypot(px, py);

  return math::XYZTLorentzVector(px, py, 0, energy);
};

}

#endif
