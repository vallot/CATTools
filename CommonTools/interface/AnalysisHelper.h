#ifndef CATTools_CommonTools_AnalysisHelper_H
#define CATTools_CommonTools_AnalysisHelper_H

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Utils/interface/PtComparator.h"

#include "TLorentzVector.h"

namespace cat {

TLorentzVector ToTLorentzVector(const reco::Candidate& t);
typedef GreaterByPt<reco::Candidate> GtByCandPt;
bool GtByTLVPt( TLorentzVector & t1, TLorentzVector & t2 );

class AnalysisHelper {
 public:
  AnalysisHelper(){triggerInfoSet_ = false;}
  AnalysisHelper(const edm::TriggerNames& triggerNames, edm::Handle<edm::TriggerResults> triggerResults, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects){triggerNames_ = triggerNames; triggerResults_ = triggerResults; triggerObjects_ = triggerObjects; triggerInfoSet_ = true;}
  ~AnalysisHelper(){}

  bool triggerNotSet();
  bool triggerFired(const std::string& trigname);
  bool triggerMatched(const std::string& trigname, const cat::Particle & recoObj, const float dR = 0.1);
  void listFiredTriggers();

 private:

  edm::TriggerNames triggerNames_;
  edm::Handle<edm::TriggerResults> triggerResults_;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  bool triggerInfoSet_;
};

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
