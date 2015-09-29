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

class AnalysisHelper {
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
  bool triggerFired(const std::string& trigname);
  bool triggerMatched(const std::string& trigname, const cat::Particle & recoObj, const float dR = 0.1);

};

}

#endif
