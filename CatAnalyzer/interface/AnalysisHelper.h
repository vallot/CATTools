#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TString.h"
#include "TLorentzVector.h"

class AnalysisHelper {

 public:
  static TLorentzVector leafToTLorentzVector(reco::LeafCandidate & leaf){return TLorentzVector(leaf.px(), leaf.py(),leaf.pz(),leaf.energy());}
  static bool ptSorting(reco::LeafCandidate & s1, reco::LeafCandidate & s2) { return ( s1.pt() > s2.pt() ); }

 private:
  
  edm::TriggerNames triggerNames_;
  edm::Handle<edm::TriggerResults> triggerResults_;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  bool triggerInfoSet_;
  
 public:
  
  AnalysisHelper(){triggerInfoSet_ = false;}
  AnalysisHelper(edm::TriggerNames triggerNames, edm::Handle<edm::TriggerResults> triggerResults, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects){triggerNames_ = triggerNames; triggerResults_ = triggerResults; triggerObjects_ = triggerObjects; triggerInfoSet_ = true;}
  ~AnalysisHelper(){}

  bool triggerNotSet();
  bool triggerFired(TString trigname);
  bool triggerMatched(TString trigname, cat::Particle & recoObj, float dR = 0.1);

};
