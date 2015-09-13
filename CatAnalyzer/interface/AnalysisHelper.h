#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "CATTools/DataFormats/interface/Particle.h"

#include "TString.h"
#include "TLorentzVector.h"

class AnalysisHelper {

 public:
  static float deltaR(float e1,float p1,float e2,float p2);
  static float deltaR(reco::LeafCandidate & l1, reco::LeafCandidate & l2){return deltaR(l1.eta(),l1.phi(),l2.eta(),l2.phi());}
  static TLorentzVector leafToTLorentzVector(reco::LeafCandidate & leaf){return TLorentzVector(leaf.px(), leaf.py(),leaf.pz(),leaf.energy());}

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
