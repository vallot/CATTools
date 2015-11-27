#include "CATTools/CommonTools/interface/AnalysisHelper.h"

using namespace std;

namespace cat
{
  TLorentzVector ToTLorentzVector(const reco::Candidate& t) { return TLorentzVector(t.px(), t.py(), t.pz(), t.energy()); }
  bool GtByTLVPt( TLorentzVector & t1, TLorentzVector & t2 ){ return t1.Pt() > t2.Pt();}

  bool AnalysisHelper::triggerNotSet()
  {
    cout <<"trigger info not set, use the constructor below to set the trigger info "<< endl;
    cout <<"AnalysisHelper(edm::TriggerNames &triggerNames, edm::Handle<edm::TriggerResults> triggerResults, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects)"<< endl;
    return false;
  }

  void AnalysisHelper::listFiredTriggers()
  {
    if (!triggerInfoSet_) return;

    const unsigned int ntrigs = triggerResults_->size();
    for (unsigned int itr=0; itr<ntrigs; itr++){
      const string& trigName = triggerNames_.triggerName(itr);
      if (!triggerResults_->accept(itr)) continue;
      std::cout << trigName << std::endl;
    }
    return;
  }

  bool AnalysisHelper::triggerFired(const std::string& trigname)
  {
    if (!triggerInfoSet_) return triggerNotSet();

    const unsigned int ntrigs = triggerResults_->size();
    for (unsigned int itr=0; itr<ntrigs; itr++){
      const string& trigName = triggerNames_.triggerName(itr);
      if (!triggerResults_->accept(itr)) continue;
      if(trigName.find(trigname) != std::string::npos)      return true;
    }
    return false;
  }

  bool AnalysisHelper::triggerMatched(const std::string& trigname, const cat::Particle & recoObj, const float dR)
  {
    if (!triggerInfoSet_) return triggerNotSet();

    for (pat::TriggerObjectStandAlone trigObj : *triggerObjects_) { // note: not "const &" since we want to call unpackPathNames
      trigObj.unpackPathNames(triggerNames_);
      std::vector<std::string> pathNamesAll  = trigObj.pathNames(false);
      for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
	if ( pathNamesAll[h].find(trigname) == 0 ){
	  if (trigObj.hasPathName( pathNamesAll[h], true, true )){
	    // found trigger
	    if ( reco::deltaR(trigObj, recoObj) < dR){
	      // found matching trigger
	      //std::cout << "\tTrigger trigObject:  pt " << trigObj.pt() << ", eta " << trigObj.eta() << ", phi " << trigObj.phi() << std::endl;
	      return true;
	    }
	  }
	}
      }
    }
    return false;
  }

  math::XYZTLorentzVector getLVFromPtPhi(const double pt, const double phi)
  {
    const double px = pt*std::cos(phi);
    const double py = pt*std::sin(phi);
    return math::XYZTLorentzVector(px, py, 0, pt);
  }

}
