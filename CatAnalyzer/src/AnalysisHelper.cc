#include "CATTools/CatAnalyzer/interface/AnalysisHelper.h"

using namespace std;

bool AnalysisHelper::triggerNotSet()
{
  cout <<"trigger info not set, use the constructor below to set the trigger info "<< endl;
  cout <<"AnalysisHelper(edm::TriggerNames &triggerNames, edm::Handle<edm::TriggerResults> triggerResults, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects)"<< endl;
  return false;
}

bool AnalysisHelper::triggerFired(TString trigname)
{
  if (!triggerInfoSet_) return triggerNotSet();
  
  const unsigned int ntrigs = triggerResults_->size();
  for (unsigned int itr=0; itr<ntrigs; itr++){
    TString trigName=triggerNames_.triggerName(itr);
    if (!triggerResults_->accept(itr)) continue;
    if(trigName.Contains(trigname))      return true;
  }
  return false;
}

bool AnalysisHelper::triggerMatched(TString trigname, cat::Particle & recoObj, float dR)
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
