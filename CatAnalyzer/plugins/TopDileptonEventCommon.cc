#include "CATTools/CatAnalyzer/interface/TopDiLeptonEventCommon.h"

#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CATTools/CatAnalyzer/interface/TopTriggerSF.h"
#include "CATTools/CatAnalyzer/interface/TopEventGlobalVar.h"
#include "CATTools/CatAnalyzer/interface/KinematicReconstructionSolution.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "CATTools/CatAnalyzer/interface/TTDiLeptonEventSelector.h"

using namespace std;
using namespace cat;
using namespace TopEventCommonGlobal;


TopDiLeptonEventCommon::TopDiLeptonEventCommon(const edm::ParameterSet& iConfig ):TopEventCommon(iConfig)
{
  eventSelect_ = new TTDiLeptonEventSelector(iConfig, consumesCollector()) ;
  paramInit(iConfig);
}

void TopDiLeptonEventCommon::showSummary() {
  TopEventInfo& evInfo_ = TopEventInfo::getInstance();
  cout <<setw(10)<<"cut flow"<<setw(10)<<"no ll"<<setw(10)<<"emu"<<setw(10)<<"ee"<<setw(10)<<"mumu"<< endl;
  for ( int i=0; i<NCutflow; ++i ) {
    cout <<setw(10)<<"step "<<i<< setw(10)<<evInfo_.cutflow_[i][0] << setw(10)<<evInfo_.cutflow_[i][1] << setw(10)<<evInfo_.cutflow_[i][2] <<setw(10)<<evInfo_.cutflow_[i][3]<< endl;
  }
}
void TopDiLeptonEventCommon::analyzeCustom(const edm::Event& iEvent, const edm::EventSetup& iSetup, int sys) {
  std::cout<<"DiLepton Analyzer"<<std::endl;
}
//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TopDiLeptonEventCommon);
