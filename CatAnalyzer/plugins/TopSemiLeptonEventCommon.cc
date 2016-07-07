#include "CATTools/CatAnalyzer/interface/TopSemiLeptonEventCommon.h"

#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CATTools/CatAnalyzer/interface/TopTriggerSF.h"
#include "CATTools/CatAnalyzer/interface/TopEventGlobalVar.h"
#include "CATTools/CatAnalyzer/interface/KinematicReconstructionSolution.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "CATTools/CatAnalyzer/interface/TTSemiLeptonEventSelector.h"

using namespace std;
using namespace cat;
using namespace TopEventCommonGlobal;


TopSemiLeptonEventCommon::TopSemiLeptonEventCommon(const edm::ParameterSet& iConfig ) :TopEventCommon(iConfig) 
{
  eventSelect_ = new TTSemiLeptonEventSelector(iConfig, consumesCollector()) ;
  paramInit(iConfig);
}

void TopSemiLeptonEventCommon::showSummary() {
  TopEventInfo& evInfo_ = TopEventInfo::getInstance();
  cout <<setw(10)<<"cut flow"<<setw(10)<<"no ll"<<setw(10)<<"mu Jet"<<setw(10)<<"e Jet"<<setw(10)<<"all"<< endl;
  for ( int i=0; i<NCutflow; ++i ) {
    cout <<setw(10)<<"step "<<i<< setw(10)<<evInfo_.cutflow_[i][0] << setw(10)<<evInfo_.cutflow_[i][1] << setw(10)<<evInfo_.cutflow_[i][2] <<setw(10)<<evInfo_.cutflow_[i][3]<< endl;
  }
}

void TopSemiLeptonEventCommon::analyzeCustom(const edm::Event& iEvent, const edm::EventSetup& iSetup, int sys) {
  std::cout<<"SemiLepton Analyzer"<<std::endl;
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TopSemiLeptonEventCommon);
