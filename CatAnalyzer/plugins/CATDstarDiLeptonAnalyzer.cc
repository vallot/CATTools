
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CATTools/CatAnalyzer/interface/TopTriggerSF.h"
#include "CATTools/CatAnalyzer/interface/TopEventGlobalVar.h"
#include "CATTools/CatAnalyzer/interface/KinematicReconstructionSolution.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "CATTools/CatAnalyzer/interface/TTDiLeptonEventSelector.h"
#include "CATTools/CatAnalyzer/interface/CATDstarSemiLeptonAnalyzer.h"

using namespace std;
using namespace cat;
using namespace TopEventCommonGlobal;

class CATDstarDiLeptonAnalyzer : public CATDstarSemiLeptonAnalyzer {
public:
  explicit CATDstarDiLeptonAnalyzer(const edm::ParameterSet&);
  virtual ~CATDstarDiLeptonAnalyzer() { showSummary(); }
protected : 
  

private:

};

CATDstarDiLeptonAnalyzer::CATDstarDiLeptonAnalyzer(const edm::ParameterSet& iConfig ) :CATDstarSemiLeptonAnalyzer(iConfig) 
{
  eventSelect_ = new TTDiLeptonEventSelector(iConfig, consumesCollector()) ;
  paramInit(iConfig);
  d0Token_  = consumes<cat::SecVertexCollection>(iConfig.getParameter<edm::InputTag>("d0s"));
  dstarToken_  = consumes<cat::SecVertexCollection>(iConfig.getParameter<edm::InputTag>("dstars"));
  mcSrc_ = consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("mcLabel"));
  matchingDeltaR_  = iConfig.getParameter<double>("matchingDeltaR");

  for (int sys = 0; sys < nsys_e; ++sys){
    auto tr = ttree_[sys];
    setBranchCustom(tr, sys);
  }
}



//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CATDstarDiLeptonAnalyzer);
