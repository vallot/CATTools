#ifndef __CATTools_TopEventCommon__
#define __CATTools_TopEventCommon__

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
/*

#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"

#include "CATTools/CommonTools/interface/TTbarModeDefs.h"
#include "CATTools/CommonTools/interface/ScaleFactorEvaluator.h"
#include "CATTools/CatAnalyzer/interface/BTagWeightEvaluator.h"
#include "CATTools/CatAnalyzer/interface/KinematicSolvers.h"

#include "CATTools/CommonTools/interface/AnalysisHelper.h"
#include "CATTools/CatAnalyzer/interface/analysisUtils.h"
#include "CATTools/CatAnalyzer/interface/KinematicReconstruction.h"
*/
#include "TTree.h"
#include "TH1D.h"

#include "CATTools/CatAnalyzer/interface/TopEventInfo.h"
#include "CATTools/CatAnalyzer/interface/TTEventSelector.h"

class TopEventCommon : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchLuminosityBlocks> {
public:
  explicit TopEventCommon(const edm::ParameterSet&);
  ~TopEventCommon();
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void analyzeCustom(const edm::Event&, const edm::EventSetup&, int sys ) ;
  void genInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup );
  void setBranch(TTree* tree, int sys);
  void setBranchCustom(TTree* tree, int sys);
  /*
  void runParameterInit(const edm::ParameterSet& iConfig) { 
    eventSelect_->parameterInit(iConfig);
  }
  */

  void resetBr() {
    evInfo_.resetBr();
    resetCustomBranch();
  }
  virtual void resetCustomBranch() {}
  int runEventSelection(const edm::Event& iEvent, const edm::EventSetup& iSetup, TTree* tree) {
   return eventSelect_->eventSelection(iEvent, iSetup, tree);
  }
  void setEventSelection( TTEventSelector* eventSelect) {
    eventSelect_ = eventSelect;
  }
  void showSummary();
// Use protect keyword for branches.
protected : 
  TopEventInfo evInfo_;
  TTEventSelector* eventSelect_;
  // I/O variables.
  std::vector<TTree*> ttree_;
  TH1D * h_nevents;

  // Token for gen Information
  edm::EDGetTokenT<float> genWeightToken_;
  edm::EDGetTokenT<std::vector<float>> pdfweightToken_, scaleupweightsToken_, scaledownweightsToken_;
  edm::EDGetTokenT<float> puweightToken_, puweightToken_up_, puweightToken_dn_, topPtWeight_;

  edm::EDGetTokenT<int>          partonTop_channel_;
  edm::EDGetTokenT<std::vector<int> > partonTop_modes_;
  edm::EDGetTokenT<reco::GenParticleCollection> partonTop_genParticles_;
  edm::EDGetTokenT<edm::View<reco::Candidate> > pseudoTop_leptons_, pseudoTop_neutrinos_, pseudoTop_jets_;

private:
  void beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&);
  void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override {};


};


#endif
