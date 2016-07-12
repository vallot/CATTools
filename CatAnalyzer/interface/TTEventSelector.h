#ifndef __CATTools_CatAnalyzer_TTEventSelector__
#define __CATTools_CatAnalyzer_TTEventSelector__

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

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
#include "TTree.h"
#include "TH1D.h"

#include "CATTools/CatAnalyzer/interface/TopTriggerSF.h"
#include "CATTools/CatAnalyzer/interface/KinematicReconstructionSolution.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "CATTools/CatAnalyzer/interface/TopEventInfo.h"
#include "CATTools/CatAnalyzer/interface/TopEventGlobalVar.h"

class TTEventSelector 
{
public :
  explicit TTEventSelector(const edm::ParameterSet&, edm::ConsumesCollector&& iC );
  explicit TTEventSelector(const edm::ParameterSet&, edm::ConsumesCollector& iC );
  ~TTEventSelector() {}
  virtual int eventSelection(const edm::Event& iEvent, const edm::EventSetup& iSetup, TTree* tree, int sys) {
    std::cout<<"No event selection"<<std::endl; return 0;
  }
  const reco::Candidate* getLast(const reco::Candidate* p) const;
  virtual void setBranch(TTree* tree, int sys) {}
  virtual void resetBranch(){} 

  /*
  virtual float selectMuons(const cat::MuonCollection& muons, cat::MuonCollection& selmuons, TopEventCommonGlobal::sys_e sys) const =0;
  virtual float selectElecs(const cat::ElectronCollection& elecs, cat::ElectronCollection& selelecs, TopEventCommonGlobal::sys_e sys) const =0 ;
  */
  virtual cat::JetCollection selectJets(const cat::JetCollection& jets, const TopEventCommonGlobal::LeptonPtrs& recolep, TopEventCommonGlobal::sys_e sys);
  virtual cat::JetCollection selectBJets(const cat::JetCollection& jets) const;

  float getMuEffSF(const cat::Lepton& p, int sys) const
  {
    const int aid = abs(p.pdgId());
    if ( aid == 13 ) {
      const double pt = p.pt(), aeta = std::abs(p.eta());
      if      ( sys == +1 ) return muonSF_(pt, aeta,  1);
      else if ( sys == -1 ) return muonSF_(pt, aeta, -1);
      else return muonSF_(pt, aeta, 0);
    }
    return 1;
  }
  float getElEffSF(const cat::Lepton& p, int sys) const
  {
    const int aid = abs(p.pdgId());
    if ( aid == 11 ) {
      const auto& el = dynamic_cast<const cat::Electron&>(p);
      const double pt = p.pt(), aeta = std::abs(el.scEta());
      if      ( sys == +1 ) return elecSF_(pt, aeta,  1);
      else if ( sys == -1 ) return elecSF_(pt, aeta, -1);
      else return elecSF_(pt, aeta, 0);
    }
    return 1;
  }
   
protected:
  edm::EDGetTokenT<int> recoFiltersToken_, nGoodVertexToken_, lumiSelectionToken_;

  edm::EDGetTokenT<cat::MuonCollection>     muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection> elecToken_;
  edm::EDGetTokenT<cat::JetCollection>      jetToken_;
  edm::EDGetTokenT<cat::METCollection>      metToken_;
  edm::EDGetTokenT<reco::VertexCollection>   vtxToken_;

  float muonPtCut_, muonEtaCut_, muonIsoCut_;
  float electronPtCut_, electronEtaCut_, electronIsoCut_;
  float jetPtCut_, jetEtaCut_;
  
  std::string muonIDCut_, electronIDCut_, bJetCSV_;

  cat::ScaleFactorEvaluator muonSF_, elecSF_;

  cat::BTagWeightEvaluator csvWeight;
  cat::BTagWeightEvaluator bTagWeightL;
  cat::BTagWeightEvaluator bTagWeightM;
  cat::BTagWeightEvaluator bTagWeightT;
};
#endif
