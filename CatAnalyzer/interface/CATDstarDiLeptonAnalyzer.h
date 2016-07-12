#ifndef __CATTools_CATDstarDiLeptonAnalyzer__
#define __CATTools_CATDstarDiLeptonAnalyzer__

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TTree.h"
#include "TH1D.h"

#include "CATTools/CatAnalyzer/interface/TTDiLeptonEventSelector.h"
#include "CATTools/CatAnalyzer/interface/TopEventCommon.h"
#include "CATTools/DataFormats/interface/SecVertex.h"
#include "TClonesArray.h"


class CATDstarDiLeptonAnalyzer : public TopEventCommon {
public:
  explicit CATDstarDiLeptonAnalyzer(const edm::ParameterSet&);
  virtual ~CATDstarDiLeptonAnalyzer() { showSummary(); }
  void analyzeCustom(const edm::Event&, const edm::EventSetup&, int sys ) final;
  virtual void setBranchCustom(TTree* , int ) final;
  virtual void resetBranchCustom() final;
  virtual void showSummary() final;
  virtual void endJob() final;
  virtual void setEventSelection(const edm::ParameterSet& iConfig) override{  
    eventSelect_ = new TTDiLeptonEventSelector(iConfig, consumesCollector()) ;
    for (int sys = 0; sys < TopEventCommonGlobal::nsys_e; ++sys){
      auto tr = ttree_[sys];
      eventSelect_->setBranch(tr, sys);
    }
  }
  std::shared_ptr<TLorentzVector> mcMatching( std::vector<TLorentzVector>& aGens, TLorentzVector& aReco) ;
// Use protect keyword for branches.
protected : 
    std::vector<bool> b_d0_true, b_d0_fit;
    //std::vector<float> b_d0_pt, b_d0_eta, b_d0_phi, b_d0_m;
    std::vector<float> b_d0_LXY, b_d0_L3D, b_d0_dRTrue, b_d0_relPtTrue, b_d0_dca;
    std::vector<float> b_d0_dau1_q;
    std::vector<float> b_d0_dau2_q;

    std::vector<bool> b_dstar_true, b_dstar_fit;
    //std::vector<float> b_dstar_pt, b_dstar_eta, b_dstar_phi, b_dstar_m;
    std::vector<float> b_dstar_q, b_dstar_LXY, b_dstar_L3D, b_dstar_dRTrue, b_dstar_relPtTrue, b_dstar_dca, b_dstar_dca2, b_dstar_dca3;
    std::vector<float> b_dstar_dau1_q;
    std::vector<float> b_dstar_dau2_q;
    std::vector<float> b_dstar_dau3_q;

    TClonesArray *b_d0,    *b_d0_dau1,    *b_d0_dau2; 
    TClonesArray *b_dstar, *b_dstar_dau1, *b_dstar_dau2, *b_dstar_dau3; 

    edm::EDGetTokenT<cat::SecVertexCollection>      d0Token_;
    edm::EDGetTokenT<cat::SecVertexCollection>      dstarToken_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> >            mcSrc_;

    double matchingDeltaR_;

  

private:

};


#endif
