#ifndef __dileptonCommon__
#define __dileptonCommon__

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

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

namespace dileptonCommonGlobal {
  enum sys_e {sys_nom,
    sys_jes_u, sys_jes_d, sys_jer_u, sys_jer_d,
    sys_mu_u, sys_mu_d, sys_el_u, sys_el_d,
	      //sys_mueff_u, sys_mueff_d, sys_eleff_u, sys_eleff_d,
	      //sys_btag_u, sys_btag_d,
    nsys_e
  };
  const int NCutflow = 12;

  const std::string sys_name[nsys_e] = {
    "nom",
    "jes_u", "jes_d", "jer_u", "jer_d",
    "mu_u", "mu_d", "el_u", "el_d",
    //    "mueff_u", "mueff_d", "eleff_u", "eleff_d",
    //    "btag_u", "btag_d"
  };
  typedef std::vector<const cat::Lepton*> LeptonPtrs;
  typedef math::XYZTLorentzVector LV;
  typedef std::vector<LV> VLV;
}

class dileptonCommon : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchLuminosityBlocks> {
public:
  explicit dileptonCommon(const edm::ParameterSet&);
  ~dileptonCommon();
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void analyzeCustom(const edm::Event&, const edm::EventSetup&, int sys ) ;
  void genInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual int eventSelection(const edm::Event& iEvent, const edm::EventSetup& iSetup, int sys);
  void parameterInit(const edm::ParameterSet& iConfig);
  virtual void resetBr();
  void resetBrCommon();
  virtual void resetBrCustom();

  virtual void setBranch(TTree* tree, int sys);
  void setBranchCommon(TTree* tree, int sys);
  virtual void setBranchCustom(TTree* tree, int sys);

  virtual float selectMuons(const cat::MuonCollection& muons, cat::MuonCollection& selmuons, dileptonCommonGlobal::sys_e sys) const;
  virtual float selectElecs(const cat::ElectronCollection& elecs, cat::ElectronCollection& selelecs, dileptonCommonGlobal::sys_e sys) const;
  virtual cat::JetCollection selectJets(const cat::JetCollection& jets, const dileptonCommonGlobal::LeptonPtrs& recolep, dileptonCommonGlobal::sys_e sys);
  virtual cat::JetCollection selectBJets(const cat::JetCollection& jets) const;
  const reco::Candidate* getLast(const reco::Candidate* p) const;
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

// Use protect keyword for branches.
protected : 
  int b_run, b_lumi, b_event;
  int b_nvertex, b_step, b_channel, b_njet, b_nbjet;
  bool b_step1, b_step2, b_step3, b_step4, b_step5, b_step6, b_step7, b_step8, b_filtered, b_keepTtbarSignal;
  float b_tri, b_tri_up, b_tri_dn;
  float b_met, b_weight, b_puweight, b_puweight_up, b_puweight_dn, b_genweight,
    b_mueffweight, b_mueffweight_up, b_mueffweight_dn,
    b_eleffweight, b_eleffweight_up, b_eleffweight_dn,
    b_btagweight, b_btagweight_up, b_btagweight_dn;
  float b_topPtWeight;
  std::vector<float> b_pdfWeights, b_scaleWeights_up, b_scaleWeights_dn, b_csvweights;
  
  int b_is3lep;
  
  int b_gen_partonChannel, b_gen_partonMode1, b_gen_partonMode2, b_gen_partonMode;
  bool b_gen_partonInPhase, b_gen_partonInPhaseLep, b_gen_partonInPhaseJet;
  TLorentzVector b_gen_partonlep1; int b_gen_partonlep1_pid;
  TLorentzVector b_gen_partonlep2; int b_gen_partonlep2_pid;
  TLorentzVector b_gen_partontop1, b_gen_partontop2, b_gen_partonnu1, b_gen_partonnu2, b_gen_partonjet1, b_gen_partonjet2, b_gen_partondilep, b_gen_partonttbar;
  float b_gen_partonttbar_dphi;
  
  int b_gen_pseudoChannel;
  bool b_gen_pseudoInPhase;
  TLorentzVector b_gen_pseudolep1; int b_gen_pseudolep1_pid;
  TLorentzVector b_gen_pseudolep2; int b_gen_pseudolep2_pid;
  TLorentzVector b_gen_pseudotop1, b_gen_pseudotop2, b_gen_pseudonu1, b_gen_pseudonu2, b_gen_pseudojet1, b_gen_pseudojet2, b_gen_pseudodilep, b_gen_pseudottbar;
  float b_gen_pseudottbar_dphi;

  TLorentzVector b_lep1, b_lep2, b_dilep;
  int b_lep1_pid, b_lep2_pid;
  
  TLorentzVector b_partonnu1, b_partonnu2, b_partonjet1, b_partonjet2, b_partontop1, b_partontop2, b_partonttbar;
  float b_partonttbar_dphi;
  float b_partonjet1_CSVInclV2, b_partonjet2_CSVInclV2;

  TLorentzVector b_pseudonu1, b_pseudonu2, b_pseudojet1, b_pseudojet2, b_pseudotop1, b_pseudotop2, b_pseudottbar;
  float b_pseudottbar_dphi;
  float b_pseudojet1_CSVInclV2, b_pseudojet2_CSVInclV2;
  
  TLorentzVector b_desyjet1, b_desyjet2, b_desytop1, b_desytop2, b_desyttbar;
  float b_desyttbar_dphi;
  float b_desyjet1_CSVInclV2, b_desyjet2_CSVInclV2;
  // I/O variables.
  std::vector<TTree*> ttree_;
  TH1D * h_nevents;
  // Exception for easy coding.
  std::vector<std::vector<int> > cutflow_;
  cat::JetCollection selectedJets ;
  cat::JetCollection selectedBJets;
  std::vector<const cat::Lepton*> recolep_;
  dileptonCommonGlobal::LV met;
  std::unique_ptr<cat::KinematicSolver> solver_;
  std::unique_ptr<cat::KinematicSolver> solverPT_;

  edm::EDGetTokenT<int> recoFiltersToken_, nGoodVertexToken_, lumiSelectionToken_;
  edm::EDGetTokenT<float> genWeightToken_;
  edm::EDGetTokenT<std::vector<float>> pdfweightToken_, scaleupweightsToken_, scaledownweightsToken_;
  edm::EDGetTokenT<float> puweightToken_, puweightToken_up_, puweightToken_dn_, topPtWeight_;
  edm::EDGetTokenT<int> trigTokenMUEL_, trigTokenMUMU_, trigTokenELEL_;

  edm::EDGetTokenT<cat::MuonCollection>     muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection> elecToken_;
  edm::EDGetTokenT<cat::JetCollection>      jetToken_;
  edm::EDGetTokenT<cat::METCollection>      metToken_;
  edm::EDGetTokenT<reco::VertexCollection>   vtxToken_;
  edm::EDGetTokenT<int>          partonTop_channel_;
  edm::EDGetTokenT<std::vector<int> > partonTop_modes_;
  edm::EDGetTokenT<reco::GenParticleCollection> partonTop_genParticles_;
  edm::EDGetTokenT<edm::View<reco::Candidate> > pseudoTop_leptons_, pseudoTop_neutrinos_, pseudoTop_jets_;

private:
  void beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&) final;
  void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override {};

  cat::ScaleFactorEvaluator muonSF_, elecSF_;

  cat::BTagWeightEvaluator csvWeight;
  cat::BTagWeightEvaluator bTagWeightL;
  cat::BTagWeightEvaluator bTagWeightM;
  cat::BTagWeightEvaluator bTagWeightT;

  //std::unique_ptr<TtFullLepKinSolver> solver;
  std::unique_ptr<KinematicReconstruction> kinematicReconstruction;


};


#endif
