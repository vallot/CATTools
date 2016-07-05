#ifndef __CATTools_TopEventInfo__
#define __CATTools_TopEventInfo__

#include "CATTools/DataFormats/interface/Lepton.h"
#include "CATTools/CommonTools/interface/TTbarModeDefs.h"
#include "CATTools/CommonTools/interface/ScaleFactorEvaluator.h"
#include "CATTools/CatAnalyzer/interface/BTagWeightEvaluator.h"
#include "CATTools/CatAnalyzer/interface/KinematicSolvers.h"

#include "CATTools/CommonTools/interface/AnalysisHelper.h"
#include "CATTools/CatAnalyzer/interface/analysisUtils.h"
#include "CATTools/CatAnalyzer/interface/KinematicReconstruction.h"
#include "TTree.h"
#include "TH1D.h"


class TopEventInfo {
public:
  /*
  std::unordered_map<std::string, int  > mapInt;
  std::unordered_map<std::string, float> mapFloat;
  std::unordered_map<std::string, bool> mapBool;
  std::unordered_map<std::string, std::vector<float> > mapVFloat;
  std::unordered_map<std::string, TLorentzVector > mapVFloat;
  */
  TopEventInfo() {
    resetBr();
  }
  void resetBr();
  int run, lumi, event;
  int nvertex, step, channel, njet, nbjet;

  // EventSelection Part
  bool step1, step2, step3, step4, step5, step6, step7, step8, filtered, keepTtbarSignal;
  // Trigger Part
  float tri, tri_up, tri_dn;
  // Weight Part
  float met, weight, puweight, puweight_up, puweight_dn, genweight,
    mueffweight, mueffweight_up, mueffweight_dn,
    eleffweight, eleffweight_up, eleffweight_dn,
    btagweight, btagweight_up, btagweight_dn;
  float topPtWeight;
  std::vector<float> pdfWeights, scaleWeights_up, scaleWeights_dn, csvweights;

  // Gen : Parton Level  
  int gen_partonChannel, gen_partonMode1, gen_partonMode2, gen_partonMode;
  bool gen_partonInPhase, gen_partonInPhaseLep, gen_partonInPhaseJet;
  TLorentzVector gen_partonlep1; int gen_partonlep1_pid;
  TLorentzVector gen_partonlep2; int gen_partonlep2_pid;
  TLorentzVector gen_partontop1, gen_partontop2, gen_partonnu1, gen_partonnu2, gen_partonjet1, gen_partonjet2, gen_partondilep, gen_partonttbar;
  float gen_partonttbar_dphi;
  
  // Gen : Particle Level  
  int gen_pseudoChannel;
  bool gen_pseudoInPhase;
  TLorentzVector gen_pseudolep1; int gen_pseudolep1_pid;
  TLorentzVector gen_pseudolep2; int gen_pseudolep2_pid;
  TLorentzVector gen_pseudotop1, gen_pseudotop2, gen_pseudonu1, gen_pseudonu2, gen_pseudojet1, gen_pseudojet2, gen_pseudodilep, gen_pseudottbar;
  float gen_pseudottbar_dphi;

  TLorentzVector lep1, lep2, dilep;
  int lep1_pid, lep2_pid;
  
  TLorentzVector partonnu1, partonnu2, partonjet1, partonjet2, partontop1, partontop2, partonttbar;
  float partonttbar_dphi;
  float partonjet1_CSVInclV2, partonjet2_CSVInclV2;

  TLorentzVector pseudonu1, pseudonu2, pseudojet1, pseudojet2, pseudotop1, pseudotop2, pseudottbar;
  float pseudottbar_dphi;
  float pseudojet1_CSVInclV2, pseudojet2_CSVInclV2;


  std::vector<std::vector<int> > cutflow_;
  cat::JetCollection selectedJets ;
  cat::JetCollection selectedBJets;
  std::vector<const cat::Lepton*> recolep_;
  math::XYZTLorentzVector metlv;

  cat::ScaleFactorEvaluator muonSF_, elecSF_;

  cat::BTagWeightEvaluator csvWeight;
  cat::BTagWeightEvaluator bTagWeightL;
  cat::BTagWeightEvaluator bTagWeightM;
  cat::BTagWeightEvaluator bTagWeightT;

  //std::unique_ptr<TtFullLepKinSolver> solver;
  //std::unique_ptr<cat::KinematicSolver> solver_;
  //std::unique_ptr<cat::KinematicSolver> solverPT_;
  //
  std::unique_ptr<KinematicReconstruction> kinematicReconstruction;

};
#endif
