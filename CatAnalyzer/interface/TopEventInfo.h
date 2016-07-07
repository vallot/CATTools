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
private :
  TopEventInfo(); 
  ~TopEventInfo(){};
public:
  void resetBr(){
    filtered = 0;
    tri = tri_up = tri_dn = 0;

    channel = 0;
    step = -1;
    step1 = step2 = step3 = step4 = step5 = step6 = step7 = step8 = 0;

    nvertex = 0;
    njet = nbjet = 0;
    met = -9;

    weight = puweight = puweight_up = puweight_dn = genweight = 1;
    mueffweight = mueffweight_up = mueffweight_dn = 1;
    eleffweight = eleffweight_up = eleffweight_dn = 1;
    btagweight = btagweight_up = btagweight_dn = 1;
    topPtWeight = 1.;

    csvweights.clear();
    pdfWeights.clear();
    scaleWeights_up.clear(); scaleWeights_dn.clear();


    gen_partonChannel = 0; gen_partonMode1 = 0; gen_partonMode2 = 0; gen_partonMode = 0;
    gen_partonInPhase = 0; gen_partonInPhaseLep = 0; gen_partonInPhaseJet = 0;
    gen_partonlep1 = TLorentzVector(); gen_partonlep1_pid = 0;
    gen_partonlep2 = TLorentzVector(); gen_partonlep2_pid = 0;
    gen_partontop1 = TLorentzVector(); gen_partontop2 = TLorentzVector();
    gen_partonjet1 = TLorentzVector(); gen_partonjet2 = TLorentzVector();
    gen_partonnu1 = TLorentzVector(); gen_partonnu2 = TLorentzVector();
    gen_partondilep = TLorentzVector(); gen_partonttbar = TLorentzVector();
    gen_partonttbar_dphi = 0;

    gen_pseudoChannel = 0;
    gen_pseudoInPhase = 0;
    gen_pseudolep1 = TLorentzVector(); gen_pseudolep1_pid = 0;
    gen_pseudolep2 = TLorentzVector(); gen_pseudolep2_pid = 0;
    gen_pseudotop1 = TLorentzVector(); gen_pseudotop2 = TLorentzVector();
    gen_pseudojet1 = TLorentzVector(); gen_pseudojet2 = TLorentzVector();
    gen_pseudonu1 = TLorentzVector(); gen_pseudonu2 = TLorentzVector();
    gen_pseudodilep = TLorentzVector(); gen_pseudottbar = TLorentzVector();
    gen_pseudottbar_dphi = 0;

    lep1 = TLorentzVector(); lep1_pid = 0;
    lep2 = TLorentzVector(); lep2_pid = 0;

    pseudonu1 = TLorentzVector(); pseudonu2 = TLorentzVector();
    pseudojet1 = TLorentzVector(); pseudojet2 = TLorentzVector();
    pseudotop1 = TLorentzVector(); pseudotop2 = TLorentzVector();
    pseudojet1_CSVInclV2 = 0; pseudojet2_CSVInclV2 = 0;
    pseudottbar = TLorentzVector();
    pseudottbar_dphi = 0;

    partonnu1 = TLorentzVector(); partonnu2 = TLorentzVector();
    partonjet1 = TLorentzVector(); partonjet2 = TLorentzVector();
    partontop1 = TLorentzVector(); partontop2 = TLorentzVector();
    partonjet1_CSVInclV2 = 0; partonjet2_CSVInclV2 = 0;
    partonttbar = TLorentzVector();
    partonttbar_dphi = 0;

    is3lep = -9;
    desyjet1 = TLorentzVector(); desyjet2 = TLorentzVector();
    desytop1 = TLorentzVector(); desytop2 = TLorentzVector();
    desyjet1_CSVInclV2 = 0; desyjet2_CSVInclV2 = 0;
    desyttbar = TLorentzVector();
    desyttbar_dphi = 0;

  }
  static TopEventInfo& getInstance() {
    static TopEventInfo ins;
    return ins;
  }
    
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

  //cat::ScaleFactorEvaluator muonSF_, elecSF_;

  /*
  cat::BTagWeightEvaluator csvWeight;
  cat::BTagWeightEvaluator bTagWeightL;
  cat::BTagWeightEvaluator bTagWeightM;
  cat::BTagWeightEvaluator bTagWeightT;
  */

  //std::unique_ptr<TtFullLepKinSolver> solver;
  //std::unique_ptr<cat::KinematicSolver> solver_;
  //std::unique_ptr<cat::KinematicSolver> solverPT_;
  //
  std::unique_ptr<KinematicReconstruction> kinematicReconstruction;

  int is3lep;
  TLorentzVector desyjet1, desyjet2, desytop1, desytop2, desyttbar;
  float desyjet1_CSVInclV2, desyjet2_CSVInclV2, desyttbar_dphi;

};

#endif
