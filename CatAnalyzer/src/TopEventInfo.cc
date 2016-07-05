#include "CATTools/CatAnalyzer/interface/TopEventInfo.h"

void TopEventInfo::resetBr(){
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

  /*
  is3lep = -9;
  desyjet1 = TLorentzVector(); desyjet2 = TLorentzVector();
  desytop1 = TLorentzVector(); desytop2 = TLorentzVector();
  desyjet1_CSVInclV2 = 0; desyjet2_CSVInclV2 = 0;
  desyttbar = TLorentzVector();
  desyttbar_dphi = 0;
  */

}
