#ifndef CATTools_CatAnalyzer_LepJetsFitterFCNH
#define CATTools_CatAnalyzer_LepJetsFitterFCNH

#include "TLorentzVector.h"
#include "TMinuit.h"
#include <iostream>
#include "CATTools/DataFormats/interface/Jet.h"

namespace fcnh {
  double SolvettbarLepJets(double &nupz, Double_t &metscale, Double_t &blscale, Double_t &bjscale, Double_t &j1scale, Double_t &j2scale);
  //void   FindHadronicTop        (TLorentzVector &lepton, std::vector<cat::ComJet> &jets, TLorentzVector &met, bool usebtaginfo, bool useCSVOrderinfo, std::vector<int> &bestindices, float &bestchi2, TLorentzVector &nusol, TLorentzVector &blrefit, TLorentzVector &bjrefit, TLorentzVector &j1refit, TLorentzVector &j2refit);
  void   FindHadronicTop        (TLorentzVector &lepton, std::vector<cat::Jet>& jets, TLorentzVector &met, bool usebtaginfo, bool useCSVOrderinfo, std::vector<int> &bestindices, float &bestchi2, TLorentzVector &nusol, TLorentzVector &blrefit, TLorentzVector &bjrefit, TLorentzVector &j1refit, TLorentzVector &j2refit);

  void fcnfull(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

  void InitMinuit();

  extern TMinuit *tm;

  extern TLorentzVector tmplep, tmpnu, tmpbl, tmpbj, tmpj1, tmpj2;
  extern float blres, bjres, j1res, j2res, metres;

  extern const float CSVWP;

}

#endif
