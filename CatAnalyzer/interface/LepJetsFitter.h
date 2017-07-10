#ifndef CATTools_CatAnalyzer_LepJetsFitter
#define CATTools_CatAnalyzer_LepJetsFitter

#include "CATTools/CatAnalyzer/interface/KinematicFitters.h"
#include "TLorentzVector.h"
#include <iostream>

namespace cat {
  // Jet Class
  class  ComJet: public TLorentzVector{
  public:
    float CSV, CvsL, CvsB;
    int Flavour, pTIndex, Mom;
  };
}

namespace ttbb{

  double SolvettbarLepJets(double &nupz, Double_t &metscale, Double_t &blscale, Double_t &bjscale, Double_t &j1scale, Double_t &j2scale);
  void   FindHadronicTop        (TLorentzVector &lepton, std::vector<cat::ComJet> &jets, TLorentzVector &met, bool usebtaginfo, bool useCSVOrderinfo, std::vector<int> &bestindices, float &bestchi2, TLorentzVector &nusol, TLorentzVector &blrefit, TLorentzVector &bjrefit, TLorentzVector &j1refit, TLorentzVector &j2refit);

  void fcnfull(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

  void InitMinuit();

  extern TMinuit *tm;

  extern TLorentzVector tmplep, tmpnu, tmpbl, tmpbj, tmpj1, tmpj2;
  extern float blres, bjres, j1res, j2res, metres;

  extern const float CSVWP;

}

#endif
