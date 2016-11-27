#ifndef LepJets_Fitter_H
#define LepJets_Fitter_H

#include "TLorentzVector.h"
#include "TMinuit.h"
#include <iostream>

namespace cat {
  // Jet Class
  class  ComJet: public TLorentzVector{
  public:
    float CSV, CvsL, CvsB;
    int Flavour, pTIndex, Mom;
  };
}
double JetEResolution         (double energy);
double METPhiResolution       (double met);
double METResolution          (double met);
double TwoObjectMassResolution(TLorentzVector &j1, double releres1, TLorentzVector &j2, double releres2);
double TwoJetMassResolution   (TLorentzVector &j1, TLorentzVector &j2);// relative mass resolution
double SolvettbarLepJets(double &nupz, Double_t &metscale, Double_t &blscale, Double_t &bjscale, Double_t &j1scale, Double_t &j2scale);
void   FindHadronicTop        (TLorentzVector &lepton, std::vector<cat::ComJet> &jets, TLorentzVector &met, bool usebtaginfo, bool useCSVOrderinfo, std::vector<int> &bestindices, float &bestchi2, TLorentzVector &nusol, TLorentzVector &blrefit, TLorentzVector &bjrefit, TLorentzVector &j1refit, TLorentzVector &j2refit);

void fcnfull(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void fcn    (Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

#endif
