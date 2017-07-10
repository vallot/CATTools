#ifndef CATTools_CatAnalyzer_KinematicFitters_H
#define CATTools_CatAnalyzer_KinematicFitters_H

#include "TLorentzVector.h"
#include "TMinuit.h"
#include <memory>

namespace cat {

struct KinematicFitter
{
  static double jetEResolution(double energy)
  {
    return std::hypot(0.05, 1.5/std::sqrt(energy));
  }

  static double metPhiResolution(double met)
  {
    return  0.05539 - 0.5183*exp(-0.01507*met);
  }

  static double metResolution(double met) { return 15.0; }

  static double twoObjectMassResolution(TLorentzVector &j1, double releres1, TLorentzVector &j2, double releres2)
  {
    // crude, but OK
    float massnominal = (j1+j2).M();
    TLorentzVector j1smeared = (1.0+releres1)*j1;
    TLorentzVector j2smeared = (1.0+releres2)*j2;

    const float deltamass1up = (j1smeared+j2).M()-massnominal;
    const float deltamass2up = (j1+j2smeared).M()-massnominal;
    return std::hypot(deltamass1up, deltamass2up)/massnominal;
  }

  static double twoJetMassResolution(TLorentzVector &j1, TLorentzVector &j2)
  {
    const float releres1 = KinematicFitter::jetEResolution(j1.E());
    const float releres2 = KinematicFitter::jetEResolution(j2.E());

    return KinematicFitter::twoObjectMassResolution(j1, releres1, j2, releres2);
  }

  std::unique_ptr<TMinuit> tm_;
};

};

#endif
