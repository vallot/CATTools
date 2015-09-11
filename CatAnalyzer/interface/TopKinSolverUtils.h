#ifndef CATTools_CatAnalysis_TopKinSolverUtils_H
#include "Math/LorentzVector.h"
#include <iostream>
#include <cmath>
#include <vector>

struct KinSolverUtils {

  constexpr static double pi = std::asin(-1);
  constexpr static double mB = 4.8;
  constexpr static double mL = 0.0;
  constexpr static double mV = 0.0;

  //class TtFullLepSolution;
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LV;

  static inline bool isZero(const double x) { return std::abs(x) < 1e-5; };
  //void print(const std::vector<TtFullLepSolution>& sols);
  static void findCoeffs(const double mT, const double mW1, const double mW2,
                         const LV& l1, const LV& l2, const LV& b1, const LV& b2,
                         const double metX, const double metY,
                         std::vector<double>& koeff, std::vector<double>& cachedPars);
  static void getNuPxPyPzE(const double nuPx, const std::vector<double>& cachedPars,
                           double nu1sol[], double nu2sol[]);
  static void solve_quartic(const std::vector<double>& h, const double a4, const double b4,
                            std::vector<double>& v);
  static void solve_cubic(const double a, const double b, const double c, const double d,
                          std::vector<double>& v);
  static void solve_quadratic(const double a, const double b, const double c,
                              std::vector<double>& v);
  static void solve_linear(const double a, const double b,
                           std::vector<double>& v);

};

#endif

