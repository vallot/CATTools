#ifndef CATTools_CatAnalysis_TopKinSolverUtils_H
#include <iostream>
#include <cmath>
#include <vector>

struct KinSolverUtils {

  constexpr static double pi = std::asin(-1);
  constexpr static double mB = 4.8;
  constexpr static double mL = 0.0;
  constexpr static double mV = 0.0;

  //class TtFullLepSolution;

  static inline bool isZero(const double x) { return std::abs(x) < 1e-5; };
  //void print(const std::vector<TtFullLepSolution>& sols);
  static void findCoeffs(const double mTop, const double mW,
                         const LorentzVector& l1, const LorentzVector& l2,
                         const LorentzVector& b1, const LorentzVector& b2,
                         double* kfs);
  static void solve_quartic(const double h0, const double h1, const double h2, const double h3, const double h4,
      const double a4, const double b4,
      std::vector<double>& v);
  static void solve_cubic(const double a, const double b, const double c, const double d,
      std::vector<double>& v);
  static void solve_quadratic(const double a, const double b, const double c,
      std::vector<double>& v);
  static void solve_linear(const double a, const double b,
      std::vector<double>& v);

};

#endif

