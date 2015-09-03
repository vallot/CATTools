#ifndef CATTools_CatAnalysis_TopKinSolverUtils_H
#include <iostream>

namespace KinSolver {

const double mB = 4.8;
const double mL = 0.0;
const double mV = 0.0;

class TtFullLepSolution;

void print(const std::vector<TtFullLepSolution>& sols);
void findCoeffs(const double* kfs);
void eqn_quartic(const double h0, const double h1, const double h2, const double h3, const double h4,
                 std::vector<double>& v);
void eqn_cubic(const double a, const double b, const double c, const double d,
               std::vector<double>& v);
void eqn_quadratic(const double a, const double b, const double c,
                   std::vector<double>& v);
void eqn_linear(const double a, const double b,
                std::vector<double>& v);

};

#endif

