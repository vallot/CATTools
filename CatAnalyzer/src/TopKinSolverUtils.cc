#include "CATTools/CatAnalyzer/interface/TopKinsolverUtils.h"
#include <iostream>
#include <cmath>

using namespace std;
using namespace KinSolver;

void print(const std::vector<TtFullLepSolution>& sols) {
  for ( int i=0, n=sols.size(); i<n; ++i ) {
    const auto& sol = sols[i];
    printf("\nSol: %d:  weight: %.3f dN: %.3f\n", i+1, sol.weight, sol.dN);
    printf(" top   pt=%.3f eta=%.3f phi=%.3f M=%.3f", sol.top().pt(), sol.top().eta(), sol.top().phi(), sol.top().mass());
    printf(" tbar  pt=%.3f eta=%.3f phi=%.3f M=%.3f", sol.tbar().pt(), sol.tbar().eta(), sol.tbar().phi(), sol.tbar().mass());
    printf(" nu    pt=%.3f eta=%.3f phi=%.3f M=%.3f", sol.nu().pt(), sol.nu().eta(), sol.nu().phi(), sol.nu().mass());
    printf(" nubar pt=%.3f eta=%.3f phi=%.3f M=%.3f", sol.nubar().pt(), sol.nubar().eta(), sol.nubar().phi(), sol.nubar().mass());
  }
}

void findCoeffs(const double* kfs) {

}

void eqn_linear(const double a, const double b,
                std::vector<double>& v) {
  v.clear();
  if ( a == 0 ) return;

  v.push_back(-1*b/a);
}

void eqn_quadratic(const double a, const double b, const double c,
                   std::vector<double>& v) {
  if ( isZero(a) ) return eqn_linear(b, c, v);

  v.clear();
  const double det = b*b - 4*a*c;
  if ( det < 0 ) return;

  if ( isZero(det) ) v.push_back(-b/2/a);
  else {
    const double rdet = sqrt(det);
    v.push_back((-b-rdet)/2/a);
    v.push_back((-b+rdet)/2/a);
  }
}

void eqn_cubic(const double a, const double b, const double c, const double d,
               std::vector<double>& v) {
  if ( isZero(a) ) return eqn_quadratic(b, c, d, v);

  const double s1 = b/a;
  const double s2 = c/a;
  const double s3 = d/a;

  const double q = (s1*s1-3*s2)/9;
  const double q3 = q*q*q;
  const double r = (2*s1*s1*s1-9*s1*s2+27*s3)/54;
  const double r2 = r*r;
  const double det = r2-q3;

  if ( det < 0 ) {
    const double F = acos(r/sqrt(q3));
    v.push_back(-2*sqrt(abs(q))*cos(F/3)-s1/3);
    v.push_back(-2*sqrt(abs(q))*cos((F+2*TMath::Pi())/3)-s1/3);
    v.push_back(-2*sqrt(abs(q))*cos((F-2*TMath::Pi())/3)-s1/3);
  }
  else {
    const long double atmp = r+sqrt(abs(r2-q3));
    const double A = -sign(atmp)*pow(abs(atmp), 1./3);
    const double B = isZero(A) ? 0 : q/A;
    if ( isZero(det) ) {
      v.push_back(A+B-s1/3);
      v.push_back(-0.5*(A+B)-s1/3);
    }
    else {
      v.push_back(A+B-s1/3);
    }
  }
}

void eqn_quartic(const double h0, const double h1, const double h2, const double h3, const double h4,
                 const double a4, const double b4,
                 std::vector<double>& v) {
  v.clear();
  if ( isZero(a4) or isZero(b4) ) return;

  if ( isZero(h0) ) return eqn_cubic(h1, h2, h3, h4, v);
  if ( isZero(h3) ) {
    eqn_cubic(h0, h1, h2, h3, v);
    v.push_back(0);
    return;
  }

  const double H1 = h1/h0;
  const double H2 = h2/h0;
  const double H3 = h3/h0;
  const double H4 = h4/h0;
  const double K1 = H2 - 3*H1*H1/8;
  const double K2 = H3 + H1*H1*H1/8 - H1*H2/2;
  const double K3 = H4 - 3*pow(H1, 4)/256 + H1*H1*H2/16 - H1*H3/4;

  if ( isZero(K3) ) {
    eqn_cubic(1, 0, K1, K2, v);
    for ( auto& x : v ) x -= H1/4;
    v.push_back(-H1/4);
    return;
  }

  std::vector<double> v_t12;
  std::vector<std::pair<double, double> > vPair_tmp;
  eqn_cubic(1, 2*K1, (K1*K1-4*K3), -K2*K2, v_t12);
  for ( auto& x : v_t12 ) {
    if ( x < 0 ) continue;
    const double sx = sqrt(x);

    vPair_tmp.push_back(make_pair( sx, (K1+x-K2/sx)/2));
    vPair_tmp.push_back(make_pair(-sx, (K1+x+K2/sx)/2));
  }
  for ( const auto& vp : vPair_tmp ) {
    std::vector<double> pre_v1;
    eqn_quadratic(1, vp.first, vp.second, pre_v1);
    for ( const auto xPre : pre_v1 ) {
      bool isOverlap = false;
      for ( const auto x : v ) {
        if ( abs(x - xPre) < 0.02 ) { isOverlap = true ; break; }
      }
      if ( !isOverlap ) v.push_back(xPre);
    }
  }
  for ( auto& x : v ) x -= H1/4;
}
