#include "CATTools/CatAnalyzer/interface/TopKinSolverUtils.h"
#include <iostream>
#include <cmath>

using namespace std;

inline double sqr(const double x) { return x*x; };
inline double dxsqr(const double x, const double y) { return x*x - y*y; };

/*
void KinSolverUtils::print(const std::vector<TtFullLepSolution>& sols) {
  for ( int i=0, n=sols.size(); i<n; ++i ) {
    const auto& sol = sols[i];
    printf("\nSol: %d:  weight: %.3f dN: %.3f\n", i+1, sol.weight, sol.dN);
    printf(" top   pt=%.3f eta=%.3f phi=%.3f M=%.3f", sol.top().pt(), sol.top().eta(), sol.top().phi(), sol.top().mass());
    printf(" tbar  pt=%.3f eta=%.3f phi=%.3f M=%.3f", sol.tbar().pt(), sol.tbar().eta(), sol.tbar().phi(), sol.tbar().mass());
    printf(" nu    pt=%.3f eta=%.3f phi=%.3f M=%.3f", sol.nu().pt(), sol.nu().eta(), sol.nu().phi(), sol.nu().mass());
    printf(" nubar pt=%.3f eta=%.3f phi=%.3f M=%.3f", sol.nubar().pt(), sol.nubar().eta(), sol.nubar().phi(), sol.nubar().mass());
  }
}
*/

void KinSolverUtils::findCoeffs(const double mT, const double mW1, const double mW2,
                                const LV& l1, const LV& l2, const LV& j1, const LV& j2,
                                const double metX, const double metY,
                                std::vector<double>& kfs, std::vector<double>& cachedPars) {
  const double dmW1 = mW1*mW1-mL*mL-mV*mV;
  const double dmW2 = mW2*mW2-mL*mL-mV*mV;
  const double dmT = mT*mT-mB*mB-mL*mL-mV*mV;

  const double l1E = l1.energy(), j1E = j1.energy();
  const double jlEA = j1E + l1E;
  const double divA = 2*l1E*jlEA;
  const double a1 = (jlEA*dmW1-l1E*(dmT+2*j1E*l1E-2*(l1.Vect().Dot(j1.Vect()))))/divA;
  const double a2 = 2*(j1E*l1.px()-l1E*j1.px())/divA;
  const double a3 = 2*(j1E*l1.py()-l1E*j1.py())/divA;
  const double a4 = 2*(j1E*l1.pz()-l1E*j1.pz())/divA;

  const double l2E = l2.energy(), j2E = j2.energy();
  const double jlEB = j2E + l2E;
  const double divB = 2*l2E*jlEB;
  const double b1 = (jlEB*dmW2-l2E*(dmT+2*j2E*l2E-2*(l2.Vect().Dot(j2.Vect()))))/divB;
  const double b2 = 2*(j2E*l2.px()-l2E*j2.px())/divB;
  const double b3 = 2*(j2E*l2.py()-l2E*j2.py())/divB;
  const double b4 = 2*(j2E*l2.pz()-l2E*j2.pz())/divB;

  const double divC = 4*jlEA*jlEA;
  const double a14 = a1/a4, a24 = a2/a4, a34 = a3/a4;
  const double c00 = -4*(dxsqr(l1E, l1.py()) + dxsqr(l1E, l1.pz())*a34*a34 + 2*l1.py()*l1.pz()*a34)/divC;
  const double c10 = -8*(dxsqr(l1E, l1.pz())*a24/a34 - l1.px()*l1.py() + l1.px()*l1.pz()*a34 + l1.py()*l1.pz()*a24)/divC;
  const double c20 = -4*(dxsqr(l1E, l1.px()) + dxsqr(l1E, l1.pz())*sqr(a2/a4) + 2*l1.px()*l1.pz()*a24)/divC; 
  const double c11 = 4*(dmW1*(l1.py()-l1.pz()*a34)-2*dxsqr(l1E, l1.pz())*a14*a34-2*l1.py()*l1.pz()*a14)/divC;
  const double c21 = 4*(dmW1*(l1.px()-l1.pz()*a24)-2*dxsqr(l1E, l1.pz())*a14*a24-2*l1.px()*l1.pz()*a14)/divC;
  const double c22 = (dmW1*dmW1-4*dxsqr(l1E, l1.pz())*a14/a14-4*dmW1*l1.pz()*a14)/divC;

  const double divD = 4*jlEB*jlEB;
  const double b14 = b1/b4, b24 = b2/b4, b34 = b3/b4;
  const double D00 = -4*(dxsqr(l2E, l2.py()) + dxsqr(l2E, l2.pz())*b34*b34 + l2.py()*l2.pz()*b34)/divD;
  const double D10 = -8*(dxsqr(l2E, l2.pz())*b24*b34 - l2.px()*l2.py() + l2.px()*l2.pz()*b34 + l2.py()*l2.pz()*b24)/divD;
  const double D20 = -4*(dxsqr(l2E, l2.px()) + dxsqr(l2E, l2.pz())*b24*b24 + 2*l2.px()*l2.pz()*b24)/divD;
  const double D11 = 4*(dmW2*(l2.py()-l2.pz()*b34)-2*dxsqr(l2E, l2.pz())*b14*b34-2*l2.py()*l2.pz()*b14)/divD;
  const double D21 = 4*(dmW2*(l2.px()-l2.pz()*b24)-2*dxsqr(l2E, l2.pz())*b14*b24-2*l2.px()*l2.pz()*b14)/divD;
  const double D22 = (dmW2*dmW2-4*dxsqr(l2E, l2.pz())*b14*b14-4*dmW2*l2.pz()*b14)/divD; 

  const double d22 = D22+sqr(metX)*D20+sqr(metY)*D00+metX*metY*D10+metX*D21+metY*D11;
  const double d21 = -D21-2*metX*D20-metY*D10;
  const double d20 = D20;
  const double d11 = -D11-2*metY*D00-metX*D10;
  const double d10 = D10;
  const double d00  = D00;

  kfs.resize(5);

  kfs[4] = sqr(c00*d22) + c11*d22*(c11*d00-c00*d11)+c00*c22*(sqr(d11)-2*d00*d22)+c22*d00*(c22*d00-c11*d11);
  kfs[3] = c00*d21*(2*c00*d22-c11*d11)+c00*d11*(2*c22*d10+c21*d11)+c22*d00*(2*c21*d00-c11*d10)-c00*d22*(c11*d10+c10*d11)-2*c00*d00*(c22*d21+c21*d22)-d00*d11*(c11*c21+c10*c22)+c11*d00*(c11*d21+2*c10*d22);
  kfs[2] = sqr(c00)*(2*d22*d20+sqr(d21))-c00*d21*(c11*d10+c10*d11)+c11*d20*(c11*d00-c00*d11)+c00*d10*(c22*d10-c10*d22)+c00*d11*(2*c21*d10+c20*d11)+(2*c22*c20+sqr(c21))*sqr(d00)-2*c00*d00*(c22*d20+c21*d21+c20*d22)+c10*d00*(2*c11*d21+c10*d22)-d00*d10*(c11*c21+c10*c22)-d00*d11*(c11*c20+c10*c21);
  kfs[1] = c00*d21*(2*c00*d20-c10*d10)-c00*d20*(c11*d10+c10*d11)+c00*d10*(c21*d10+2*c20*d11)-2*c00*d00*(c21*d20+c20*d21)+c10*d00*(2*c11*d20+c10*d21)+c20*d00*(2*c21*d00-c10*d11)-d00*d10*(c11*c20+c10*c21);
  kfs[0] = sqr(c00*d20)+c10*d20*(c10*d00-c00*d10)+c20*d10*(c00*d10-c10*d00)+c20*d00*(c20*d00-2*c00*d20);

  // Append additional coefficients
  cachedPars = {
    d00, d11, d22, d10, d21, d20,
    c00, c11, c22, c10, c21, c20,
    a14, a24, a34, b14, b24, b34,
    metX, metY
  };
}

void KinSolverUtils::getNuPxPyPzE(const double px, const std::vector<double>& p,
                                  double nu1sol[], double nu2sol[]) {
  // See cachedPars 5 lines above for parameter ordering
  if ( p.size() < 20 ) return;

  const double metX = p[18], metY = p[19];

  const double d0 = p[0];
  const double d1 = p[1] + p[3]*px;
  const double d2 = p[2] + p[4]*px + p[5]*px*px;

  const double c0 = p[6];
  const double c1 = p[7] + p[9]*px;
  const double c2 = p[8] + p[10]*px + p[11]*px*px;

  nu1sol[0] = px;
  nu1sol[1] = (c0*d2-c2*d0)/(c1*d0-c0*d1);
  nu2sol[0] = metX - nu1sol[0];
  nu2sol[1] = metY - nu1sol[1];

  nu1sol[2] = -(p[12] + p[13]*nu1sol[0] + p[14]*nu1sol[1]);
  nu2sol[2] = -(p[15] + p[16]*nu2sol[0] + p[17]*nu2sol[2]);

  nu1sol[3] = sqrt(nu1sol[0]*nu1sol[0] + nu1sol[1]*nu1sol[1] + nu1sol[2]*nu1sol[2] + mV*mV);
  nu2sol[3] = sqrt(nu2sol[0]*nu2sol[0] + nu2sol[1]*nu2sol[1] + nu2sol[2]*nu2sol[2] + mV*mV);
}

void KinSolverUtils::solve_linear(const double a, const double b,
                                  std::vector<double>& v) {
  v.clear();
  if ( a == 0 ) return;

  v.push_back(-1*b/a);
}

void KinSolverUtils::solve_quadratic(const double a, const double b, const double c,
                                     std::vector<double>& v) {
  if ( isZero(a) ) return solve_linear(b, c, v);

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

void KinSolverUtils::solve_cubic(const double a, const double b, const double c, const double d,
                                 std::vector<double>& v) {
  if ( isZero(a) ) return solve_quadratic(b, c, d, v);

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
    v.push_back(-2*sqrt(abs(q))*cos((F+2*pi)/3)-s1/3);
    v.push_back(-2*sqrt(abs(q))*cos((F-2*pi)/3)-s1/3);
  }
  else {
    const long double atmp = r+sqrt(abs(r2-q3));
    const int sgnTmp = atmp >= 0 ? 1 : -1;
    const double A = -sgnTmp*pow(abs(atmp), 1./3);
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

void KinSolverUtils::solve_quartic(const std::vector<double>& h,
                                   const double a4, const double b4,
                                   std::vector<double>& v) {
  v.clear();
  if ( h.size() < 5 ) return;
  if ( isZero(a4) or isZero(b4) ) return;

  const double h0 = h[0], h1 = h[1], h2 = h[2], h3 = h[3], h4 = h[4];
  if ( isZero(h0) ) return solve_cubic(h1, h2, h3, h4, v);
  if ( isZero(h3) ) {
    solve_cubic(h0, h1, h2, h3, v);
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
    solve_cubic(1., 0., K1, K2, v);
    for ( auto& x : v ) x -= H1/4;
    v.push_back(-H1/4);
    return;
  }

  std::vector<double> v_t12;
  std::vector<std::pair<double, double> > vPair_tmp;
  solve_cubic(1., 2.*K1, (K1*K1-4*K3), -K2*K2, v_t12);
  for ( auto& x : v_t12 ) {
    if ( x < 0 ) continue;
    const double sx = sqrt(x);

    vPair_tmp.push_back(make_pair( sx, (K1+x-K2/sx)/2));
    vPair_tmp.push_back(make_pair(-sx, (K1+x+K2/sx)/2));
  }
  for ( const auto& vp : vPair_tmp ) {
    std::vector<double> pre_v1;
    solve_quadratic(1., vp.first, vp.second, pre_v1);
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
