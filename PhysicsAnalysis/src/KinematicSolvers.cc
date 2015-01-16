#include "CATTools/PhysicsAnalysis/interface/KinematicSolvers.h"
#include "TopQuarkAnalysis/TopKinFitter/interface/TtFullLepKinSolver.h"
#include <boost/assign/std/vector.hpp> // for 'operator+=()'
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

using namespace cat;
using namespace boost::assign;
using namespace std;

void TTDileptonSolver::solve(const KinematicSolver::LorentzVector input[])
{
  quality_ = -1e9; // default value
  values_.clear();

  const auto& met = input[0];
  const auto& l1 = input[1];
  const auto& l2 = input[2];
  const auto& j1 = input[3];
  const auto& j2 = input[4];

  nu1_ = met/2;
  nu2_ = met/2;

  const double mW = 80.4;
  if ( abs((nu1_+l1).mass()-mW)+abs((nu2_+l2).mass()-mW) <
       abs((nu1_+l2).mass()-mW)+abs((nu2_+l1).mass()-mW) ) swap(nu1_, nu2_);
  //if ( abs((nu1_+l1).mass()-mW) > 20 or abs((nu2_+l2).mass()-mW) > 20 ) return;

  quality_ = abs((nu1_+l1+j1).mass()-172.5) + abs((nu2_+l2+j2).mass()-172.5);
  quality_ = 0.01/(0.01+quality_);

  //const auto& lb1 = l1+j1;
  //const auto& lb2 = l2+j2;
  //if ( lb1.mass() > 150 or lb2.mass() > 150 ) quality_ = 0;

  const auto& vsum = l1+l2+j1+j2+met;
  if ( vsum.mass() < 300 ) quality_ = 0;

  values_.push_back(vsum.mass());
}

namespace CATMT2
{

const double mbsqr = 4.8*4.8;
const double mwsqr = 80.4*80.4;

double fconstraint(double y, void* p)
{
  double* pars = (double*)p;

  const double px1 = pars[1], py1 = pars[2];
  const double px2 = pars[3], py2 = pars[4];
  const double metx = pars[5], mety = pars[6];
  const double e1 = sqrt(mbsqr + px1*px1 + py1*py1);
  const double e2 = sqrt(mbsqr + px2*px2 + py2*py2);

  const double x1 = pars[0], y1 = y;
  const double x2 = metx-x1, y2 = mety-y1;

  const double retVal = e1*sqrt(x1*x1 + y1*y1 + mwsqr)
                      - e2*sqrt(x2*x2 + y2*y2 + mwsqr)
                      - (px1+px2)*x1 - (py1+py2)*y1
                      + px2*metx + py2*mety;

  return retVal;
}

double solveConstraint(double qx1, double* p0)
{
  // Solve constraint equation to find y
  double pars[] = {qx1, p0[0], p0[1], p0[2], p0[3], p0[4], p0[5]};
  const double xmin = p0[6];
  const double xmax = p0[7];

  //gsl_root_fsolver* solver = gsl_root_fsolver_alloc(gsl_root_fsolver_bisection);
  gsl_root_fsolver* solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
  gsl_function fgsl_constraint;
  fgsl_constraint.function = &fconstraint;
  fgsl_constraint.params = &pars;

  // Check solution exists in this range and expand if needed
  double ymin = xmin, ymax = xmax;
  while ( fconstraint(ymin, pars)*fconstraint(ymax, pars) >= 0 )
  {
    const double dy = (ymax-ymin)/2;
    ymin -= dy;
    ymax += dy;
  }
  gsl_root_fsolver_set(solver, &fgsl_constraint, ymin, ymax);

  double qy1;
  for ( int iter = 0; iter < 100; ++iter )
  {
    gsl_root_fsolver_iterate(solver);
    qy1 = gsl_root_fsolver_root(solver);
    const double xlo = gsl_root_fsolver_x_lower(solver);
    const double xhi = gsl_root_fsolver_x_upper(solver);
    int status = gsl_root_test_interval(xlo, xhi, 0, 1e-3);
    if ( status != GSL_CONTINUE ) break;
  }
  gsl_root_fsolver_free(solver);

  return qy1;
}

double mtsqr(const double qx, const double qy,
             const double px, const double py)
{
  const double eb = sqrt(mbsqr + px*px + py*py);
  const double ew = sqrt(mwsqr + qx*qx + qy*qy);
  return mwsqr + mbsqr + 2*eb*ew - 2*px*qx - 2*py*qy;
}

double fminimize(double qx1, void* p)
{
  double* p0 = (double*)p;
  const double qy1 = solveConstraint(qx1, p0);

  // Return transverse mass squared
  const double px1 = p0[0], py1 = p0[1];
  return mtsqr(qx1, qy1, px1, py1);
}

}

void MT2Solver::solve(const KinematicSolver::LorentzVector input[])
{
  using namespace CATMT2;

  quality_ = -1e9; // default value
  values_.clear();

  const auto& met = input[0];
  const auto& l1 = input[1];
  const auto& l2 = input[2];
  const auto& j1 = input[3];
  const auto& j2 = input[4];

  // pars : px1, py1, px2, py2, metx, mety, xmin, xmax
  double pars[] = {j1.px(), j1.py(), j2.px(), j2.py(), met.px(), met.py(), -l1.pt(), l1.py()};
  double& xmin = pars[6];
  double& xmax = pars[7];

  // Do the minimization
  gsl_min_fminimizer* minimizer = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
  //gsl_min_fminimizer* minimizer = gsl_min_fminimizer_alloc(gsl_min_fminimizer_quad_golden);
  gsl_function fgsl_minimize;
  fgsl_minimize.function = &fminimize;
  fgsl_minimize.params = &pars;
  // Check minimum exists in this range and expand if needed
  while ( fminimize(xmin, pars) <= fminimize((xmin+xmax)/2, pars) or
          fminimize(xmax, pars) <= fminimize((xmin+xmax)/2, pars) )
  {
    const double dx = (xmax-xmin)/2;
    xmax += dx; xmin -= dx;
  }
  gsl_min_fminimizer_set(minimizer, &fgsl_minimize, (xmin+xmax)/2, xmin, xmax);

  double qx1 = 1e9;
  for ( int iter=0; iter<100; ++iter )
  {
    gsl_min_fminimizer_iterate(minimizer);
    qx1 = gsl_min_fminimizer_x_minimum(minimizer);
    const double xlo = gsl_min_fminimizer_x_lower(minimizer);
    const double xhi = gsl_min_fminimizer_x_upper(minimizer);

    int status = gsl_min_test_interval(xlo, xhi, 0.001, 0);
    if ( status != GSL_CONTINUE ) break;
  }
  gsl_min_fminimizer_free(minimizer);

  const double qy1 = solveConstraint(qx1, pars);
  const double qx2 = met.px()-qx1;
  const double qy2 = met.py()-qy1;
  const double mt2 = sqrt(mtsqr(qx1, qy1, pars[0], pars[1]));

  const double ew1 = sqrt(mwsqr + qx1*qx1 + qy1*qy1);
  const double ew2 = sqrt(mwsqr + qx2*qx2 + qy2*qy2);
  KinematicSolver::LorentzVector w1(qx1, qy1, 0, ew1);
  KinematicSolver::LorentzVector w2(qx2, qy2, 0, ew2);
  nu1_ = w1-l1;
  nu2_ = w2-l2;

  values_.push_back(mt2);

}

void MAOSSolver::solve(const KinematicSolver::LorentzVector input[])
{
  MT2Solver::solve(input); // run the MT2 solver at the first step.

  // recover neutrino z component by assuming MT2 value as the top mass

}

CMSKinSolver::CMSKinSolver()
{
  std::vector<double> nuPars;
  nuPars += 30.7137,56.2880,23.0744,59.1015,24.9145;
  solver_ = new TtFullLepKinSolver(100, 300, 1, nuPars);
}

CMSKinSolver::~CMSKinSolver()
{
  if ( solver_ ) delete solver_;
}

void CMSKinSolver::solve(const KinematicSolver::LorentzVector input[])
{
  quality_ = -1e9; // default value
  values_.clear();

  const double metX = input[0].px();
  const double metY = input[0].py();
  const TLorentzVector l1(input[1].px(), input[1].py(), input[1].pz(), input[1].e());
  const TLorentzVector l2(input[2].px(), input[2].py(), input[2].pz(), input[2].e());
  const TLorentzVector j1(input[3].px(), input[3].py(), input[3].pz(), input[3].e());
  const TLorentzVector j2(input[4].px(), input[4].py(), input[4].pz(), input[4].e());

  solver_->SetConstraints(metX+l1.Px()+j1.Px()+l2.Px()+j2.Px(),
                          metY+l1.Py()+j1.Py()+l2.Py()+j2.Py());
  auto nuSol = solver_->getNuSolution(l1, l2, j1, j2);
  quality_  = nuSol.weight;
  if ( quality_ < 0 ) return;

  nu1_ = nuSol.neutrino.p4();
  nu2_ = nuSol.neutrinoBar.p4();

  const LorentzVector t1 = input[1]+input[3]+nu1_;
  const LorentzVector t2 = input[2]+input[4]+nu2_;
  values_.push_back((t1+t2).mass());
  values_.push_back(t1.mass());
  values_.push_back(t2.mass());
}

void NuWeightSolver::solve(const KinematicSolver::LorentzVector input[])
{
//cout << "STARTING NUSOLVER" << endl;
  quality_ = -1e9;
  values_.clear();

  const double metX = input[0].px();
  const double metY = input[0].py();
  const auto& l1 = input[1];
  const auto& l2 = input[2];
  const auto& j1 = input[3];
  const auto& j2 = input[4];

  const double sigmaEXsqr = 1, sigmaEYsqr = 1;

  double bestWeight = -1e9;
  double bestPx1 = 1e9, bestPx2 = 1e9, bestPy1 = 1e9, bestPy2 = -1e9;
  double bestEta1 = -1e9, bestEta2 = -1e9;
  double tmpSolX1[2], tmpSolY1[2], tmpSolX2[2], tmpSolY2[2];
  for ( double mt = 100; mt < 300; mt += 1 )
  {
    for ( double nu1Eta = -5.0; nu1Eta <= 5.0; nu1Eta += 0.05 )
    {
      if ( !computeNuPxPy(l1, j1, mt, nu1Eta, tmpSolX1[0], tmpSolY1[0], tmpSolX1[1], tmpSolY1[1]) ) continue;
      for ( double nu2Eta = -5.0; nu2Eta <= 5.0; nu2Eta += 0.05 )
      {
        if ( !computeNuPxPy(l2, j2, mt, nu2Eta, tmpSolX2[0], tmpSolY2[0], tmpSolX2[1], tmpSolY2[1]) ) continue;
        for ( int i=0; i<2; ++i )
        {
          for ( int j=0; j<2; ++j )
          {
            const double dpx = metX-tmpSolX1[i]-tmpSolX2[j];
            const double dpy = metY-tmpSolY1[i]-tmpSolY2[j];
            const double weight = exp(-dpx*dpx/2/sigmaEXsqr)*exp(-dpy*dpy/2/sigmaEYsqr);
            if ( weight > bestWeight )
            {
              bestWeight = weight;
              bestPx1 = tmpSolX1[i];
              bestPy1 = tmpSolY1[i];
              bestPx2 = tmpSolX2[j];
              bestPy2 = tmpSolY2[j];
              bestEta1 = nu1Eta;
              bestEta2 = nu2Eta;
cout << "New mT=" << mt << " w=" << weight << " pT1=" << hypot(bestPx1, bestPy1) << " eta1=" << bestEta1 << " pT2=" << hypot(bestPx2, bestPy2) << " eta2=" << bestEta2 << endl;
            }
          }
        }
      }
    }
  }

  quality_ = bestWeight;
  const double bestPz1 = hypot(bestPx1, bestPy1)*sinh(bestEta1);
  const double bestPz2 = hypot(bestPx2, bestPy2)*sinh(bestEta2);
  nu1_ = LorentzVector(bestPx1, bestPy1, bestPz1, sqrt(bestPx1*bestPx1+bestPy1*bestPy1+bestPz1*bestPz1));
  nu2_ = LorentzVector(bestPx2, bestPy2, bestPz2, sqrt(bestPx2*bestPx2+bestPy2*bestPy2+bestPz2*bestPz2));

  const LorentzVector t1 = l1+j1+nu1_;
  const LorentzVector t2 = l2+j2+nu2_;
  values_.push_back((t1+t2).mass());
  values_.push_back(t1.mass());
  values_.push_back(t2.mass());
if ( quality_ > 0 )
{
cout << "BEST W=" << quality_ << ' ' << hypot(bestPx1, bestPy1) << ' ' << hypot(bestPx2, bestPy2) << endl;
cout << " Mtop=" << t1.mass() << ' ' << t2.mass() << endl;
cout << " MW  =" << (l1+nu1_).mass() << ' ' << (l2+nu2_).mass() << endl;
}
}

bool NuWeightSolver::computeNuPxPy(const KinematicSolver::LorentzVector& lep,
                                   const KinematicSolver::LorentzVector& jet,
                                   const double mT, const double nuEta,
                                   double& nuPx1, double& nuPy1, double& nuPx2, double& nuPy2) const
{
  const double mTsqr = mT*mT;
  const double mWsqr = 80.4*80.4;
  const double mBsqr = 4.8*4.8;
  const double mLsqr = 0;

  const double alpha1 = (mTsqr - mBsqr - mWsqr)/2; // eqn B.6
  const double alpha2 = (mWsqr - mLsqr)/2; // eqn B.7

  const double pxl = lep.px(), pyl = lep.py(), pzl = lep.pz(), el = lep.energy();
  const double pxb = jet.px(), pyb = jet.py(), pzb = jet.pz(), eb = jet.energy();

  const double denB = eb*sinh(nuEta) - pzb*cosh(nuEta); // eqn B.14
  const double denL = el*sinh(nuEta) - pzl*cosh(nuEta); // eqn B.14

  const double ab = pxb/denB; // eqn B.14
  const double bb = pyb/denB; // eqn B.14
  const double cb = alpha1/denB; // eqn B.14

  const double al = pxl/denL; // eqn B.14
  const double bl = pyl/denL; // eqn B.14
  const double cl = alpha2/denL; // eqn B.14

  const double kpa = (bl-bb)/(ab-al); // eqn B.16
  const double eps = (cl-cb)/(ab-al); // eqn B.16

  // Now solve the equation B.17
  // 0 = (k^2 (a^2-1)+b^2-1) x^2+x (2 e k (a^2-1)+2 a c k+2 b c)+(a^2-1) e^2+2 e a c+c^2
  // This is simple 2nd order polynomial.
  // -> 0 = A x^2 + 2 B x + C
  // where A = k^2 (a^2-1)+b^2-1
  //       B = ek (a^2-1)+ack+bc
  //       C = (a^2-1)e^2+2eac+c^2
  // -> x = (-B +- sqrt[B^2 -AC])/A
  const double a = kpa*kpa*(ab*ab-1)+bb*bb-1;
  if ( std::abs(a) < 1e-9 ) return false; // Avoid divergence
  const double b = eps*kpa*(ab*ab-1) + ab*cb*kpa + bb*cb;
  const double c = (ab*ab-1)*eps*eps + 2*eps*ab*cb + cb*cb;

  // Check determinant
  double det = b*b - a*c;
  if ( det < 0 ) return false; //-1e-3 ) return false; // No real solution
  //else if ( det < 0 ) det = 0; // Ignore negative determinant due to numerical error

  nuPy1 = (-b + sqrt(det))/a;
  nuPy2 = (-b - sqrt(det))/a;
  nuPx1 = kpa*nuPy1 + eps;
  nuPx2 = kpa*nuPy2 + eps;

  return true;
}

