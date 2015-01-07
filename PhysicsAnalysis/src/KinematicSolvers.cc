#include "CATTools/PhysicsAnalysis/interface/KinematicSolvers.h"
#include "TopQuarkAnalysis/TopKinFitter/interface/TtFullLepKinSolver.h"
#include <boost/assign/std/vector.hpp> // for 'operator+=()'

using namespace cat;
using namespace boost::assign;
using namespace std;

void TTDileptonSolver::solve(const KinematicSolver::LorentzVector input[])
{
  quality_ = -1e9; // default value
  values_.clear();

  const auto& l1 = input[0];
  const auto& l2 = input[1];
  const auto& j1 = input[2];
  const auto& j2 = input[3];
  const auto& met = input[4];

  nu1_ = met/2;
  nu2_ = met/2;

  const auto& lb1 = l1+j1;
  const auto& lb2 = l2+j2;
  if ( lb1.mass() > 150 or lb2.mass() > 150 ) quality_ = 0;

  const auto& vsum = l1+l2+j1+j2+met;
  if ( vsum.mass() < 300 ) quality_ = 0;

  values_.push_back(vsum.mass());
}

void MT2Solver::solve(const KinematicSolver::LorentzVector input[])
{
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
  solver_ = new TtFullLepKinSolver(150, 200, 1, nuPars, 80.4, 4.8);
}

CMSKinSolver::~CMSKinSolver()
{
  if ( solver_ ) delete solver_;
}

void CMSKinSolver::solve(const KinematicSolver::LorentzVector input[])
{
  quality_ = -1e9; // default value
  values_.clear();

  const TLorentzVector l1(input[0].px(), input[0].py(), input[0].pz(), input[0].e());
  const TLorentzVector l2(input[1].px(), input[1].py(), input[1].pz(), input[1].e());
  const double metX = input[2].px();
  const double metY = input[2].py();
  const TLorentzVector j1(input[3].px(), input[3].py(), input[3].pz(), input[3].e());
  const TLorentzVector j2(input[4].px(), input[4].py(), input[4].pz(), input[4].e());

  solver_->SetConstraints(metX, metY);
  auto nuSol = solver_->getNuSolution(l1, l2, j1, j2);
  quality_  = nuSol.weight;
  nu1_ = nuSol.neutrino.p4();
  nu2_ = nuSol.neutrinoBar.p4();
}

void NuWeightSolver::solve(const KinematicSolver::LorentzVector input[])
{
//cout << "STARTING NUSOLVER" << endl;
  quality_ = -1e9;
  values_.clear();

  const auto& l1 = input[0];
  const auto& l2 = input[1];
  const double metX = input[2].px();
  const double metY = input[2].py();
  const auto& j1 = input[3];
  const auto& j2 = input[4];

  const double sigmaEXsqr = 1, sigmaEYsqr = 1;

  double bestWeight = -1e9;
  double bestPx1 = 1e9, bestPx2 = 1e9, bestPy1 = 1e9, bestPy2 = -1e9;
  double bestEta1 = -1e9, bestEta2 = -1e9;
  double tmpSolX1[2], tmpSolY1[2], tmpSolX2[2], tmpSolY2[2];
  for ( double mt = 140; mt < 220; mt += 1 )
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
//cout << "New MTop " << mt << ' ' << weight << ' ' << hypot(bestPx1, bestPy1) << ' ' << bestEta1 << ':' << hypot(bestPx2, bestPy2) << ' ' << bestEta2 << endl;
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
}

bool NuWeightSolver::computeNuPxPy(const KinematicSolver::LorentzVector& lep,
                                   const KinematicSolver::LorentzVector& jet,
                                   const double mT, const double nuEta,
                                   double& nuPx1, double& nuPy1, double& nuPx2, double& nuPy2) const
{
  const double mTsqr = mT*mT;
  const static double mWsqr = 80.4*80.4;
  const static double mBsqr = 4.8*4.8;
  const static double mLsqr = 0;

  const double alpha1 = (mTsqr - mBsqr - mWsqr)/2; // eqn B.6
  const static double alpha2 = (mWsqr - mLsqr)/2; // eqn B.7

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
  // -> 0 = C x^2 + 2 B x^2 + C
  // -> x = (-b +- sqrt[b^2 - ac])/a
  // This is simple 2nd order polynomial.
  const double a = kpa*kpa*(ab*ab-1)+bb*bb-1;
  if ( std::abs(a) < 1e-9 ) return false; // Avoid divergence
  const double b = eps*kpa*(ab*ab-1) + ab*cb*kpa + bb*cb;
  const double c = (ab*ab-1)*eps*eps + 2*eps*ab*cb + cb*cb;

  // Check determinant
  double det = b*b - a*c;
  if ( det < -1e-3 ) return false; // No real solution
  else if ( det < 0 ) det = 0; // Ignore negative determinant due to numerical error

  nuPx1 = -b/a + sqrt(det);
  nuPx2 = -b/a - sqrt(det);
  nuPy1 = kpa*nuPx1 + eps;
  nuPy2 = kpa*nuPx2 + eps;

  return true;
}

