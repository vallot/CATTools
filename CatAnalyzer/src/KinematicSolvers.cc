#include "CATTools/CatAnalyzer/interface/KinematicSolvers.h"
#include "TopQuarkAnalysis/TopKinFitter/interface/TtFullLepKinSolver.h"
#include "CATTools/CatAnalyzer/interface/TopKinSolverUtils.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_errno.h>
#include "TFile.h"
#include "TRandom3.h"

using namespace cat;
using namespace std;

void TTDileptonSolver::solve(const LV input[])
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

const double mwsqr = 80.4*80.4;
const double mbsqr = 4.8*4.8;

double mtsqr(const double qx, const double qy,
             const double px, const double py)
{
  const double ew = sqrt(mwsqr + qx*qx + qy*qy);
  const double eb = sqrt(mbsqr + px*px + py*py);
  return mwsqr + mbsqr + 2*eb*ew - 2*px*qx - 2*py*qy;
}

double fminimize(const gsl_vector* xx, void* p)
{
  const double qx1 = gsl_vector_get(xx, 0);
  const double qy1 = gsl_vector_get(xx, 1);

  double* p0 = (double*)p;
  const double qx2 = p0[4] - qx1;
  const double qy2 = p0[5] - qy1;

  const double mtA = mtsqr(qx1, qy1, p0[0], p0[1]);
  const double mtB = mtsqr(qx2, qy2, p0[2], p0[3]);

  return max(mtA, mtB);
}

}

void MT2Solver::solve(const LV input[])
{
  using namespace CATMT2;

  quality_ = -1e9; // default value
  values_.clear();

  const auto& met = input[0];
  const auto& l1 = input[1];
  const auto& l2 = input[2];
  const auto& j1 = input[3];
  const auto& j2 = input[4];

  // pars : bx1, by1, bx2, by2, Kx, Ky, xmin, xmax
  //        b : b jet px, py
  //        K : MET + lep1 + lep2
  double pars[] = {j1.px(), j1.py(), j2.px(), j2.py(),
                   met.px()+l1.px()+l2.px(), met.py()+l1.py()+l2.py()};

  // Do the minimization
  gsl_multimin_fminimizer* minimizer = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex, 2);
  gsl_vector* xx = gsl_vector_alloc(2);
  gsl_vector_set_all(xx, 0.0); // Initial value
  gsl_vector* ss = gsl_vector_alloc(2);
  gsl_vector_set_all(ss, 0.1); // Step size

  gsl_multimin_function fgsl_minimize;
  fgsl_minimize.n = 2;
  fgsl_minimize.f = fminimize;
  fgsl_minimize.params = pars;
  gsl_multimin_fminimizer_set(minimizer, &fgsl_minimize, xx, ss);

  int status = 0;
  for ( int iter=0; iter<100; ++iter )
  {
    status = gsl_multimin_fminimizer_iterate(minimizer);
    if ( status ) break;
    auto size = gsl_multimin_fminimizer_size(minimizer);
    status = gsl_multimin_test_size(size, 1e-3);

    if ( status != GSL_CONTINUE ) break;
  }
  gsl_vector_free(xx);
  gsl_vector_free(ss);

  const double qx1 = gsl_vector_get(minimizer->x, 0);
  const double qy1 = gsl_vector_get(minimizer->x, 1);

  gsl_multimin_fminimizer_free(minimizer);
  if ( status != GSL_SUCCESS ) return;

  const double qx2 = met.px()+l1.px()+l2.px()-qx1;
  const double qy2 = met.py()+l1.py()+l2.py()-qy1;
  const double mt2 = sqrt(mtsqr(qx1, qy1, pars[0], pars[1]));
  quality_ = 1;
  // For debug : print out mt2(qx1, qy1) vs mt2(qx2, qy2).
  // If minimization was failed, mt2_A != mt2_B
  // Of course this is not full test of minimization, but this must be the very first observable
  // cout << sqrt(mtsqr(qx1, qy1, pars[0], pars[1]))-sqrt(mtsqr(qx2, qy2, pars[2], pars[3])) << endl;

  const double ew1 = sqrt(mwsqr + qx1*qx1 + qy1*qy1);
  const double ew2 = sqrt(mwsqr + qx2*qx2 + qy2*qy2);
  LV w1(qx1, qy1, 0, ew1);
  LV w2(qx2, qy2, 0, ew2);
  nu1_ = w1-l1;
  nu2_ = w2-l2;

  values_.push_back(mt2);

}

void MAOSSolver::solve(const LV input[])
{
  MT2Solver::solve(input); // run the MT2 solver at the first step.

  // recover neutrino z component by assuming MT2 value as the top mass

}

CMSKinSolver::CMSKinSolver()
{
  //std::vector<double> nuPars = {30.7137,56.2880,23.0744,59.1015,24.9145}; // 7TeV pythia, from CMSSW default
  //std::vector<double> nuPars = {30.641,57.941,22.344,57.533,22.232}; // ??, from DESY
  std::vector<double> nuPars = {27.23,53.88,19.92,53.89,19.9}; // 13TeV 50ns powheg, by Youn.
  solver_.reset(new TtFullLepKinSolver(100, 300, 0.5, nuPars));
}

void CMSKinSolver::solve(const LV input[])
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

  const LV t1 = input[1]+input[3]+nu1_;
  const LV t2 = input[2]+input[4]+nu2_;
  values_.push_back((t1+t2).mass());
  values_.push_back(t1.mass());
  values_.push_back(t2.mass());
}

void DESYMassLoopSolver::solve(const LV input[])
{
  quality_ = -1e9;
  values_.clear();

  const double metX = input[0].px(), metY = input[0].py();
  const auto& l1 = input[1], l2 = input[2];
  const auto& j1 = input[3], j2 = input[4];
  const auto visSum = l1+l2+j1+j2;

  const double l1E = l1.energy(), j1E = j1.energy();
  const double l2E = l2.energy(), j2E = j2.energy();
  const double a4 = (j1E*l1.pz()-l1E*j1.pz())/l1E/(j1E+l1E);
  const double b4 = (j2E*l2.pz()-l2E*j2.pz())/l2E/(j2E+l2E);

  std::vector<double> koef, cache, sols;
  for ( double mTop = 100; mTop < 300.5; mTop += 1 ) {
    KinSolverUtils::findCoeffs(mTop, 80.4, 80.4, l1, l2, j1, j2, metX, metY, koef, cache);
    KinSolverUtils::solve_quartic(koef, a4, b4, sols);
    //const int nSol = sols.size();
    double nu1sol[4], nu2sol[4];
    for ( const double& sol : sols ) {
      // Recompute neutrino four momentum
      KinSolverUtils::getNuPxPyPzE(sol, cache, nu1sol, nu2sol);

      // Compute weight
      // This was originally Landau2D(nuE, nubarE)
      // but most recent code returns 1/mTTbar
      const double ttX = visSum.px()+nu1sol[0]+nu2sol[0];
      const double ttY = visSum.py()+nu1sol[1]+nu2sol[1];
      const double ttZ = visSum.pz()+nu1sol[2]+nu2sol[2];
      const double ttE = visSum.E()+nu1sol[3]+nu2sol[3];
      const double weight = 1/sqrt(ttE*ttE-ttX*ttX-ttY*ttY-ttZ*ttZ);

      if ( weight < quality_ ) continue;

      quality_ = weight;
      nu1_.SetPxPyPzE(nu1sol[0], nu1sol[1], nu1sol[2], nu1sol[3]);
      nu2_.SetPxPyPzE(nu2sol[0], nu2sol[1], nu2sol[2], nu2sol[3]);
    }
  }
}

DESYSmearedSolver::DESYSmearedSolver()
{
  TFile* f = TFile::Open(edm::FileInPath("CATTools/CatAnalyzer/data/desyKinRecoInput.root").fullPath().c_str());

  h_jetEres_.reset(dynamic_cast<TH1*>(f->Get("KinReco_fE_jet_step7")));
  h_jetAres_.reset(dynamic_cast<TH1*>(f->Get("KinReco_d_angle_jet_step7")));
  h_lepEres_.reset(dynamic_cast<TH1*>(f->Get("KinReco_fE_lep_step7")));
  h_lepAres_.reset(dynamic_cast<TH1*>(f->Get("KinReco_d_angle_lep_step7")));
  h_wmass_.reset(dynamic_cast<TH1*>(f->Get("KinReco_W_mass_step0")));
  h_mbl_w_.reset(dynamic_cast<TH1*>(f->Get("KinReco_mbl_true_step0")));

  h_jetEres_->SetDirectory(0);
  h_jetAres_->SetDirectory(0);
  h_lepEres_->SetDirectory(0);
  h_lepAres_->SetDirectory(0);
  h_wmass_->SetDirectory(0);
  h_mbl_w_->SetDirectory(0);

  f->Close();
}

void DESYSmearedSolver::solve(const LV input[])
{
  quality_ = -1e9;
  values_.clear();

  const double metX = input[0].px(), metY = input[0].py();
  const auto& l1 = input[1], l2 = input[2];
  const auto& j1 = input[3], j2 = input[4];

  // Skip if L+B mass is large in Smeared solution case
  const double mbl1 = (l1+j1).mass();
  const double mbl2 = (l2+j2).mass();
  if ( mbl1 > 180. or mbl2 > 180. ) return;
  const auto visSum = l1+l2+j1+j2;

  const double l1E = l1.energy(), j1E = j1.energy();
  const double l2E = l2.energy(), j2E = j2.energy();
  const double a4 = (j1E*l1.pz()-l1E*j1.pz())/l1E/(j1E+l1E);
  const double b4 = (j2E*l2.pz()-l2E*j2.pz())/l2E/(j2E+l2E);

  std::vector<double> koef, cache, sols;
  // Try 100 times with energy/angle smearing. Take weighted average of solutions
  const double mTopInput = 172.5;
  //double sumW = 0., maxSumW = 0;
  // Set random seed using kinematics (follow DESY code)
  //const unsigned int seed = std::abs(static_cast<int>(1E6*j1.pt()/j2.pt() * sin(1E6*(l1.pt() + 2.*l2.pt()))));
  //gsl_rng_mt19937
  //gRandom->SetSeed(seed);
  for ( int i=0; i<100; ++i )
  {
    // Generate smearing factors for jets and leptons
    const auto newj1 = getSmearedLV(j1, h_jetEres_->GetRandom(), h_jetAres_->GetRandom());
    const auto newj2 = getSmearedLV(j2, h_jetEres_->GetRandom(), h_jetAres_->GetRandom());
    const auto newl1 = getSmearedLV(l1, h_lepEres_->GetRandom(), h_lepAres_->GetRandom());
    const auto newl2 = getSmearedLV(l2, h_lepEres_->GetRandom(), h_lepAres_->GetRandom());

    // new MET = old MET + MET correction; MET correction = - (Vis. Pxy sum)
    const double newmetX = metX - (-visSum.px() + newj1.px() + newj2.px() + newl1.px() + newl2.px());
    const double newmetY = metY - (-visSum.py() + newj1.py() + newj2.py() + newl1.py() + newl2.py());

    KinSolverUtils::findCoeffs(mTopInput, h_wmass_->GetRandom(), h_wmass_->GetRandom(),
                               newl1, newl2, newj1, newj2, newmetX, newmetY, koef, cache);
    KinSolverUtils::solve_quartic(koef, a4, b4, sols);
    //const int nSol = sols.size();
    double nu1sol[4], nu2sol[4];
    for ( const double& sol : sols ) {
      // Recompute neutrino four momentum
      KinSolverUtils::getNuPxPyPzE(sol, cache, nu1sol, nu2sol);

      // Compute weight by m(B,L)
      const double w1 = h_mbl_w_->GetBinContent(h_mbl_w_->FindBin(mbl1));
      const double w2 = h_mbl_w_->GetBinContent(h_mbl_w_->FindBin(mbl2));
      const double weight = w1*w2;

      if ( weight < quality_ ) continue;

      quality_ = weight;
      nu1_.SetPxPyPzE(nu1sol[0], nu1sol[1], nu1sol[2], nu1sol[3]);
      nu2_.SetPxPyPzE(nu2sol[0], nu2sol[1], nu2sol[2], nu2sol[3]);
    }
  }
}

LV DESYSmearedSolver::getSmearedLV(const LV& lv0,
                                   const double fE, const double dRot)
{
  // Rescale at the first step
  const double e = fE*lv0.energy();
  const double p = sqrt(std::max(0., e*e-lv0.M2()));
  if ( KinSolverUtils::isZero(e) or KinSolverUtils::isZero(p) ) return LV();

  // Apply rotation
  const double localPhi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());
  const double theta = lv0.Theta() + dRot*cos(localPhi);
  const double phi = lv0.Phi() + dRot*sin(localPhi);

  const double pz = p*cos(theta);
  const double px = p*sin(theta)*cos(phi);
  const double py = p*sin(theta)*sin(phi);

  return LV(px, py, pz, e);
}


namespace CATNuWeight
{

bool computeNuPxPy(const LV& lep,
                   const LV& jet,
                   const double mT, const double nuEta,
                   double& nuPx1, double& nuPy1, double& nuPx2, double& nuPy2)
{
  const double mTsqr = mT*mT;
  const double mWsqr = 80.4*80.4;
  const double mBsqr = 4.8*4.8;
  const double mLsqr = 0;

  const double alpha1 = (mTsqr - mBsqr - mWsqr)/2; // eqn B.6
  const double alpha2 = (mWsqr - mLsqr)/2; // eqn B.7

  const double pxl = lep.px(), pyl = lep.py(), pzl = lep.pz(), el = lep.energy();
  const double pxb = jet.px(), pyb = jet.py(), pzb = jet.pz(), eb = jet.energy();

  const double denB = eb*cosh(nuEta) - pzb*sinh(nuEta); // eqn B.14
  const double denL = el*cosh(nuEta) - pzl*sinh(nuEta); // eqn B.14

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

double nuWeight(const LV& l1,
                const LV& l2,
                const LV& j1,
                const LV& j2,
                const double metX, const double metY,
                const double sigmaEXsqr, const double sigmaEYsqr,
                const double mt, const double nu1Eta, const double nu2Eta,
                double& bestPx1, double& bestPy1, double& bestPx2, double& bestPy2)
{
  double tmpSolX1[2], tmpSolY1[2], tmpSolX2[2], tmpSolY2[2];
  if ( !computeNuPxPy(l1, j1, mt, nu1Eta, tmpSolX1[0], tmpSolY1[0], tmpSolX1[1], tmpSolY1[1]) ) return -1;
  if ( !computeNuPxPy(l2, j2, mt, nu2Eta, tmpSolX2[0], tmpSolY2[0], tmpSolX2[1], tmpSolY2[1]) ) return -1;

  double bestWeight = -1e9;
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
        bestPy1 = tmpSolX1[i];
        bestPx2 = tmpSolX2[j];
        bestPy2 = tmpSolX2[j];
      }
    }
  }

  return bestWeight;
}

double fminimize(const gsl_vector* xx, void* p)
{
  const double mt = gsl_vector_get(xx, 0);
  const double nu1Eta = gsl_vector_get(xx, 1);
  const double nu2Eta = gsl_vector_get(xx, 2);

  double* p0 = (double*) p;
  LV l1(p0[ 0], p0[ 1], p0[ 2], p0[ 3]);
  LV l2(p0[ 4], p0[ 5], p0[ 6], p0[ 7]);
  LV j1(p0[ 8], p0[ 9], p0[10], p0[11]);
  LV j2(p0[12], p0[13], p0[14], p0[15]);
  const double metX = p0[16], metY = p0[17];
  const double sigmaEXsqr = p0[18], sigmaEYsqr = p0[19];

  double bestPx1, bestPy1, bestPx2, bestPy2;
  const double weight = nuWeight(l1, l2, j1, j2, metX, metY, sigmaEXsqr, sigmaEYsqr,
                                 mt, nu1Eta, nu2Eta,
                                 bestPx1, bestPy1, bestPx2, bestPy2);
  if ( weight <= 0 ) return 1e9;
  return 1/(0.01+weight);
}

}

void NuWeightSolver::solve(const LV input[])
{
  using namespace CATNuWeight;

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
  double bestMt = 0;
  double bestPx1 = 0, bestPx2 = 0, bestPy1 = 0, bestPy2 = 0;
  double bestEta1 = 0, bestEta2 = 0;
  for ( double mt = 100; mt < 300; mt += 2 )
  {
    for ( double nu1Eta = -5.0; nu1Eta <= 5.0; nu1Eta += 0.5 )
    {
      for ( double nu2Eta = -5.0; nu2Eta <= 5.0; nu2Eta += 0.5 )
      {
        double nu1Px, nu1Py, nu2Px, nu2Py;
        const double weight = nuWeight(l1, l2, j1, j2, metX, metY,  sigmaEXsqr, sigmaEYsqr,
                                       mt, nu1Eta, nu2Eta,
                                       nu1Px, nu1Py, nu2Px, nu2Py);
        if ( weight <= 0 ) continue;
        if ( weight > bestWeight )
        {
          bestWeight = weight;
          bestMt = mt;
          bestPx1 = nu1Px; bestPy1 = nu1Py;
          bestPx2 = nu2Px; bestPy2 = nu2Py;
          bestEta1 = nu1Eta; bestEta2 = nu2Eta;
        }
      }
    }
  }

  double pars[] = {
    l1.px(), l1.py(), l1.pz(), l1.e(),
    l2.px(), l2.py(), l2.pz(), l2.e(),
    j1.px(), j1.py(), j1.pz(), j1.e(),
    j2.px(), j2.py(), j2.pz(), j2.e(),
    metX, metY,
    sigmaEXsqr, sigmaEYsqr
  };

  gsl_multimin_fminimizer* minimizer = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex, 3);
  gsl_vector* xx = gsl_vector_alloc(3);
  gsl_vector_set(xx, 0, bestMt);
  gsl_vector_set(xx, 1, bestEta1);
  gsl_vector_set(xx, 2, bestEta2);
  gsl_vector* ss = gsl_vector_alloc(3);
  gsl_vector_set_all(ss, 0.1);

  gsl_multimin_function fgsl_minimize;
  fgsl_minimize.n = 3;
  fgsl_minimize.f = fminimize;
  fgsl_minimize.params = pars;
  gsl_multimin_fminimizer_set(minimizer, &fgsl_minimize, xx, ss);

  int status = 0;
  for ( int iter=0; iter<100; ++iter )
  {
    status = gsl_multimin_fminimizer_iterate(minimizer);
    if ( status ) break;
    auto size = gsl_multimin_fminimizer_size(minimizer);
    status = gsl_multimin_test_size(size, 1e-3);
    if ( status != GSL_CONTINUE ) break;
  }
  gsl_vector_free(xx);
  gsl_vector_free(ss);
  if ( status == GSL_SUCCESS )
  {
    bestMt = gsl_vector_get(minimizer->x, 0);
    bestEta1 = gsl_vector_get(minimizer->x, 1);
    bestEta2 = gsl_vector_get(minimizer->x, 2);
    quality_ = nuWeight(l1, l2, j1, j2, metX, metY, sigmaEXsqr, sigmaEYsqr,
                        bestMt, bestEta1, bestEta2,
                        bestPx1, bestPy1, bestPx2, bestPy2);
  }
  gsl_multimin_fminimizer_free(minimizer);

  const double bestPz1 = hypot(bestPx1, bestPy1)*sinh(bestEta1);
  const double bestPz2 = hypot(bestPx2, bestPy2)*sinh(bestEta2);
  nu1_ = LV(bestPx1, bestPy1, bestPz1, sqrt(bestPx1*bestPx1+bestPy1*bestPy1+bestPz1*bestPz1));
  nu2_ = LV(bestPx2, bestPy2, bestPz2, sqrt(bestPx2*bestPx2+bestPy2*bestPy2+bestPz2*bestPz2));

  const LV t1 = l1+j1+nu1_;
  const LV t2 = l2+j2+nu2_;
  values_.push_back((t1+t2).mass());
  values_.push_back(t1.mass());
  values_.push_back(t2.mass());
}


