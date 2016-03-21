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
  const auto& met = input[0];
  const auto& l1 = input[1];
  const auto& l2 = input[2];
  const auto& j1 = input[3];
  const auto& j2 = input[4];
  sol_.setVisible(l1, l2, j1, j2);

  auto nu1 = met/2;
  auto nu2 = met/2;

  const double mW = 80.4;
  if ( abs((nu1+l1).mass()-mW)+abs((nu2+l2).mass()-mW) <
       abs((nu1+l2).mass()-mW)+abs((nu2+l1).mass()-mW) ) swap(nu1, nu2);

  double quality = abs((nu1+l1+j1).mass()-172.5) + abs((nu2+l2+j2).mass()-172.5);
  quality = 0.01/(0.01+quality);

  //const auto& lb1 = l1+j1;
  //const auto& lb2 = l2+j2;
  //if ( lb1.mass() > 150 or lb2.mass() > 150 ) sol_.quality_ = 0;

  const auto& vsum = l1+l2+j1+j2+met;
  if ( vsum.mass() < 300 ) quality = 0;

  vector<double> values = {vsum.mass()};
  sol_.setSolution(quality, nu1, nu2, values);
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

MT2Solver::MT2Solver(const edm::ParameterSet& pset): KinematicSolver(pset)
{
}

void MT2Solver::solve(const LV input[])
{
  using namespace CATMT2;

  const auto& met = input[0];
  const auto& l1 = input[1];
  const auto& l2 = input[2];
  const auto& j1 = input[3];
  const auto& j2 = input[4];

  sol_.setVisible(l1, l2, j1, j2);

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
  double quality = 1;
  // For debug : print out mt2(qx1, qy1) vs mt2(qx2, qy2).
  // If minimization was failed, msol_.t2_A != msol_.t2_B
  // Of course this is not full test of minimization, but this must be the very first observable
  // cout << sqrt(mtsqr(qx1, qy1, pars[0], pars[1]))-sqrt(mtsqr(qx2, qy2, pars[2], pars[3])) << endl;

  const double ew1 = sqrt(mwsqr + qx1*qx1 + qy1*qy1);
  const double ew2 = sqrt(mwsqr + qx2*qx2 + qy2*qy2);
  LV w1(qx1, qy1, 0, ew1);
  LV w2(qx2, qy2, 0, ew2);

  vector<double> values = {mt2};

  sol_.setSolution(quality, w1-l1, w2-l2, values);

}

void MAOSSolver::solve(const LV input[])
{
  MT2Solver::solve(input); // run the MT2 solver at the first step.

  // recover neutrino z component by assuming MT2 value as the top mass

}

CMSKinSolver::CMSKinSolver(const edm::ParameterSet& pset): KinematicSolver(pset)
{
  const double tMassBegin = pset.getParameter<double>("tMassBegin");
  const double tMassEnd   = pset.getParameter<double>("tMassEnd");
  const double tMassStep  = pset.getParameter<double>("tMassStep");
  std::vector<double> nuPars = pset.getParameter<std::vector<double> >("nuPars");
  solver_.reset(new TtFullLepKinSolver(tMassBegin, tMassEnd, tMassStep, nuPars));
}

void CMSKinSolver::solve(const LV input[])
{
  sol_.setVisible(input[1], input[2], input[3], input[4]);

  const double metX = input[0].px();
  const double metY = input[0].py();
  const TLorentzVector l1(input[1].px(), input[1].py(), input[1].pz(), input[1].e());
  const TLorentzVector l2(input[2].px(), input[2].py(), input[2].pz(), input[2].e());
  const TLorentzVector j1(input[3].px(), input[3].py(), input[3].pz(), input[3].e());
  const TLorentzVector j2(input[4].px(), input[4].py(), input[4].pz(), input[4].e());

  solver_->SetConstraints(metX+l1.Px()+j1.Px()+l2.Px()+j2.Px(),
                          metY+l1.Py()+j1.Py()+l2.Py()+j2.Py());
  auto nuSol = solver_->getNuSolution(l1, l2, j1, j2);
  const double quality  = nuSol.weight;
  if ( quality < 0 ) return;

  auto nu1 = nuSol.neutrino.p4();
  auto nu2 = nuSol.neutrinoBar.p4();

  const LV t1 = input[1]+input[3]+nu1;
  const LV t2 = input[2]+input[4]+nu2;
  vector<double> values = {(t1+t2).mass(), t1.mass(), t2.mass()};

  sol_.setSolution(quality, nu1, nu2, values);
}

DESYMassLoopSolver::DESYMassLoopSolver(const edm::ParameterSet& pset):
  KinematicSolver(pset),
  tMassBegin_(pset.getParameter<double>("tMassBegin")),
  tMassEnd_(pset.getParameter<double>("tMassEnd")),
  tMassStep_(pset.getParameter<double>("tMassStep"))
{
}

void DESYMassLoopSolver::solve(const LV input[])
{
  const double metX = input[0].px(), metY = input[0].py();
  const auto& l1 = input[1], l2 = input[2];
  const auto& j1 = input[3], j2 = input[4];
  sol_.setVisible(l1, l2, j1, j2);

  const auto visSum = l1+l2+j1+j2;

  const double l1E = l1.energy(), j1E = j1.energy();
  const double l2E = l2.energy(), j2E = j2.energy();
  const double a4 = (j2E*l2.pz()-l2E*j2.pz())/l2E/(j2E+l2E);
  const double b4 = (j1E*l1.pz()-l1E*j1.pz())/l1E/(j1E+l1E);

  double quality = sol_.quality(); // Just take the default quality value for initial input
  math::XYZTLorentzVector nu1, nu2;
  std::vector<double> koef, cache, sols;
  for ( double mTop = tMassBegin_; mTop < tMassEnd_+0.5*tMassStep_; mTop += tMassStep_ ) {
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

      if ( weight < quality ) continue;

      quality = weight;
      nu1.SetPxPyPzE(nu1sol[0], nu1sol[1], nu1sol[2], nu1sol[3]);
      nu2.SetPxPyPzE(nu2sol[0], nu2sol[1], nu2sol[2], nu2sol[3]);
    }
  }
  sol_.setSolution(quality, nu1, nu2);
}

DESYSmearedSolver::DESYSmearedSolver(const edm::ParameterSet& pset):
  KinematicSolver(pset),
  nTrial_(pset.getParameter<int>("nTrial")),
  maxLBMass_(pset.getParameter<double>("maxLBMass")),
  mTopInput_(pset.getParameter<double>("mTopInput"))
{
  rng_ = 0;

  const auto filePath = pset.getParameter<string>("inputTemplatePath");
  TFile* f = TFile::Open(edm::FileInPath(filePath).fullPath().c_str());

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
  const double metX = input[0].px(), metY = input[0].py();
  auto l1 = input[1], l2 = input[2];
  auto j1 = input[3], j2 = input[4];
  sol_.setVisible(l1, l2, j1, j2);

  // Skip if L+B mass is large in Smeared solution case
  const double mbl1 = (l1+j1).mass();
  const double mbl2 = (l2+j2).mass();
  if ( mbl1 > maxLBMass_ or mbl2 > maxLBMass_ ) return;
  const auto visSum = l1+l2+j1+j2;

  // Continue to get smeared solution if exact solver fails
  // Try 100 times with energy/angle smearing. Take weighted average of solutions
  //double sumW = 0., maxSumW = 0;
  // Set random seed using kinematics (follow DESY code)
  //const unsigned int seed = std::abs(static_cast<int>(1E6*j1.pt()/j2.pt() * sin(1E6*(l1.pt() + 2.*l2.pt()))));
  //gsl_rng_mt19937
  //gRandom->SetSeed(seed);
  double sumW = 0;
  double sumP[6][3] = {{0,},};
  for ( int i=0; i<nTrial_; ++i )
  {
    std::vector<double> koef, cache, sols;

    // Generate smearing factors for jets and leptons
    const auto newl1 = getSmearedLV(l1, getRandom(h_lepEres_.get()), getRandom(h_lepAres_.get()));
    const auto newl2 = getSmearedLV(l2, getRandom(h_lepEres_.get()), getRandom(h_lepAres_.get()));
    const auto newj1 = getSmearedLV(j1, getRandom(h_jetEres_.get()), getRandom(h_jetAres_.get()));
    const auto newj2 = getSmearedLV(j2, getRandom(h_jetEres_.get()), getRandom(h_jetAres_.get()));
    const double newl1E = newl1.E(), newl2E = newl2.E();
    const double newj1E = newj1.E(), newj2E = newj2.E();
    const double a4 = (newj2E*newl2.pz()-newl2E*newj2.pz())/newl2E/(newj2E+newl2E);
    const double b4 = (newj1E*newl1.pz()-newl1E*newj1.pz())/newl1E/(newj1E+newl1E);

    // new MET = old MET + MET correction; MET correction = - (Vis. Pxy sum)
    const auto newVisSum = newl1 + newl2 + newj1 + newj2;
    const double newmetX = metX + visSum.px() - newVisSum.px();
    const double newmetY = metY + visSum.py() - newVisSum.py();

    // Compute weight by m(B,L)
    const double w1 = h_mbl_w_->GetBinContent(h_mbl_w_->FindBin((newl1+newj1).mass()));
    const double w2 = h_mbl_w_->GetBinContent(h_mbl_w_->FindBin((newl2+newj2).mass()));
    double weight = w1*w2/h_mbl_w_->Integral()/h_mbl_w_->Integral();
    if ( weight <= 0 ) continue;

//#define KINRECONOSM
#ifdef KINRECONOSM
    KinSolverUtils::findCoeffs(mTopInput_, 80.4, 80.4,
                               newl1, newl2, newj1, newj2, newmetX, newmetY, koef, cache);
#else
    KinSolverUtils::findCoeffs(mTopInput_, getRandom(h_wmass_.get()), getRandom(h_wmass_.get()),
                               newl1, newl2, newj1, newj2, newmetX, newmetY, koef, cache);
#endif
    KinSolverUtils::solve_quartic(koef, a4, b4, sols);
    double nu1sol[4] = {}, nu2sol[4] = {};
    // Choose one solution with minimal mass of top pair
    double maxWeightSol = 0;
    for ( const double& sol : sols ) {
      // Recompute neutrino four momentum
      double nu1solTmp[4], nu2solTmp[4];
      KinSolverUtils::getNuPxPyPzE(sol, cache, nu1solTmp, nu2solTmp);
      bool hasNan = false;
      for ( int i=0; i<4; ++i ) {
        if ( std::isnan(nu1solTmp[i]) or std::isnan(nu2solTmp[i]) ) {
          hasNan = true;
          break;
        }
      }
      if ( hasNan ) continue;
      const LV nu1(nu1solTmp[0], nu1solTmp[1], nu1solTmp[2], nu1solTmp[3]);
      const LV nu2(nu2solTmp[0], nu2solTmp[1], nu2solTmp[2], nu2solTmp[3]);

/*
      const double nw1 = h_mbl_w_->GetBinContent(h_mbl_w_->FindBin((nu1+newj1).mass()));
      const double nw2 = h_mbl_w_->GetBinContent(h_mbl_w_->FindBin((nu2+newj2).mass()));
      const double nw = nw1*nw2/h_mbl_w_->Integral()/h_mbl_w_->Integral();
      if ( nw <= 0 ) continue;
*/

      const double nw = 1./(nu1+nu2+newVisSum).mass();
      if ( nw > maxWeightSol ) {
        maxWeightSol = nw;
        std::copy(nu1solTmp, nu1solTmp+4, nu1sol);
        std::copy(nu2solTmp, nu2solTmp+4, nu2sol);
      }
    }
    if ( maxWeightSol <= 0 ) continue;

    weight = sqrt(weight);
    sumW += weight;
    sumP[0][0] += weight*newl1.px(); sumP[0][1] += weight*newl1.py(); sumP[0][2] += weight*newl1.pz();
    sumP[1][0] += weight*newl2.px(); sumP[1][1] += weight*newl2.py(); sumP[1][2] += weight*newl2.pz();
    sumP[2][0] += weight*newj1.px(); sumP[2][1] += weight*newj1.py(); sumP[2][2] += weight*newj1.pz();
    sumP[3][0] += weight*newj2.px(); sumP[3][1] += weight*newj2.py(); sumP[3][2] += weight*newj2.pz();
    sumP[4][0] += weight*nu1sol[0]; sumP[4][1] += weight*nu1sol[1]; sumP[4][2] += weight*nu1sol[2];
    sumP[5][0] += weight*nu2sol[0]; sumP[5][1] += weight*nu2sol[1]; sumP[5][2] += weight*nu2sol[2];
  }
  if ( sumW <= 0 ) return;
  for ( int i=0; i<6; ++i ) for ( int j=0; j<3; ++j ) sumP[i][j] /= sumW;

  l1.SetXYZT(sumP[0][0], sumP[0][1], sumP[0][2], KinSolverUtils::computeEnergy(sumP[0], l1.M()));//KinSolverUtils::mL));
  l2.SetXYZT(sumP[1][0], sumP[1][1], sumP[1][2], KinSolverUtils::computeEnergy(sumP[1], l2.M()));//KinSolverUtils::mL));
  j1.SetXYZT(sumP[2][0], sumP[2][1], sumP[2][2], KinSolverUtils::computeEnergy(sumP[2], j1.M()));//KinSolverUtils::mB));
  j2.SetXYZT(sumP[3][0], sumP[3][1], sumP[3][2], KinSolverUtils::computeEnergy(sumP[3], j2.M()));//KinSolverUtils::mB));

  sol_.setVisible(l1, l2, j1, j2);
  const math::XYZTLorentzVector nu1(sumP[4][0], sumP[4][1], sumP[4][2], KinSolverUtils::computeEnergy(sumP[4], KinSolverUtils::mV));
  const math::XYZTLorentzVector nu2(sumP[5][0], sumP[5][1], sumP[5][2], KinSolverUtils::computeEnergy(sumP[5], KinSolverUtils::mV));
  sol_.setSolution(sumW, nu1, nu2);
  sol_.t1_ = l1+j1+nu1;
  sol_.t2_ = l2+j2+nu2;
  sol_.t1_.SetE(sqrt(sol_.t1_.P2() + mTopInput_*mTopInput_));
  sol_.t2_.SetE(sqrt(sol_.t2_.P2() + mTopInput_*mTopInput_));

}

LV DESYSmearedSolver::getSmearedLV(const LV& lv0,
                                   const double fE, const double dRot)
{
#ifdef KINRECONOSM
  return lv0;
#endif
  // Rescale at the first step
  const double e = fE*lv0.energy();
  const double p = std::sqrt(std::max(0., e*e-lv0.M2()));
  if ( KinSolverUtils::isZero(e) or KinSolverUtils::isZero(p) ) return LV();

  const double px0 = std::abs(lv0.px()) < 0.001 ? 0 : lv0.px();
  const double py0 = std::abs(lv0.py()) < 0.001 ? 0 : lv0.py();
  const double pz0 = std::abs(lv0.pz()) < 0.001 ? 0 : lv0.pz();
  if ( p == 0 ) return LV(0, 0, 0, e);

  // Apply rotation
  const double localPhi = 2*TMath::Pi()*rng_->flat();
  const double px1 = -p*sin(dRot)*sin(localPhi);
  const double py1 =  p*sin(dRot)*cos(localPhi);
  const double pz1 =  p*cos(dRot);

  if ( py0 == 0 and pz0 == 0 ) return LV(pz1, px1, py1, e);

  const double d = hypot(pz0, py0);

  const double x1 =  d/p        , y1 =  0    , z1 = px0/p;
  const double x2 = -px0*py0/d/p, y2 =  pz0/d, z2 = py0/p;
  const double x3 = -px0*pz0/d/p, y3 = -py0/d, z3 = pz0/p;

  const double px = x1*px1 + y1*py1 + z1*pz1;
  const double py = x2*px1 + y2*py1 + z2*pz1;
  const double pz = x3*px1 + y3*py1 + z3*pz1;

  return LV(px, py, pz, e);
}

double DESYSmearedSolver::getRandom(TH1* h)
{
  if ( !h ) return 0;

  //h->GetRandom();
  const int n = h->GetNbinsX();
  const double* fIntegral = h->GetIntegral();
  const double integral = fIntegral[n];

  if (integral == 0) return 0;

  const double r1 = rng_->flat();
  const int ibin = TMath::BinarySearch(n, fIntegral, r1);
  double x = h->GetBinLowEdge(ibin+1);
  const double binW = h->GetBinWidth(ibin+1);
  const double y1 = fIntegral[ibin], y2 = fIntegral[ibin+1];
  if ( r1 > y1 ) x += binW*(r1-y1)/(y2-y1);

  return x;
}

