#include "CATTools/PhysicsAnalysis/interface/KinematicSolvers.h"
#include "TopQuarkAnalysis/TopKinFitter/interface/TtFullLepKinSolver.h"
#include <boost/assign/std/vector.hpp> // for 'operator+=()'

using namespace cat;
using namespace boost::assign;

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
  const TLorentzVector j1(input[2].px(), input[2].py(), input[2].pz(), input[2].e());
  const TLorentzVector j2(input[3].px(), input[3].py(), input[3].pz(), input[3].e());
  const double metX = input[4].px();
  const double metY = input[4].py();

  solver_->SetConstraints(metX, metY);
  auto nuSol = solver_->getNuSolution(l1, l2, j1, j2);
  quality_  = nuSol.weight;
  nu1_ = nuSol.neutrino.p4();
  nu2_ = nuSol.neutrinoBar.p4();
}
