#include "CATTools/PhysicsAnalysis/interface/KinematicSolvers.h"

using namespace cat;

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

void CMSKinSolver::solve(const KinematicSolver::LorentzVector input[])
{
}

