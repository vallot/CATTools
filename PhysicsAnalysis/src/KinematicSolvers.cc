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

  const auto& vsum = l1+l2+j1+j2+met;

  values_.push_back(vsum.mass());
}
