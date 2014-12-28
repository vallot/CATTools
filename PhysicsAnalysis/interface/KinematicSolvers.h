#ifndef CATTools_PhysicsAnalysis_KinematicSolvers_H
#define CATTools_PhysicsAnalysis_KinematicSolvers_H

#include "DataFormats/Candidate/interface/Candidate.h"

namespace cat {

class KinematicSolver
{
public:
  typedef reco::Candidate::LorentzVector LorentzVector;
  virtual void solve(const LorentzVector input[]) = 0;

  double quality() const { return quality_; };
  const LorentzVector& nu1() const { return nu1_; };
  const LorentzVector& nu2() const { return nu2_; };
  double value(int i) const { return values_.at(i); };

protected:
  double quality_;
  std::vector<double> values_;
  LorentzVector nu1_, nu2_;
};

class TTDileptonSolver : public KinematicSolver
{
public:
  TTDileptonSolver() {};
  void solve(const LorentzVector input[]) override;
};

}

#endif

