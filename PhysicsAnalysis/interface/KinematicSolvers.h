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

protected:
  double quality_;
  LorentzVector nu1_, nu2_;
  std::vector<double> values_;
};

class TTDileptonSolver : public KinematicSolver // A dummy solver for now
{
public:
  TTDileptonSolver() {};
  void solve(const LorentzVector input[]) override;
};

class MT2Solver : public KinematicSolver
{
public:
  MT2Solver() {};
  void solve(const LorentzVector input[]) override;
  double mt2();

protected:
  
};

class MAOSSolver : public MT2Solver
{
public:
  MAOSSolver(): MT2Solver() {};
  void solve(const LorentzVector input[]) override;

protected:
};

}

#endif

