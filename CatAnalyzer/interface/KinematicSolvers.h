#ifndef CATTools_CatAnalyzer_KinematicSolvers_H
#define CATTools_CatAnalyzer_KinematicSolvers_H

#include "DataFormats/Candidate/interface/Candidate.h"
#include <memory>
#include "TH1.h"

class TtFullLepKinSolver;

namespace cat {

class KinematicSolver
{
public:
  virtual ~KinematicSolver() {};
  typedef reco::Candidate::LorentzVector LorentzVector;
  virtual void solve(const LorentzVector input[]) = 0;

  double quality() const { return quality_; };
  const LorentzVector& nu1() const { return nu1_; };
  const LorentzVector& nu2() const { return nu2_; };
  double aux(size_t i) const { return values_.at(i); };
  const std::vector<double>& aux() const { return values_; };

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
};

class CMSKinSolver : public KinematicSolver
{
public:
  CMSKinSolver();
  void solve(const LorentzVector input[]) override;

protected:
  std::unique_ptr<TtFullLepKinSolver> solver_;
};

class DESYMassLoopSolver : public KinematicSolver
{
public:
  DESYMassLoopSolver();
  void solve(const LorentzVector input[]) override;
};

class DESYSmearedSolver : public KinematicSolver
{
public:
  DESYSmearedSolver();
  void solve(const LorentzVector input[]) override;

protected:
  math::XYZTLorentzVector getSmearedLV(const math::XYZTLorentzVector& v, const double fE, const double dRot);

  std::unique_ptr<TH1> h_jetEres_, h_jetAres_;
  std::unique_ptr<TH1> h_lepEres_, h_lepAres_;
  std::unique_ptr<TH1> h_wmass_;
  std::unique_ptr<TH1> h_mbl_w_;
};

// Neutrino weighting method (from thesis by Temple)
class NuWeightSolver : public KinematicSolver
{
public:
  void solve(const LorentzVector input[]) override;
};

}

#endif

