#ifndef CATTools_CatAnalyzer_KinematicSolvers_H
#define CATTools_CatAnalyzer_KinematicSolvers_H

#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CLHEP/Random/RandomEngine.h"
#include "Math/LorentzVector.h"
#include <string>
#include <memory>
#include "TH1.h"

class TtFullLepKinSolver;

namespace cat {

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LV;
class KinematicSolver;

struct KinematicSolution
{
  KinematicSolution() { quality_ = -1e9; }

  double quality() const { return quality_; }
  const LV l1() const { return l1_; }
  const LV l2() const { return l2_; }
  const LV j1() const { return j1_; }
  const LV j2() const { return j2_; }
  const LV nu1() const { return nu1_; }
  const LV nu2() const { return nu2_; }
  double aux(size_t i) const { return values_.at(i); }
  const std::vector<double>& aux() const { return values_; }
  const LV w1() const { return l1()+nu1(); }
  const LV w2() const { return l2()+nu2(); }
  const LV tt() const { return t1()+t2(); }
  const LV t1() const { return t1_; }
  const LV t2() const { return t2_; }

  void setVisible(const LV& l1, const LV& l2, const LV& j1, const LV& j2)
  {
    quality_ = -1e9;
    values_.clear();
    l1_ = l1; l2_ = l2; j1_ = j1; j2_ = j2;
  }

  void setSolution(const double quality, const LV& nu1, const LV& nu2,
                   const std::vector<double> values = {})
  {
    quality_ = quality;
    values_ = values;

    nu1_ = nu1; nu2_ = nu2;

    t1_ = l1_+j1_+nu1_;
    t2_ = l2_+j2_+nu2_;
  }

  double quality_;
  LV l1_, l2_, j1_, j2_;
  LV nu1_, nu2_;
  LV t1_, t2_;
  std::vector<double> values_;

};

class KinematicSolver
{
public:
  KinematicSolver(const edm::ParameterSet& pset) {}
  virtual ~KinematicSolver() {}
  void solve(const LV met, const LV l1, const LV l2, const LV j1, const LV j2)
  {
    LV input[5] = {met, l1, l2, j1, j2};
    solve(input);
  }
  virtual void solve(const LV input[]) = 0;
  virtual std::string algoName() = 0;

  KinematicSolution solution() const { return sol_; }
  // Legacy interfaces
  double quality() const { return sol_.quality_; }
  const LV l1() const { return sol_.l1_; }
  const LV l2() const { return sol_.l2_; }
  const LV j1() const { return sol_.j1_; }
  const LV j2() const { return sol_.j2_; }
  const LV nu1() const { return sol_.nu1_; }
  const LV nu2() const { return sol_.nu2_; }
  double aux(size_t i) const { return sol_.values_.at(i); }
  const std::vector<double>& aux() const { return sol_.values_; }
  const LV w1() const { return l1()+nu1(); }
  const LV w2() const { return l2()+nu2(); }
  const LV tt() const { return t1()+t2(); }
  const LV t1() const { return sol_.t1_; }
  const LV t2() const { return sol_.t2_; }

protected:
  KinematicSolution sol_;
};

class TTDileptonSolver : public KinematicSolver // A dummy solver for now
{
public:
  TTDileptonSolver(const edm::ParameterSet& pset): KinematicSolver(pset) {}
  void solve(const LV input[]) override;
  std::string algoName() override { return "DUMMY"; }
};

class MT2Solver : public KinematicSolver
{
public:
  MT2Solver(const edm::ParameterSet& pset);
  void solve(const LV input[]) override;
  std::string algoName() override { return "MT2"; }
  double mt2();

protected:

};

class MAOSSolver : public MT2Solver
{
public:
  MAOSSolver(const edm::ParameterSet& pset): MT2Solver(pset) {}
  std::string algoName() override { return "MAOS"; }
  void solve(const LV input[]) override;
};

class CMSKinSolver : public KinematicSolver
{
public:
  CMSKinSolver(const edm::ParameterSet& pset);
  std::string algoName() override { return "CMSKIN"; }
  void solve(const LV input[]) override;

protected:
  std::unique_ptr<TtFullLepKinSolver> solver_;
};

class DESYMassLoopSolver : public KinematicSolver
{
public:
  DESYMassLoopSolver(const edm::ParameterSet& pset);
  void solve(const LV input[]) override;
  std::string algoName() override { return "DESYMassLoop"; }
protected:
  const double tMassBegin_, tMassEnd_, tMassStep_;
};

class DESYSmearedSolver : public KinematicSolver
{
public:
  DESYSmearedSolver(const edm::ParameterSet& pset);
  void solve(const LV input[]) override;
  std::string algoName() override { return "DESYSmeared"; }
  void setRandom(CLHEP::HepRandomEngine* rng) { rng_ = rng; }

protected:
  LV getSmearedLV(const LV& v, const double fE, const double dRot);
  double getRandom(TH1* h);

  CLHEP::HepRandomEngine* rng_;
  std::unique_ptr<TH1> h_jetEres_, h_jetAres_;
  std::unique_ptr<TH1> h_lepEres_, h_lepAres_;
  std::unique_ptr<TH1> h_wmass_;
  std::unique_ptr<TH1> h_mbl_w_;

  const int nTrial_;
  const double maxLBMass_, mTopInput_;
};

}

#endif

