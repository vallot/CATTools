#ifndef TTLJKinFit_cxx
#define TTLJKinFit_cxx

#include <Minuit2/Minuit2Minimizer.h>
#include <TLorentzVector.h>
#include <memory>
#include <vector>

struct TTLJKinFitFtn : public ROOT::Math::IBaseFunctionMultiDim
{
  TTLJKinFitFtn(const double hadW_m, const double hadT_m,
                const double lepW_m, const double lepT_m);
  IBaseFunctionMultiDim* Clone() const override;

  unsigned int NDim() const override;
  double DoEval(const double* x) const override;

  TLorentzVector solveLepTopNu(TLorentzVector nu, const TLorentzVector& lep, const TLorentzVector& lb) const;
  std::vector<TLorentzVector> getSolution(const double* x) const;

  const double hadW_m_, hadT_m_, lepW_m_, lepT_m_;
  TLorentzVector lvs_[6];
};

struct TTLJKinFit
{
  TTLJKinFit(const double hadW_m = 80.4, const double hadT_m = 172.5,
         const double lepW_m = 80.4, const double lepT_m = 172.5);

  double compute(const TLorentzVector metP4, const TLorentzVector l1, const TLorentzVector b1,
                 const TLorentzVector wj1, const TLorentzVector wj2, const TLorentzVector hb);

  std::vector<TLorentzVector> getSolution() const;

  std::unique_ptr<ROOT::Math::Minimizer> min0_, min_;
  TTLJKinFitFtn ftn_;
  double chi2_;

};

#endif
