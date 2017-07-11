#include "CATTools/CatAnalyzer/interface/TTLJKinFit.h"

#include <TLorentzVector.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <vector>
#include <cmath>

using namespace std;

TTLJKinFitFtn::TTLJKinFitFtn(const double hadW_m, const double hadT_m,
                             const double lepW_m, const double lepT_m):
  ROOT::Math::IBaseFunctionMultiDim(),
  hadW_m_(hadW_m), hadT_m_(hadT_m), lepW_m_(lepW_m), lepT_m_(lepT_m)
{
}

TTLJKinFitFtn::IBaseFunctionMultiDim* TTLJKinFitFtn::Clone() const
{
  TTLJKinFitFtn* obj = new TTLJKinFitFtn(hadW_m_, hadT_m_, lepW_m_, lepT_m_);
  for ( size_t i=0; i<6; ++i ) obj->lvs_[i] = this->lvs_[i];
  return obj;
}

unsigned int TTLJKinFitFtn::NDim() const { return 1; }

double TTLJKinFitFtn::DoEval(const double* x) const
{
  const auto sol = getSolution(x);
  const auto nu = sol[0], lep = sol[1], lb = sol[2];
  const auto wj1 = sol[3], wj2 = sol[4], hb = sol[5];

  const double wLepDm = (lep+nu).M()-lepW_m_;
  const double tLepDm = (lep+nu+lb).M()-lepT_m_;
  const double wHadDm = (wj1+wj2).M()-hadW_m_;
  const double tHadDm = (wj1+wj2+hb).M()-hadT_m_;

  return wLepDm*wLepDm + tLepDm*tLepDm + wHadDm*wHadDm + tHadDm*tHadDm;
}

TLorentzVector TTLJKinFitFtn::solveLepTopNu(TLorentzVector nu, const TLorentzVector& lep, const TLorentzVector& lb) const
{
  // Solve neutrino pz with M=80.4 constraint, M^2 = (E_l+E_v)^2 - (P_l+P_v)^2
  // There are two fold ambiguity, resolve one of them by minimizing the chi2
  const double pp = lepW_m_*lepW_m_/2 + (nu.Px()*lep.Px()+nu.Py()*lep.Py());
  const double a = pp*lep.Pz()/lep.Perp2();
  const double det = a*a + (pp*pp - lep.P()*lep.P()*nu.Perp2())/lep.Perp2();
  const double b = sqrt(max(0., det));
  //const double wpx = nu.Px()+lep.Px();
  //const double wpy = nu.Py()+lep.Py();
  TLorentzVector nuPos, nuNeg;
  nuPos.SetXYZM(nu.Px(), nu.Py(), a+b, 0);
  nuNeg.SetXYZM(nu.Px(), nu.Py(), a-b, 0);
  const auto wPos = lep+nuPos, wNeg = lep+nuNeg;

  const double tLepDmPos = abs((wPos+lb).M() - lepT_m_);
  const double tLepDmNeg = abs((wNeg+lb).M() - lepT_m_);
  const double nupz = (tLepDmPos < tLepDmNeg ) ? a+b : a-b;

  nu.SetXYZM(nu.Px(), nu.Py(), nupz, 0);
  return nu;
}

std::vector<TLorentzVector> TTLJKinFitFtn::getSolution(const double* x) const
{
  const double scale_jet = x[0];
  std::vector<TLorentzVector> out;
  out.push_back(lvs_[0]); // MET
  out.push_back(lvs_[1]); // Lepton
  double dx = 0, dy = 0;
  for ( size_t i=2; i<6; ++i ) {
    out.push_back(lvs_[i]*scale_jet);
    dx += out[i].X() - lvs_[i].X();
    dy += out[i].Y() - lvs_[i].Y();
  }
  out[0].SetXYZM(out[0].Px()-dx, out[0].Py()-dy, 0, 0); // MET correction by JEC FIXME: this is an approximation, does not consider additional jets
  out[0] = solveLepTopNu(out[0], out[1], out[2]);

  return out;
}

TTLJKinFit::TTLJKinFit(const double hadW_m, const double hadT_m, const double lepW_m, const double lepT_m):
  ftn_(hadW_m, hadT_m, lepW_m, lepT_m)
{
  min_.reset(new ROOT::Minuit2::Minuit2Minimizer(ROOT::Minuit2::kMigrad));
  //min_.reset(new ROOT::Minuit2::Minuit2Minimizer(ROOT::Minuit2::kCombined));
  //min_.reset(new ROOT::Minuit2::Minuit2Minimizer(ROOT::Minuit2::kScan));
  min_->SetVariable(0, "jes", 1, 0.01);
  //min_->SetLimitedVariable(0, "jes", 1, 0.1, 0.5, 1.5);
  min_->SetFunction(ftn_);
}

std::vector<TLorentzVector> TTLJKinFit::getSolution() const
{
  return ftn_.getSolution(min_->X());
}

double TTLJKinFit::compute(const TLorentzVector metP4, const TLorentzVector l1, const TLorentzVector b1,
                           const TLorentzVector wj1, const TLorentzVector wj2, const TLorentzVector hb)
{
  ftn_.lvs_[0] = metP4;
  ftn_.lvs_[1] = l1;
  ftn_.lvs_[2] = b1;
  ftn_.lvs_[3] = wj1;
  ftn_.lvs_[4] = wj2;
  ftn_.lvs_[5] = hb;

  min_->SetVariableValue(0, 1.0);
  //min_->SetPrecision(1e-2);
  bool isValid = false;
  for ( int i=0; i<4; ++i ) {
    isValid = min_->Minimize();
    if ( isValid ) break;
  }
  chi2_ = !isValid ? 1e9 : ftn_.DoEval(min_->X());

  return chi2_;
}
