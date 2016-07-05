#ifndef __CATTools_TopEventGlobalVar__
#define __CATTools_TopEventGlobalVar__
namespace TopEventCommonGlobal {
  enum sys_e {sys_nom,
    sys_jes_u, sys_jes_d, sys_jer_u, sys_jer_d,
    sys_mu_u, sys_mu_d, sys_el_u, sys_el_d,
    nsys_e
  };

  static std::vector<std::string> sys_name = {
    "nom",
    "jes_u", "jes_d", "jer_u", "jer_d",
    "mu_u", "mu_d", "el_u", "el_d",
  };
  static const int NCutflow = sys_name.size();
  typedef std::vector<const cat::Lepton*> LeptonPtrs;
  typedef math::XYZTLorentzVector LV;
  typedef std::vector<LV> VLV;
}
#endif
