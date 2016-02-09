#ifndef BTagScaleFactorEvaluaor_H
#define BTagScaleFactorEvaluaor_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CATTools/CatAnalyzer/interface/BTagCalibrationStandalone.h"
#include "CATTools/DataFormats/interface/Jet.h"

#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TH1D.h"

namespace cat {

class BTagWeightEvaluator
{
public:
  BTagWeightEvaluator() {};
  void init(const edm::ParameterSet& pset);
  double computeWeight(const cat::JetCollection& jets, const int unc) const;

  // For per-jet SF evaluation
  double getSF(const cat::Jet& jet, const int unc) const {
    if ( method_ == CSVWEIGHT ) return computeCSVWeightFromROOT(jet, unc);
    return computeCSVWeightFromCSV(jet, unc);
  };

private:
  int method_;
  std::string btagAlgo_;
  std::vector<std::string> uncNames_;
  std::map<int, BTagCalibrationReader> readers_;

  // For the CSV weights
  double computeCSVWeightFromCSV(const cat::Jet& jet, const int unc) const;
  double computeCSVWeightFromROOT(const cat::Jet& jet, const int unc) const;
  void fillCSVHistos(TFile *fileHF, TFile *fileLF, int nHFptBins);

  enum METHOD {
    INCL, MUJET, ITERATIVEFIT, CSVWEIGHT
  };

  enum CSVUNC {
    CENTRAL,
    JES_UP, JES_DN,
    LF_UP, LF_DN,
    HF_UP, HF_DN,
    HFSTAT1_UP, HFSTAT1_DN, HFSTAT2_UP, HFSTAT2_DN,
    LFSTAT1_UP, LFSTAT1_DN, LFSTAT2_UP, LFSTAT2_DN,
    CFERR1_UP, CFERR1_DN, CFERR2_UP, CFERR2_DN
  };

  TH1D *h_csv_wgt_hf[9][6];
  TH1D *c_csv_wgt_hf[9][6];
  TH1D *h_csv_wgt_lf[9][4][3];
  int nHFptBins_;
};

}

#endif
