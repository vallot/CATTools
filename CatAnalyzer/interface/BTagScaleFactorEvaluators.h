#ifndef BTagScaleFactorEvaluaor_H
#define BTagScaleFactorEvaluaor_H

#include "CATTools/CatAnalyzer/interface/BTagCalibrationStandalone.h"
#include "CATTools/DataFormats/interface/Jet.h"

// original code : https://github.com/cms-ttH/MiniAOD/blob/master/MiniAODHelper/interface/CSVHelper.h
#include "TFile.h"
#include "TH1D.h"

#include <vector>
#include <utility>
#include <algorithm>

namespace cat {

class CSVWeightEvaluator
{
public:
  CSVWeightEvaluator(const int inputType = ROOT);

  double operator()(const cat::Jet& jet, const int unc) const;

  enum InputType { ROOT, CSV };

  enum UNC {
    CENTRAL,
    JES_UP, JES_DN,
    LF_UP, LF_DN,
    HF_UP, HF_DN,
    HFSTAT1_UP, HFSTAT1_DN, HFSTAT2_UP, HFSTAT2_DN,
    LFSTAT1_UP, LFSTAT1_DN, LFSTAT2_UP, LFSTAT2_DN,
    CFERR1_UP, CFERR1_DN, CFERR2_UP, CFERR2_DN
  };
  
private:
  // For the CSV reader
  double computeWeightFromCSV(const cat::Jet& jet, const int unc) const;
  std::map<int, BTagCalibrationReader> readers_;

  // For the ROOT file reader
  double computeWeightFromROOT(const cat::Jet& jet, const int unc) const;
  void fillCSVHistos(TFile *fileHF, TFile *fileLF, int nHFptBins);

  // CSV reweighting
  TH1D *h_csv_wgt_hf[9][6];
  TH1D *c_csv_wgt_hf[9][6];
  TH1D *h_csv_wgt_lf[9][4][3];
  int nHFptBins_;
};

}

#endif
