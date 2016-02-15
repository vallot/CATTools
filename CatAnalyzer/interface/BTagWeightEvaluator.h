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
  void init(const int combineMethod,
            const std::string btagName, BTagEntry::OperatingPoint operationPoint,
            const int minNbjet);
  void initCSVWeight(const bool isFromROOT, const std::string btagName);

  double eventWeight(const cat::JetCollection& jets, const int unc) const;

  // For per-jet SF evaluation - useful for CSV weight
  double getSF(const cat::Jet& jet, const int unc) const;

  enum CSVUNC {
    CENTRAL, JES_UP, JES_DN,
    LF_UP, LF_DN, HF_UP, HF_DN,
    HFSTAT1_UP, HFSTAT1_DN, HFSTAT2_UP, HFSTAT2_DN,
    LFSTAT1_UP, LFSTAT1_DN, LFSTAT2_UP, LFSTAT2_DN,
    CFERR1_UP, CFERR1_DN, CFERR2_UP, CFERR2_DN
  };
  enum TYPE {STANDARD, ITERATIVEFIT, CSVWEIGHT};

private:
  int type_; // Measurement type
  int method_; // Combination method

  std::string btagAlgo_;
  std::vector<std::string> uncNames_;
  std::map<int, BTagCalibrationReader> readers_;

  // For the standard methods
  int minNbjet_;

  // For the iterative fit method
  double getCSVWeightSFFromROOT(const double pt, const double aeta,
                                const int flav, const double discr,
                                const int unc) const;
  void fillCSVHistos(TFile *fileHF, TFile *fileLF, int nHFptBins);

  TH1D *h_csv_wgt_hf[9][6];
  TH1D *c_csv_wgt_hf[9][6];
  TH1D *h_csv_wgt_lf[9][4][3];
  int nHFptBins_;
};

}

#endif
