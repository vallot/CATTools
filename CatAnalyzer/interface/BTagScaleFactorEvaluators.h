#ifndef BTagScaleFactorEvaluaor_H
#define BTagScaleFactorEvaluaor_H

#include <vector>
#include <utility>
#include <algorithm>

namespace cat {

class CSVWeightEvaluator
{
public:
  CSVWeightEvaluator();
  double operator()(const cat::Jet& jet, const int unc) const;

  enum UNC {
    CENTRAL,
    JES_UP, JES_DN,
    LF_UP, LF_DN,
    HF_UP, HF_DN,
    HFSTATS1_UP, HFSTATS1_DN, HFSTATS2_UP, HFSTATS2_DN, 
    LFSTATS1_UP, LFSTATS1_DN, LFSTATS2_UP, LFSTATS2_DN, 
    CFERR1_UP, CFERR1_DN, CFERR2_UP, CFERR2_DN;
  };
  
private:
  std::vector<int, BTagCalibrationReader> readers_;
};

}

#endif
