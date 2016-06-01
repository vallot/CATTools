#ifndef TopTriggerSF_H
#define TopTriggerSF_H
#include "CATTools/CommonTools/interface/TTbarModeDefs.h"

// Scale factors from AN2015/309
double computeTrigSF(const int channel, const double eta1, const double eta2, int direction=0)
{
  const double aeta1 = std::abs(eta1);
  const double aeta2 = std::abs(eta2);

  if ( channel == cat::CH_ELEL ) {
    if ( aeta1 < 1.2 ) {
      if      ( aeta2 < 1.2 ) return 0.953 + direction*0.023;
      else if ( aeta2 < 2.4 ) return 0.957 + 0.5*direction*0.027;
    }
    else if ( aeta1 < 2.4 ) {
      if      ( aeta2 < 1.2 ) return 0.967 + 0.5*direction*0.029;
      else if ( aeta2 < 2.4 ) return 0.989 - direction*0.031;
    }
  }
  else if ( channel == cat::CH_MUMU ) {
    if ( aeta1 < 1.2 ) {
      if      ( aeta2 < 1.2 ) return 0.926 + direction*0.022;
      else if ( aeta2 < 2.4 ) return 0.943 + 0.5*direction*0.026;
    }
    else if ( aeta1 < 2.4 ) {
      if      ( aeta2 < 1.2 ) return 0.958 + 0.5*direction*0.027;
      else if ( aeta2 < 2.4 ) return 0.926 - direction*0.027;
    }
  }
  else if ( channel == cat::CH_MUEL ) {
    if ( aeta1 < 1.2 ) {
      if      ( aeta2 < 1.2 ) return 0.966 + direction*0.022;
      else if ( aeta2 < 2.4 ) return 0.973 + 0.5*direction*0.025;
    }
    else if ( aeta1 < 2.4 ) {
      if      ( aeta2 < 1.2 ) return 0.981 + 0.5*direction*0.025;
      else if ( aeta2 < 2.4 ) return 0.984 - direction*0.030;
    }
  }

  return 1;
};

// Scale factors are from AN16-025 (v4) http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2016_025_v4.pdf
/*
double computeTrigSF(const int channel, const double eta1, const double eta2, int direction=0)
{
  if      ( channel == cat::CH_ELEL ) return 0.953 + direction*0.009;
  else if ( channel == cat::CH_MUMU ) return 0.948 + direction*0.002;
  else if ( channel == cat::CH_MUEL ) return 0.975 + direction*0.004;
}
*/

#endif
