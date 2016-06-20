#ifndef TopTriggerSF_H
#define TopTriggerSF_H

#include "CATTools/DataFormats/interface/Lepton.h"

double computeTrigSF(const cat::Lepton& lep, int direction=0);
double computeTrigSF(const cat::Lepton& lep1, const cat::Lepton& lep2, int direction=0);
double computeTrigSFInclusive(const cat::Lepton& lep1, const cat::Lepton& lep2, int direction=0);

#endif
