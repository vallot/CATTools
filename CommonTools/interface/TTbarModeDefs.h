#ifndef CATTools_CommonTools_TTbarModeDefs_H
#define CATTools_CommonTools_TTbarModeDefs_H

namespace cat
{

enum TTbarMode { CH_NONE = -1, CH_FULLHADRON = 0, CH_SEMILEPTON, CH_FULLLEPTON };
enum DecayMode { CH_HADRON = 0, CH_MUON, CH_ELECTRON,
                 CH_TAU_HADRON, CH_TAU_MUON, CH_TAU_ELECTRON };

};

#endif

