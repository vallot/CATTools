#ifndef CATTools_CommonTools_TTbarModeDefs_H
#define CATTools_CommonTools_TTbarModeDefs_H

namespace cat
{

enum TTChannel { CH_NOTT = -1, CH_FULLHADRON = 0, CH_SEMILEPTON, CH_FULLLEPTON };
enum WMode { CH_HADRON = 0, CH_MUON, CH_ELECTRON,
             CH_TAU_HADRON, CH_TAU_MUON, CH_TAU_ELECTRON };
enum TTLLChannel { CH_NOLL = 0, CH_MUEL, CH_ELEL, CH_MUMU };
enum TTLJChannel { CH_NOLEP = 0, CH_MUJET, CH_ELJET };

};

#endif

