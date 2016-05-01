import FWCore.ParameterSet.Config as cms

## Muon SF reference https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffsRun2
## Retrieve data from the cmsdoc web page:
##   https://cmsdoc.cern.ch/cms/Physics/muon/ReferenceEfficiencies/Run2015/25ns/MuonID_Z_RunCD_Reco76X_Feb15.json
##   https://cmsdoc.cern.ch/cms/Physics/muon/ReferenceEfficiencies/Run2015/25ns/MuonIso_Z_RunCD_Reco76X_Feb15.json
muonSFTight = cms.PSet(
    # Values of "MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1 + MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1"
    pt_bins = cms.vdouble(20, 25, 30, 40, 50, 60, 120),
    abseta_bins = cms.vdouble(0, 0.9, 1.2, 2.1, 2.4),
    values = cms.vdouble(
        0.975283, 0.980129, 0.986676, 0.985789, 0.980877, 0.985445,
        0.984484, 0.970783, 0.980673, 0.980057, 0.977943, 0.984907,
        0.986784, 0.990845, 0.993367, 0.990324, 0.988209, 0.991915,
        0.971296, 0.971883, 0.976009, 0.977, 0.973655, 0.963135,
    ),
    errors = cms.vdouble(
        0.0157468, 0.0136165, 0.0120878, 0.0118247, 0.0125879, 0.0135921,
        0.0181504, 0.0154891, 0.0127943, 0.0121561, 0.0134922, 0.0152142,
        0.0145121, 0.0132486, 0.0120305, 0.0116593, 0.0124603, 0.0137997,
        0.0176168, 0.0153647, 0.0131657, 0.0126734, 0.0148934, 0.0197265,
    ),

)

## Electorn SF from http://fcouderc.web.cern.ch/fcouderc/EGamma/scaleFactors/moriond2016_76X/eleCutBasedID/CutBasedID_MediumWP_76X_18Feb.txt_egammaPlots.pdf
electronSFCutBasedIDMediumWP = cms.PSet(
    # Values of "CutBasedID_MediumWP"
    pt_bins = cms.vdouble(10, 20, 30, 40, 50, 200),
    abseta_bins = cms.vdouble(0, 0.8, 1.444, 1.566, 2, 2.5),
    values = cms.vdouble(
        1.00726, 0.971182, 0.984906, 0.985899, 0.988598,
        1.0902, 0.983359, 0.987179, 0.98692, 0.986159,
        1.08642, 0.963054, 0.949123, 0.981612, 0.997257,
        0.984444, 0.936809, 0.975066, 0.992806, 1.00579,
        1.03557, 0.986446, 0.963351, 1.00611, 1.00949,
    ),
    errors = cms.vdouble(
        0.0656124, 0.0279776, 0.0291216, 0.0175478, 0.00603364,
        0.0391175, 0.0256741, 0.0171527, 0.0277844, 0.00540994,
        0.0800687, 0.261884 , 0.0511487, 0.036227 , 0.113466  ,
        0.0447765, 0.0186922, 0.0274652, 0.0404665, 0.0125727 ,
        0.0303603, 0.0263447, 0.0445603, 0.0103732, 0.0137829 ,
    ),
)

## Retrieve data from the cmsdoc web page:
##   https://cmsdoc.cern.ch/cms/Physics/muon/ReferenceEfficiencies/Run2015/25ns/MuonID_Z_RunCD_Reco74X_Dec1.json
muonSFTight74X = cms.PSet(
    # Values of "NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1 + NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1"
    pt_bins = cms.vdouble(20, 25, 30, 40, 50, 60, 120),
    abseta_bins = cms.vdouble(0, 0.9, 1.2, 2.1, 2.4),
    values = cms.vdouble(
        0.979479, 0.984375, 0.987068, 0.986086, 0.983674, 0.985021,
        0.978054, 0.980927, 0.984458, 0.981568, 0.975983, 0.984876,
        0.99541, 0.991193, 0.994704, 0.992697, 0.988925, 0.995536,
        0.983608, 0.974278, 0.976037, 0.980598, 0.969666, 0.969941,
    ),
    errors = cms.vdouble(
        0.0191201, 0.0168928, 0.0151444, 0.0148111, 0.015527, 0.0162596,
        0.021688, 0.0189092, 0.0159175, 0.0150692, 0.0164092, 0.0178394,
        0.0179975, 0.0166027, 0.0151404, 0.0147407, 0.0154785, 0.0167413,
        0.0219747, 0.0194149, 0.0165014, 0.0159522, 0.0180924, 0.0243999,
    ),
)

electronSFWP9074X = cms.PSet(
    ## SF from https://indico.cern.ch/event/369259/contribution/3/attachments/1204731/1755463/update_SFs_triggering_MVA_ID_systematics.pdf
    eta_bins = cms.vdouble(-2.5, -1.5, -1.0, 0, 1.0, 1.5, 2.5),
    pt_bins = cms.vdouble(15, 25, 35, 45, 55, 1e9),
    values = cms.vdouble(
        0.96, 0.95, 0.98, 0.99, 0.99, 0.97,
        0.98, 0.97, 0.97, 0.99, 0.99, 0.98,
        0.98, 0.99, 0.99, 0.99, 0.99, 0.98,
        0.98, 0.99, 0.99, 0.99, 0.99, 0.99,
        0.99, 0.99, 1.00, 1.00, 0.99, 0.99,
    ),
    errors = cms.vdouble(
        0.02, 0.06, 0.01, 0.01, 0.02, 0.01,
        0.01, 0.01, 0.02, 0.01, 0.01, 0.01,
        0.00, 0.00, 0.00, 0.00, 0.00, 0.01,
        0.01, 0.00, 0.01, 0.00, 0.01, 0.01,
        0.01, 0.02, 0.01, 0.00, 0.02, 0.02,
    ),
)

electronSFCutBasedIDMediumWP74X = cms.PSet(
    ## SF from https://indico.cern.ch/event/491528/contribution/2/attachments/1231399/1805319/CutBasedId_EGM_19Feb.pdf
    ## Actual numbers are taken from https://arun-afs.web.cern.ch/arun-afs/Final_Fits_Data_74X/CutBasedID_MediumWP_fromTemplates_withSyst_Final.txt
    # Values of "CutBasedID_MediumWP"
    pt_bins = cms.vdouble(10, 20, 30, 40, 50, 200),
    eta_bins = cms.vdouble(-2.5, -2, -1.566, -1.444, -0.8, 0, 0.8, 1.444, 1.566, 2, 2.5),
    values = cms.vdouble(
        1.25131, 1.0536, 1.01235, 1.00875, 1.006,
        1.09254, 0.969453, 0.971279, 0.983313, 0.985109,
        1.29412, 1.04043, 1.02482, 0.990155, 0.948856,
        1.32673, 1.05939, 0.990897, 0.974057, 0.958097,
        1.23371, 1.02262, 0.978616, 0.969838, 0.967489,
        1.22889, 1.01799, 0.986164, 0.974448, 0.970819,
        1.30827, 1.07742, 1.00654, 0.976359, 0.964773,
        1.35498, 1.18801, 0.990991, 0.956153, 0.974394,
        1.06283, 0.980392, 0.961892, 0.978469, 0.994253,
        1.20769, 1.02649, 1.0096, 1.00125, 1.00481,
    ),
    errors = cms.vdouble(
        0.712863, 0.0517635, 0.0242599, 0.0222102, 0.0337348,
        0.662784, 0.0757762, 0.575897, 0.0392023, 0.0671983,
        3.75487, 1.17119, 0.124767, 0.0150145, 0.0555887,
        1.2277, 0.20447, 0.0280944, 0.00855398, 0.0355446,
        0.813571, 0.108281, 0.00853765, 0.0217419, 0.0302874,
        0.761126, 0.580145, 0.0313306, 0.0196618, 0.0253676,
        0.874693, 0.237932, 0.650881, 0.0103168, 0.0286228,
        4.59432, 1.84669, 0.1764, 0.0506619, 0.161535,
        0.797724, 0.0632952, 0.0377029, 0.0313089, 0.0300417,
        0.617985, 0.0693681, 0.0293348, 0.024152, 0.0894748,
    ),
)

