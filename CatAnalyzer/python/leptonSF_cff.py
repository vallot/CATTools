import FWCore.ParameterSet.Config as cms

## Muon SF reference https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffsRun2
## Retrieve data from the cmsdoc web page:
##   https://cmsdoc.cern.ch/cms/Physics/muon/ReferenceEfficiencies/Run2015/25ns/MuonID_Z_RunCD_Reco74X_Dec1.json
muonSFTight = cms.PSet(
    # Values of "NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1 + NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1"
    abseta_bins = cms.vdouble(0, 0.9, 1.2, 2.1, 2.4),
    pt_bins = cms.vdouble(20, 25, 30, 40, 50, 60, 120),
    values = cms.vdouble(
        0.979479, 0.978054, 0.99541, 0.983608, 
        0.984375, 0.980927, 0.991193, 0.974278, 
        0.987068, 0.984458, 0.994704, 0.976037, 
        0.986086, 0.981568, 0.992697, 0.980598, 
        0.983674, 0.975983, 0.988925, 0.969666, 
        0.985021, 0.984876, 0.995536, 0.969941, 
    ),
    errors = cms.vdouble(
        0.0191201, 0.021688, 0.0179975, 0.0219747, 
        0.0168928, 0.0189092, 0.0166027, 0.0194149, 
        0.0151444, 0.0159175, 0.0151404, 0.0165014, 
        0.0148111, 0.0150692, 0.0147407, 0.0159522, 
        0.015527, 0.0164092, 0.0154785, 0.0180924, 
        0.0162596, 0.0178394, 0.0167413, 0.0243999, 
    ),
)

electronSFWP90 = cms.PSet(
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

electronSFCutMedium = cms.PSet(
    ## SF from https://indico.cern.ch/event/369259/contribution/1/attachments/1204682/1754936/EGM_14Dec.pdf
    ## Actual numbers are taken from /afs/cern.ch/work/a/arun/public/forEGM/CutBased_textFiles_withSyst/CutBasedID_MediumWP_fromTemplates_withSyst_v1.txt
    ## Documentations are available in https://twiki.cern.ch/twiki/bin/viewauth/CMS/HWW2015TriggerResults
    ## FIXME : NUMBERS SEEM TO BE PRELIMINARY, SCALE FACTORS ARE GREATER THAN 1
    ## FIXME : I'M SETTING UNCERTAINTY TO ZERO!!! NOT EASY TO GUESS FROM THE UNDOCUMENTED TABLE
    pt_bins = cms.vdouble(10.000, 20.000, 30.000, 40.000, 50.000, 200.000),
    scEta_bins = cms.vdouble(-2.500, -2.000, -1.566, -1.444, -0.800, 0.000, 0.800, 1.444, 1.566, 2.000, 2.500),
    values = cms.vdouble(
        0.478/0.382, 0.629/0.597, 0.738/0.729, 0.807/0.800, 0.838/0.833,
        0.425/0.389, 0.603/0.622, 0.744/0.766, 0.825/0.839, 0.860/0.873,
        0.308/0.238, 0.386/0.371, 0.578/0.564, 0.704/0.711, 0.705/0.743,
        0.536/0.404, 0.660/0.623, 0.762/0.769, 0.826/0.848, 0.846/0.883,
        0.549/0.445, 0.678/0.663, 0.778/0.795, 0.836/0.862, 0.863/0.892,
        0.553/0.450, 0.679/0.667, 0.784/0.795, 0.839/0.861, 0.865/0.891,
        0.522/0.399, 0.668/0.620, 0.770/0.765, 0.826/0.846, 0.849/0.880,
        0.313/0.231, 0.436/0.367, 0.550/0.555, 0.676/0.707, 0.723/0.742,
        0.406/0.382, 0.600/0.612, 0.732/0.761, 0.818/0.836, 0.865/0.870,
        0.471/0.390, 0.620/0.604, 0.736/0.729, 0.799/0.798, 0.835/0.831,
    ),
    errors = cms.vdouble(
        0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0.,
    ),
)
