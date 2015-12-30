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
    pt_bins = cms.vdouble(15, 25, 35, 45, 55, 1e9),
    eta_bins = cms.vdouble(-2.5, -1.5, -1.0, 0, 1.0, 1.5, 2.5),
    values = cms.vdouble(
        0.96, 0.95, 0.98, 0.99, 0.99, 0.97,
        0.98, 0.97, 0.97, 0.99, 0.99, 0.98,
        0.98, 0.99, 0.99, 0.99, 0.99, 0.98,
        0.98, 0.99, 0.99, 0.99, 0.99, 0.99,
        0.99, 0.99, 0.99, 1.00, 0.99, 0.99,
    ),
    errors = cms.vdouble(
        0.02, 0.06, 0.01, 0.01, 0.02, 0.01,
        0.01, 0.01, 0.02, 0.01, 0.01, 0.01,
        0.00, 0.00, 0.00, 0.00, 0.00, 0.01,
        0.01, 0.00, 0.01, 0.00, 0.01, 0.01,
        0.01, 0.02, 0.01, 0.00, 0.02, 0.02,
    ),
)

