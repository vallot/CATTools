import FWCore.ParameterSet.Config as cms

filterLumi = cms.EDFilter("EventWeightThresholdFilter",
    isLessThan = cms.bool(False),
    threshold = cms.double(0.5), # Require >= 0.5 to select events with bit=1
    type = cms.string("int"),
    src = cms.InputTag("lumiMask"),
)

filterLumiSilver = filterLumi.clone(src = cms.InputTag("lumiMaskSilver"))

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
filterRECO = cms.EDFilter("CATTriggerBitCombiner",
    src = cms.InputTag("catTrigger"),
    combineBy = cms.string("and"),
    triggersToMatch = cms.vstring(
        "Flag_HBHENoiseFilter",
        "Flag_HBHENoiseIsoFilter",
        "Flag_EcalDeadCellTriggerPrimitiveFilter",
        "Flag_goodVertices",
        "Flag_eeBadScFilter",
        "Flag_globalTightHalo2016Filter",

        "Flag_badPFMuon",
        "Flag_badChargedHadron",
    ),
    doFilter = cms.bool(False),
)

filterRECOMC = filterRECO.clone(
    triggersToMatch = cms.vstring(
        "Flag_HBHENoiseFilter",
        "Flag_HBHENoiseIsoFilter",
        "Flag_goodVertices",
        "Flag_eeBadScFilter",
        "Flag_globalTightHalo2016Filter",

        "Flag_badPFMuon",
        "Flag_badChargedHadron",
    ),
)

filterTrigMUEL = cms.EDFilter("CATTriggerBitCombiner",
    src = cms.InputTag("catTrigger"),
    combineBy = cms.string("or"),
    triggersToMatch = cms.vstring(
#        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
#        "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
#        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
#        "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4 OR HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v4"
    ),
    doFilter = cms.bool(False),
)

filterTrigMUELMC = filterTrigMUEL.clone(triggersToMatch = cms.vstring(
    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",))

filterTrigELEL = filterTrigMUEL.clone(triggersToMatch = cms.vstring(
    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",))

filterTrigMUMU = filterTrigMUEL.clone(triggersToMatch = cms.vstring(
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",))

filterTrigMU = filterTrigMUEL.clone(triggersToMatch = cms.vstring(
    "HLT_IsoMu24_v", "HLT_IsoTkMu24_v",))

filterTrigEL = filterTrigMUEL.clone(triggersToMatch = cms.vstring(
    "HLT_Ele32_eta2p1_WPTight_Gsf_v",))

removeLumisWithL1TCert = cms.EDFilter("LumiMaskFilter",
    LumiSections = cms.untracked.VLuminosityBlockRange(
        ## LS to apply L1T certification
        "274241:297", "274315:376", "274388:1731", "274442:257", "274969:481-274969:484",
        "275067:234", "275375:843-275375:864",
    ),
    acceptOnFail = cms.untracked.bool(True),
    doFilter = cms.untracked.bool(True)
)
