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
    triggerResults = cms.InputTag("TriggerResults::PAT"),
    secondaryTriggerResults = cms.InputTag("TriggerResults::RECO"),
    triggerPrescales = cms.InputTag("patTrigger"),
    combineBy = cms.string("and"),
    triggersToMatch = cms.vstring(
        "HBHENoiseFilter",
        "HBHENoiseIsoFilter",
        "EcalDeadCellTriggerPrimitiveFilter",
        "goodVertices",
        "eeBadScFilter",
        "globalTightHalo2016Filter",
    ),
    doFilter = cms.bool(False),
)

filterTrigMUEL = cms.EDFilter("CATTriggerBitCombiner",
    triggerResults = cms.InputTag("TriggerResults::HLT2"),
    secondaryTriggerResults = cms.InputTag("TriggerResults::HLT"),
    triggerPrescales = cms.InputTag("patTrigger"),
    combineBy = cms.string("or"),
    triggersToMatch = cms.vstring(
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
        "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
        "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
    ),
    doFilter = cms.bool(False),
)

filterTrigMuELMC = filterTrigMUEL.clone(triggersToMatch = cms.vstring(
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
