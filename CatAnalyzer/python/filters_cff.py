import FWCore.ParameterSet.Config as cms

filterLumi = cms.EDFilter("EventWeightThresholdFilter",
    isLessThan = cms.bool(False),
    threshold = cms.double(0.5), # Require >= 0.5 to select events with bit=1
    type = cms.string("int"),
    src = cms.InputTag("lumiMask"),
)

filterLumiSilver = filterLumi.clone(src = cms.InputTag("lumiMaskSilver"))

filterRECO = cms.EDFilter("CATTriggerBitCombiner",
    triggerResults = cms.InputTag("TriggerResults::PAT"),
    secondaryTriggerResults = cms.InputTag("TriggerResults::RECO"),
    triggerPrescales = cms.InputTag("patTrigger"),
    combineBy = cms.string("and"),
    triggersToMatch = cms.vstring(
      "HBHENoiseFilter",
      "HBHENoiseIsoFilter",
      "CSCTightHalo2015Filter",
      "EcalDeadCellTriggerPrimitiveFilter",
      "goodVertices",
      "eeBadScFilter"
    ),
    doFilter = cms.bool(False),
)

filterTrigMUEL = cms.EDFilter("CATTriggerBitCombiner",
    triggerResults = cms.InputTag("TriggerResults::HLT"),
    triggerPrescales = cms.InputTag("patTrigger"),
    combineBy = cms.string("or"),
    triggersToMatch = cms.vstring(
        "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
        "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v",
    ),
    doFilter = cms.bool(False),
)

filterTrigELEL = filterTrigMUEL.clone(
    triggersToMatch = cms.vstring(
        "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
    ),
)

filterTrigMUMU = filterTrigMUEL.clone(
    triggersToMatch = cms.vstring(
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
    ),
)
filterTrigMUJET = filterTrigMUEL.clone(
    triggersToMatch = cms.vstring(
        "HLT_IsoMu18_v",
    ),
)
filterTrigELJET = filterTrigMUEL.clone(
    triggersToMatch = cms.vstring(
        "HLT_Ele23_WPLoose_Gsf_v",
    ),
)


removeLumisWithBadBS = cms.EDFilter("LumiMaskFilter",
    LumiSections = cms.untracked.VLuminosityBlockRange(
        ## LS with bad beamspot problems
        ## which was announced from the certification HN https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2556.html
        ## JSON file was taken from /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/BeamSpotIssue_JSON.txt
        "259626:83-259626:106", "259626:108-259626:111", "259626:115-259626:166", "259626:169-259626:215", "259626:218-259626:437", 
        "259637:1-259637:221", "259681:64-259681:98", "259682:1-259682:4", "259683:3-259683:19", "259683:22-259683:23", 
        "259683:25-259683:94", "259685:1-259685:209", "259685:213-259685:240", "259685:242-259685:290", "259685:292-259685:544", 
        "259685:546-259685:630"),
    acceptOnFail = cms.untracked.bool(True),
    doFilter = cms.untracked.bool(True)
)

removeUncheckedLumis765Prod2015 = cms.EDFilter("LumiMaskFilter",
    LumiSections = cms.untracked.VLuminosityBlockRange(
        ## This is a temporary filter since the cat765 prod was run with old JSON file v1.
        ## {"256801": [[73, 74]], "256842": [[131, 132]]} are to be removed to include HO DCS information
        ## https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2611.html
        "256801:73-256801:74", "256842:131-256842:132"),
    acceptOnFail = cms.untracked.bool(True),
    doFilter = cms.untracked.bool(True)
)
