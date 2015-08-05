import FWCore.ParameterSet.Config as cms

catTriggerInfo = cms.EDProducer("CATTriggerInfoProducer",
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    prescales = cms.InputTag("patTrigger"),
    PreScaled = cms.vstring(),
    unPreScaled = cms.vstring(
    "HLT_Mu17_Mu8_DZ_v*","HLT_Mu17_TkMu8_DZ_v*","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*","HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*","HLT_Ele12_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*", "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*"
    )
)
