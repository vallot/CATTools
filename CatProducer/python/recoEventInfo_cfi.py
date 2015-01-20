import FWCore.ParameterSet.Config as cms

recoEventInfo = cms.EDProducer("RecoEventInfoProducer",
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    HLT = cms.PSet(
        DoubleMu = cms.vstring(
            "HLT_Mu17_Mu8_v*", "HLT_Mu17_TkMu8_v*"
        ),
    ),
)

