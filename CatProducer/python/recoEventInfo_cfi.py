import FWCore.ParameterSet.Config as cms

recoEventInfo = cms.EDProducer("RecoEventInfoProducer",
    vertex = cms.InputTag("goodOfflinePrimaryVertices"),
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    HLT = cms.PSet(
        DoubleMu = cms.vstring(
            "HLT_Mu17_Mu8_v*", "HLT_Mu17_TkMu8_v*"
        ),
        DoubleElectron = cms.vstring(
            "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
        ),
        MuEG = cms.vstring(
            "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
            "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
        ),
        MuJet = cms.vstring(
            "HLT_IsoMu24_eta2p1_v*",
        ),
        ElJet = cms.vstring(
            "HLT_Ele25_CaloIdVL_*", "HLT_Ele25_CaloIdVT_*",
        ),
    ),
)

