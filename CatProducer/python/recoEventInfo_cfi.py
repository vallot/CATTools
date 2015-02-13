import FWCore.ParameterSet.Config as cms

recoEventInfo = cms.EDProducer("RecoEventInfoProducer",
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    HLT = cms.PSet(
        DoubleMu = cms.vstring(
            "HLT_Mu17_Mu8_v*", "HLT_Mu17_TkMu8_v*","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsonVVL_v*"
        ),
        DoubleEl = cms.vstring(
            "HLT_Ele23_Ele12_CaloId_TrackId_Iso_v*"
        ),
        MuEl = cms.vstring(
            "HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v*","HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v*"
        ),
        MuJet = cms.vstring(
            "HLT_IsoMu24_eta2p1_v*"
        ),
        ELJet = cms.vstring(
            "HLT_Ele27_WP80_v*"
        ),
    ),
)

