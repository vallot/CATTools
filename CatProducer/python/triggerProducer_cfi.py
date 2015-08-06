import FWCore.ParameterSet.Config as cms

catTrigger = cms.EDProducer("CATTriggerProducer",
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    prescales = cms.InputTag("patTrigger"),
    PreScaled = cms.vstring("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v*","HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v*","HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*"),
    unPreScaled = cms.vstring(
    "HLT_Mu17_Mu8_DZ_v*","HLT_Mu17_TkMu8_DZ_v*","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*","HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*","HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*","HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*","HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet30_v*"
    )
)
