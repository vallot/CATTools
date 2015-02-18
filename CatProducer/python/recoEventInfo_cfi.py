import FWCore.ParameterSet.Config as cms

recoEventInfo = cms.EDProducer("RecoEventInfoProducer",
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    HLT = cms.PSet(
        ## double Mu
        HLT_Mu17_Mu8 = cms.vstring("HLT_Mu17_Mu8_v*"),
        HLT_Mu17_TkMu8 = cms.vstring("HLT_Mu17_TkMu8_v*"),
        HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsonVVL = cms.vstring("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsonVVL_v*"),
        ## double El
        HLT_Ele23_Ele12_CaloId_TrackId_Iso = cms.vstring("HLT_Ele23_Ele12_CaloId_TrackId_Iso_v*"),
        ## mu+el
        HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP = cms.vstring("HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v*"),
        HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP = cms.vstring("HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v*"),
        ## 8TeV
        HLT_IsoMu20_eta2p1 = cms.vstring("HLT_IsoMu20_eta2p1_v*"),
        HLT_Ele27_WP80 = cms.vstring("HLT_Ele27_WP80_v*"),
        ## 5.3E33 [50ns]
        HLT_IsoMu24_eta2p1 = cms.vstring("HLT_IsoMu24_eta2p1_v*"),
        HLT_Ele27_WP75_eta2p1 = cms.vstring("HLT_Ele27_WP75_eta2p1_v*"),
        ## 7E33 [25ns] OOT+in time <mu> = 20
        HLT_IsoMu20_eta2p1_IterTrk02 = cms.vstring("HLT_IsoMu20_eta2p1_IterTrk02_v*"),
        ## 1.4E34 [25ns]
        HLT_IsoMu24_eta2p1_IterTrk02 = cms.vstring("HLT_IsoMu24_eta2p1_IterTrk02_v*"),
        HLT_Ele32_WP75_eta2p1 = cms.vstring("HLT_Ele32_WP75_eta2p1_v*"),
        HLT_Ele32_WP75_Gsf = cms.vstring("HLT_Ele32_WP75_Gsf_v*"),
    ),
)

