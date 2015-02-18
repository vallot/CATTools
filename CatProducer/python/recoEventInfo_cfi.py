import FWCore.ParameterSet.Config as cms

recoEventInfo = cms.EDProducer("RecoEventInfoProducer",
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    HLT = cms.PSet(
        ## double Mu
        Mu17Mu8 = cms.vstring("HLT_Mu17_Mu8_v*"),
        Mu17TkMu8 = cms.vstring("HLT_Mu17_TkMu8_v*"),
        Mu17TrkIsoVVLTkMu8TrkIsonVVL = cms.vstring("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsonVVL_v*"),
        ## double El
        Ele23Ele12CaloIdTrackIdIso = cms.vstring("HLT_Ele23_Ele12_CaloId_TrackId_Iso_v*"),
        ## mu+el
        Mu23TrkIsoVVLEle12GsfCaloIdTrackIdIsoMediumWP = cms.vstring("HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v*"),
        Mu8TrkIsoVVLEle23GsfCaloIdTrackIdIsoMediumWP = cms.vstring("HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v*"),
        ## 8TeV
        IsoMu20eta2p1 = cms.vstring("HLT_IsoMu20_eta2p1_v*"),
        Ele27WP80 = cms.vstring("HLT_Ele27_WP80_v*"),
        ## 5.3E33 [50ns]
        IsoMu24eta2p1 = cms.vstring("HLT_IsoMu24_eta2p1_v*"),
        Ele27WP75eta2p1 = cms.vstring("HLT_Ele27_WP75_eta2p1_v*"),
        ## 7E33 [25ns] OOT+in time <mu> = 20
        IsoMu20eta2p1IterTrk02 = cms.vstring("HLT_IsoMu20_eta2p1_IterTrk02_v*"),
        ## 1.4E34 [25ns]
        IsoMu24eta2p1IterTrk02 = cms.vstring("HLT_IsoMu24_eta2p1_IterTrk02_v*"),
        Ele32WP75eta2p1 = cms.vstring("HLT_Ele32_WP75_eta2p1_v*"),
        Ele32WP75Gsf = cms.vstring("HLT_Ele32_WP75_Gsf_v*"),
    ),
)

