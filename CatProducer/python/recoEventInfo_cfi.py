import FWCore.ParameterSet.Config as cms

recoEventInfo = cms.EDProducer("RecoEventInfoProducer",
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    HLT = cms.PSet(
        ## old triggers
        #Mu17Mu8 = cms.vstring("HLT_Mu17_Mu8_v*"),
        #Mu17TkMu8 = cms.vstring("HLT_Mu17_TkMu8_v*"),
        ## Phys14 sample
        DoubleMu = cms.vstring("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsonVVL_v*"),
        DoubleEl = cms.vstring("HLT_Ele23_Ele12_CaloId_TrackId_Iso_v*"),
        MuEl = cms.vstring("HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v*"),
        ElMu = cms.vstring("HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v*"),
        MuJet = cms.vstring("HLT_IsoMu24_eta2p1_v*"),
        ElJet = cms.vstring("HLT_Ele27_WP80_v*"),
        
        ## 5.3E33 [50ns] && 7E33[25ns]
        ##IsoMu20eta2p1 = cms.vstring("HLT_IsoMu20_eta2p1_IterTrk02_v*"),
        ##Ele27WP75eta2p1 = cms.vstring("HLT_Ele27_eta2p1_WP75_v*"),
        ## 1.4E34 [25ns]
        ##IsoMu24eta2p1IterTrk02 = cms.vstring("HLT_IsoMu24_eta2p1_IterTrk02_v*"),
        ##Ele32WP75eta2p1 = cms.vstring("HLT_Ele32_eta2p1_WP75_v*"),
    ),
)

