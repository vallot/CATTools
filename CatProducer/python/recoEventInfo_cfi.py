import FWCore.ParameterSet.Config as cms

recoEventInfo = cms.EDProducer("RecoEventInfoProducer",
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("selectedPatTrigger"),
    HLT = cms.PSet(
        ## from https://github.com/cms-sw/cmssw/blob/CMSSW_7_2_X/HLTriggerOffline/Top/python
        ## for PHYS14
        DoubleMu = cms.vstring('HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*','HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*','HLT_Mu17_TkMu8_v*'),
        DoubleEl = cms.vstring('HLT_Ele23_Ele12_CaloId_TrackId_iso_v*'),
        MuEl = cms.vstring('HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v*','HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v*'),
        SingleMu = cms.vstring('HLT_IsoMu24_IterTrk02_v*', 'HLT_IsoMu24_IterTrk02_TriCentralPFJet60_50_35_v*', 'HLT_IsoMu24_IterTrk02_TriCentralPFJet40_v*', 'HLT_IsoMu24_IterTrk02_CentralPFJet30_BTagCSV_v*', 'HLT_IsoMu20_eta2p1_IterTrk02_CentralPFJet30_BTagCSV_v*', 'HLT_IsoMu20_eta2p1_IterTrk02_TriCentralPFJet40_v*', 'HLT_IsoMu20_eta2p1_IterTrk02_TriCentralPFJet60_50_35_v*', 'HLT_IsoMu20_eta2p1_IterTrk02_v*', 'HLT_IsoTkMu20_eta2p1_IterTrk02_v*', 'HLT_IsoMu24_eta2p1_IterTrk02_CentralPFJet30_BTagCSV_v*', 'HLT_IsoMu24_eta2p1_IterTrk02_TriCentralPFJet40_v*', 'HLT_IsoMu24_eta2p1_IterTrk02_TriCentralPFJet60_50_35_v*', 'HLT_IsoMu24_eta2p1_IterTrk02_v*', 'HLT_IsoTkMu24_eta2p1_IterTrk02_v*'),
        SingleEl = cms.vstring('HLT_Ele27_eta2p1_WP85_Gsf_v*', 'HLT_Ele27_eta2p1_WP85_Gsf_TriCentralPFJet40_v*', 'HLT_Ele27_eta2p1_WP85_Gsf_TriCentralPFJet60_50_35_v*', 'HLT_Ele27_eta2p1_WP85_Gsf_CentralPFJet30_BTagCSV_v*', 'HLT_Ele32_eta2p1_WP85_Gsf_CentralPFJet30_BTagCSV_v*', 'HLT_Ele32_eta2p1_WP85_Gsf_TriCentralPFJet40_v*', 'HLT_Ele32_eta2p1_WP85_Gsf_TriCentralPFJet60_50_35_v*', 'HLT_Ele32_eta2p1_WP85_Gsf_v*')
    ),
)

