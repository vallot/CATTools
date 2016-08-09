import FWCore.ParameterSet.Config as cms

catTrigger = cms.EDProducer("CATTriggerProducer",
    triggerBits = cms.VInputTag(
        cms.InputTag("TriggerResults","","HLT2"),# due to reHLT, this is the first choice
        cms.InputTag("TriggerResults","","HLT"),# if above is not found, falls to default
    ),
    metFilterBits = cms.VInputTag(
        cms.InputTag("TriggerResults","","PAT"),
        cms.InputTag("TriggerResults","","RECO"),
    ),
    triggerObjects = cms.InputTag("selectedPatTrigger"),
    triggerPrescales = cms.InputTag("patTrigger"),
    selectTrigObjects = cms.vstring("HLT_Ele", "HLT_DoubleEle", "HLT_Mu", "HLT_TkMu", "HLT_IsoMu", "HLT_IsoTkMu", "HLT_DoubleIsoMu", "HLT_PFJet","HLT_DoublePhoton", "HLT_Photon"),
    metFilterNames = cms.vstring(),
    hltPathNames = cms.vstring(),
##     metFilterNames = cms.vstring(
##     "HBHENoiseFilter",
##     "CSCTightHaloFilter",
##     "goodVertices",
##     "eeBadScFilter",
##     "EcalDeadCellTriggerPrimitiveFilter",
##     ),
##     hltPathNames = cms.vstring(
##     "HLT_Mu17_Mu8_DZ_v*","HLT_Mu17_TkMu8_DZ_v*","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*","HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*","HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*","HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*","HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet30_v*",
## "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v*","HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v*"
##     ),
)
