import FWCore.ParameterSet.Config as cms

catTrigger = cms.EDProducer("CATTriggerProducer",
    trigger = cms.PSet(
        triggerResults = cms.VInputTag(
            cms.InputTag("TriggerResults","","HLT2"),# due to reHLT, this is the first choice
            cms.InputTag("TriggerResults","","HLT"),# if above is not found, falls to default
        ),
        objects = cms.InputTag("selectedPatTrigger"),
        prescales = cms.InputTag("patTrigger"),
        prefix = cms.vstring(
            "HLT_Ele", "HLT_DoubleEle",
            "HLT_Mu", "HLT_TkMu", "HLT_IsoMu", "HLT_IsoTkMu",
            "HLT_DiMu", "HLT_DoubleIsoMu",
            "HLT_PFJet",
            "HLT_DoublePhoton", "HLT_Photon"
        ),
    ),
    flags = cms.PSet(
        triggerResults = cms.VInputTag(
            cms.InputTag("TriggerResults","","RECO"),
            cms.InputTag("TriggerResults","","PAT"),
        ),
        names = cms.vstring(
            "Flag_goodVertices",
            "Flag_globalTightHalo2016Filter",
            "Flag_HBHENoiseFilter",
            "Flag_HBHENoiseIsoFilter",
            "Flag_EcalDeadCellTriggerPrimitiveFilter",
            "Flag_eeBadScFilter",
        ),
        bools = cms.VInputTag(
            cms.InputTag("BadChargedCandidateFilter"),
            cms.InputTag("BadPFMuonFilter"),
        ),
    ),
)
