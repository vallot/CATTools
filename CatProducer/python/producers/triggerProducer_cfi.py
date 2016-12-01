import FWCore.ParameterSet.Config as cms

catTrigger = cms.EDProducer("CATTriggerProducer",
    HLT = cms.PSet(
        triggerResults = cms.VInputTag(
            cms.InputTag("TriggerResults","","HLT2"),# due to reHLT, this is the first choice
            cms.InputTag("TriggerResults","","HLT"),# if above is not found, falls to default
        ),
        objects = cms.InputTag("selectedPatTrigger"),
        prescales = cms.InputTag("patTrigger"),
        select = cms.vstring(
            "HLT_Ele", "HLT_DoubleEle", 
            "HLT_Mu", "HLT_TkMu", "HLT_IsoMu", "HLT_IsoTkMu", 
            "HLT_DiMu", "HLT_DoubleIsoMu",
            "HLT_PFJet",
            "HLT_DoublePhoton", "HLT_Photon"
        ),
    ),

    Flag = cms.PSet(
        Flags = cms.VInputTag(
            ## Keep TriggerResults with prefix "Flag_"
            cms.InputTag("TriggerResults","","PAT"),
            cms.InputTag("TriggerResults","","RECO"),
        ),
        bools = cms.VInputTag(
            ## Keep EDProducer output with boolean type (Flag_ will be added in the front of the names)
            cms.InputTag("BadChargedCandidateFilter"),
            cms.InputTag("BadPFMuonFilter"),
        ),
        names = cms.vstring(
            ## Flags to keep
            "HBHENoiseFilter", "HBHENoiseIsoFilter",
            "BadPFMuonFilter", "BadChargedCandidateFilter",
            "EcalDeadCellTriggerPrimitiveFilter",
            "goodVertices",
            "eeBadScFilter",
            "globalTightHalo2016Filter",
        ),
    ),

)
