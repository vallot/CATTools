import FWCore.ParameterSet.Config as cms

catMuons = cms.EDProducer("CATMuonProducer",
    src = cms.InputTag("slimmedMuons"),
    mcLabel = cms.InputTag("prunedGenParticles"),
    vertexLabel = cms.InputTag("catVertex"),
    pfSrc  = cms.InputTag("packedPFCandidates"),
    beamLineSrc = cms.InputTag("offlineBeamSpot")
)
