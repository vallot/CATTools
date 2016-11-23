import FWCore.ParameterSet.Config as cms

catMuons = cms.EDProducer("CATMuonProducer",
    src = cms.InputTag("slimmedMuons"),
    mcLabel = cms.InputTag("prunedGenParticles"),
    vertexLabel = cms.InputTag("catVertex"),
    beamLineSrc = cms.InputTag("offlineBeamSpot")
)
