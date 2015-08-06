import FWCore.ParameterSet.Config as cms

catMuons = cms.EDProducer("CATMuonProducer",
    src = cms.InputTag("slimmedMuons"),
    shiftedEnDownSrc = cms.InputTag("shiftedSlimmedMuonsEnDown"),
    shiftedEnUpSrc = cms.InputTag("shiftedSlimmedMuonsEnUp"),
    mcLabel = cms.InputTag("prunedGenParticles"),
    vertexLabel = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamLineSrc = cms.InputTag("offlineBeamSpot")
)
