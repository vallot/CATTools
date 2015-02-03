import FWCore.ParameterSet.Config as cms

catMuons = cms.EDProducer("CATMuonProducer",
    src = cms.InputTag("patMuons"),
    shiftedEnDownSrc = cms.InputTag("shiftedSlimmedMuonsEnDown"),
    shiftedEnUpSrc = cms.InputTag("shiftedSlimmedMuonsEnUp"),
    mcLabel = cms.InputTag("genParticles"),
    vertexLabel = cms.InputTag("offlinePrimaryVertices"),
    beamLineSrc = cms.InputTag("offlineBeamSpot")
)
