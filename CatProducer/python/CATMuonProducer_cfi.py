import FWCore.ParameterSet.Config as cms

catMuons = cms.EDProducer("CATMuonProducer",
    # input collection
    src = cms.InputTag("selectedPatMuonsPFlow"),
    mcLabel = cms.InputTag("genParticles"),
    vertexLabel = cms.InputTag("offlinePrimaryVertices"),
    beamLineSrc = cms.InputTag("offlineBeamSpot"),
    runOnMC = cms.bool(True)
)

