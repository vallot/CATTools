import FWCore.ParameterSet.Config as cms

catMuons = cms.EDProducer("CATMuonProducer",
    # input collection
    src = cms.InputTag("patMuons"),
    mcLabel = cms.InputTag("genParticles"),
    vertexLabel = cms.InputTag('offlinePrimaryVertices'),
)

