import FWCore.ParameterSet.Config as cms

catMuons = cms.EDProducer("CATMuonProducer",
    # input collection
    src = cms.InputTag("patMuons"),
    vertexLabel = cms.InputTag('offlinePrimaryVertices'),
)

