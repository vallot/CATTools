import FWCore.ParameterSet.Config as cms

catVertexInfo = cms.EDProducer("CATVertexInfoProducer",
    vertexLabel = cms.InputTag("goodOfflinePrimaryVertices"),
)
