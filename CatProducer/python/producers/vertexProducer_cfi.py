import FWCore.ParameterSet.Config as cms

catVertex = cms.EDProducer("CATVertexProducer",
    vertexLabel = cms.InputTag("offlineSlimmedPrimaryVertices"),
    minimumNDOF = cms.uint32(4),
    maxAbsZ = cms.double(24),	
    maxd0 = cms.double(2)	
)
