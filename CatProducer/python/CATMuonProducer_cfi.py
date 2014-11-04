import FWCore.ParameterSet.Config as cms

catMuonsSource = "selectedPatMuonsPFlow"
catVertexSource = "offlinePrimaryVertices"
catMCsource = "genParticles"
catBeamSpot = "offlineBeamSpot"

catMuons = cms.EDProducer("CATMuonProducer",
    # input collection
    src = cms.InputTag(catMuonsSource),
    mcLabel = cms.InputTag(catMCsource),
    vertexLabel = cms.InputTag(catVertexSource),
    beamLineSrc = cms.InputTag(catBeamSpot),
)

