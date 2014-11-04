import FWCore.ParameterSet.Config as cms

catElectronsSource = "selectedPatElectronsPFlow"
catVertexSource = "offlinePrimaryVertices"

catElectrons = cms.EDProducer("CATElectronProducer",
    # input collection
    src = cms.InputTag(catElectronsSource),
    vertexLabel = cms.InputTag(catVertexSource)
)

