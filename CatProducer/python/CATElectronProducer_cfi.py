import FWCore.ParameterSet.Config as cms

catElectrons = cms.EDProducer("CATElectronProducer",
    # input collection
    src = cms.InputTag("patElectrons"),
    vertexLabel = cms.InputTag('offlinePrimaryVertices')
    )

