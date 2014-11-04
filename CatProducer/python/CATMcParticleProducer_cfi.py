import FWCore.ParameterSet.Config as cms

catMCsource = "genParticles"

catMCParticles = cms.EDProducer("CATMCParticleProducer",
    # input collection
    src = cms.InputTag(catMCsource),
    pt = cms.double(10),
    eta = cms.double(2.5)
)

