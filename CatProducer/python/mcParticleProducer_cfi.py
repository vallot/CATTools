import FWCore.ParameterSet.Config as cms

catMCParticles = cms.EDProducer("CATMCParticleProducer",
    # input collection
    src = cms.InputTag("genParticles"),
    pt = cms.double(10),
    eta = cms.double(2.5)
)

